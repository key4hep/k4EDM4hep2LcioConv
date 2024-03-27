#include "ObjectMapping.h"

#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/TrackerHitPlane.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/TPCHit.h"
#include "EVENT/Cluster.h"
#include "EVENT/Vertex.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ParticleID.h"
#include "UTIL/LCIterator.h"
#include "UTIL/PIDHandler.h"

#include "edm4hep/TrackCollection.h"
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
  using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
} // namespace edm4hep
#endif
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/ParticleIDCollection.h"

#include "podio/Frame.h"

#include "k4EDM4hep2LcioConv/MappingUtils.h"

#include <vector>
#include <string>
#include <unordered_map>

template<typename LcioT, typename EDM4hepT, typename MapT>
void fillMap(MapT& map, EVENT::LCCollection* lcioColl, const EDM4hepT& edmColl)
{
  // Simply assume that the collections have the same size here
  UTIL::LCIterator<LcioT> lcioIt(lcioColl);

  for (const auto edmElem : edmColl) {
    const auto* lcioElem = lcioIt.next();
    map.emplace(lcioElem, edmElem.getObjectID());
  }
}

#define FILL_MAP(Type, mapName)                                  \
  if (type == #Type) {                                           \
    using namespace edm4hep;                                     \
    auto& edm4hepColl = edmEvt.get<Type::collection_type>(name); \
    fillMap<EVENT::Type>(mapName, lcioColl, edm4hepColl);        \
  }

void fillRecoPIDMaps(
  ObjectMappings::Map<const EVENT::ReconstructedParticle*>& recoMap,
  std::unordered_map<const EVENT::ParticleID*, edm4hep::ParticleID>& pidMap,
  const std::string& recoName,
  EVENT::LCEvent* lcEvt,
  const podio::Frame& edmEvt)
{
  UTIL::LCIterator<EVENT::ReconstructedParticle> lcioIt(lcEvt, recoName);
  UTIL::PIDHandler pidHandler(lcioIt());
  const auto& edmColl = edmEvt.get<edm4hep::ReconstructedParticleCollection>(recoName);

  // Simply need to assume here that per algorithm the pid "collections" run in
  // parallel.
  std::unordered_map<std::string, unsigned> algoIdcs {};
  for (const auto edmElem : edmColl) {
    const auto* lcioElem = lcioIt.next();
    recoMap.emplace(lcioElem, edmElem.getObjectID());

    for (const auto* lcioPid : lcioElem->getParticleIDs()) {
      const auto& algoName = pidHandler.getAlgorithmName(lcioPid->getAlgorithmType());
      auto idx = algoIdcs[algoName]++;
      // No need to do any fancy caching, because the Frame does that already
      const auto& pidColl = edmEvt.get<edm4hep::ParticleIDCollection>(recoName + "_PID_" + algoName);
      pidMap.emplace(lcioPid, pidColl[idx]);
    }
  }
}

std::string getRecoName(const std::string& pidName)
{
  const auto pos = pidName.find("_PID_");
  if (pos != std::string::npos) {
    return pidName.substr(0, pos);
  }
  return "";
}

ObjectMappings ObjectMappings::fromEvent(EVENT::LCEvent* lcEvt, const podio::Frame& edmEvt)
{
  ObjectMappings mapping {};
  for (const auto& name : *(lcEvt->getCollectionNames())) {
    // We only use non subset collections here, because we want the "real"
    // objects
    const auto lcioColl = lcEvt->getCollection(name);
    if (lcioColl->isSubset()) {
      continue;
    }
    const auto type = lcioColl->getTypeName();
    if (type == "LCRelation") {
      continue;
    }
    FILL_MAP(Track, mapping.tracks);
    FILL_MAP(TrackerHitPlane, mapping.trackerHitPlanes);
    FILL_MAP(SimTrackerHit, mapping.simTrackerHits);
    FILL_MAP(Cluster, mapping.clusters);
    FILL_MAP(CalorimeterHit, mapping.caloHits);
    FILL_MAP(RawCalorimeterHit, mapping.rawCaloHits);
    FILL_MAP(SimCalorimeterHit, mapping.simCaloHits);
    FILL_MAP(MCParticle, mapping.mcParticles);
    FILL_MAP(Vertex, mapping.vertices);
    // Need special treatment for TPCHit type mismatch
    if (type == "TPCHit") {
      auto& edm4hepColl = edmEvt.get<edm4hep::RawTimeSeriesCollection>(name);
      fillMap<EVENT::TPCHit>(mapping.tpcHits, lcioColl, edm4hepColl);
    }
    if (type == "TrackerHit") {
      auto& edm4hepColl = edmEvt.get<edm4hep::TrackerHit3DCollection>(name);
      fillMap<EVENT::TrackerHit>(mapping.trackerHits, lcioColl, edm4hepColl);
    }
    // Special treatment for reco particles and ParticleIDs because of
    // conceptual differences
    if (type == "ReconstructedParticle") {
      fillRecoPIDMaps(mapping.recoParticles, mapping.particleIDs, name, lcEvt, edmEvt);
    }
  }

  return mapping;
}
