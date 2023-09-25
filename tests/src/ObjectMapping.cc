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
#include "EVENT/ParticleID.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "UTIL/LCIterator.h"

#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/VertexCollection.h"

#include "podio/ObjectID.h"
#include "podio/Frame.h"

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
    FILL_MAP(TrackerHit, mapping.trackerHits);
    FILL_MAP(TrackerHitPlane, mapping.trackerHitPlanes);
    FILL_MAP(SimTrackerHit, mapping.simTrackerHits);
    FILL_MAP(Cluster, mapping.clusters);
    FILL_MAP(CalorimeterHit, mapping.caloHits);
    FILL_MAP(RawCalorimeterHit, mapping.rawCaloHits);
    FILL_MAP(SimCalorimeterHit, mapping.simCaloHits);
    FILL_MAP(MCParticle, mapping.mcParticles);
    FILL_MAP(ReconstructedParticle, mapping.recoParticles);
    FILL_MAP(Vertex, mapping.vertices);
    // Need special treatment for TPCHit type mismatch
    if (type == "TPCHit") {
      auto& edm4hepColl = edmEvt.get<edm4hep::RawTimeSeriesCollection>(name);
      fillMap<EVENT::TPCHit>(mapping.tpcHits, lcioColl, edm4hepColl);
    }
  }

  return mapping;
}
