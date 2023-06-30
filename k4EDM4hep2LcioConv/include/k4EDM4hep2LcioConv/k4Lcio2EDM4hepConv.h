#ifndef K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H
#define K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H

#include "k4EDM4hep2LcioConv/MappingUtils.h"

// EDM4hep
#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/MCRecoCaloParticleAssociationCollection.h"
#include "edm4hep/MCRecoClusterParticleAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/MCRecoTrackParticleAssociationCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoTrackerHitPlaneAssociationCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/RecoParticleVertexAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/VertexCollection.h"

// LCIO
#include <EVENT/CalorimeterHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ParticleID.h>
#include <EVENT/RawCalorimeterHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TPCHit.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <EVENT/Vertex.h>
#include <EVENT/LCIntVec.h>
#include <EVENT/LCFloatVec.h>
#include <UTIL/LCIterator.h>
#include <lcio.h>

#include "podio/Frame.h"
#include "podio/UserDataCollection.h"

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <tuple>
#include <vector>

namespace LCIO2EDM4hepConv {
  template<typename LcioT, typename EdmT>
  using TypeMapT = k4EDM4hep2LcioConv::MapT<LcioT, EdmT>;

  /**
   * Maping holding all the original and converted objects in a 1:1 mapping in a
   * way that makes the lookup from LCIO to EDM4hep easy.
   */
  struct LcioEdmTypeMapping {
    TypeMapT<lcio::Track*, edm4hep::MutableTrack> tracks {};
    TypeMapT<lcio::TrackerHit*, edm4hep::MutableTrackerHit> trackerHits {};
    TypeMapT<lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit> simTrackerHits {};
    TypeMapT<lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit> caloHits {};
    TypeMapT<lcio::RawCalorimeterHit*, edm4hep::MutableRawCalorimeterHit> rawCaloHits {};
    TypeMapT<lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit> simCaloHits {};
    TypeMapT<lcio::TPCHit*, edm4hep::MutableRawTimeSeries> tpcHits {};
    TypeMapT<lcio::Cluster*, edm4hep::MutableCluster> clusters {};
    TypeMapT<lcio::Vertex*, edm4hep::MutableVertex> vertices {};
    TypeMapT<lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle> recoParticles {};
    TypeMapT<lcio::MCParticle*, edm4hep::MutableMCParticle> mcParticles {};
    TypeMapT<lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane> trackerHitPlanes {};
    TypeMapT<lcio::ParticleID*, edm4hep::MutableParticleID> particleIDs {};
  };

  using CollNamePair = std::tuple<std::string, std::unique_ptr<podio::CollectionBase>>;

  /*
   * Convert a LCRunHeader to EDM4hep as a frame.
   */
  podio::Frame convertRunHeader(EVENT::LCRunHeader* rheader);

  /**
   * Convert a complete LCEvent from LCIO to EDM4hep.
   *
   * A second, optional argument can be passed to limit the collections to
   * convert to the subset that is passed. NOTE: There is an implicit assumption
   * here that collsToConvert only contains collection names that are present in
   * the passed evt. There is no exception handling internally to guard against
   * collections that are missing.
   */
  podio::Frame convertEvent(EVENT::LCEvent* evt, const std::vector<std::string>& collsToConvert = {});

  /**
   * Convert an LCIOCollection by dispatching to the specific conversion
   * function for the corresponding type (after querying the input collection).
   * Populates the correct object mapping along the way.
   *
   * Returns a vector of names and collections (since some LCIO collections will
   * result in more than one EDM4hep collection)
   */
  template<typename ObjectMappingT = LcioEdmTypeMapping>
  std::vector<CollNamePair>
  convertCollection(const std::string& name, EVENT::LCCollection* LCCollection, ObjectMappingT& typeMapping);

  /**
   * Resolve all relations in all converted objects that are held in the map.
   * Dispatch to the correpsonding implementation for all the types that have
   * relations
   */
  template<typename ObjectMappingT = LcioEdmTypeMapping>
  void resolveRelations(ObjectMappingT& typeMapping);

  /**
   * Convert LCRelation collections into the corresponding Association collections in EDM4hep
   */
  template<typename ObjectMappingT = LcioEdmTypeMapping>
  std::vector<CollNamePair> createAssociations(
    const ObjectMappingT& typeMapping,
    const std::vector<std::pair<std::string, EVENT::LCCollection*>>& LCRelation);

  /**
   * Convert a subset collection, dispatching to the correct function for the
   * type of the input collection
   */
  template<typename ObjectMappingT = LcioEdmTypeMapping>
  std::unique_ptr<podio::CollectionBase>
  fillSubset(EVENT::LCCollection* LCCollection, const ObjectMappingT& typeMapping, const std::string& type);

  /*
   * Converts a LCIntVec or LCFloatVec Collection into a podio::UserDataCollection of the appropriate type.
   *
   * NOTE: LC[Int|Float]Vec are nested, but podio::UserDataCollection are flat. Hence, this will put all
   * contents into one collection, and the [begin, end) indices in this collection into a second (flat)
   * collection (with the suffix "_VecLengths" added to its name), such that the elements at position i,
   * resp. (i + 1) form the [begin, end) indices for each of the original vector collections.
   */
  template<typename LCVecType>
  std::vector<CollNamePair> convertLCVec(const std::string& name, EVENT::LCCollection* LCCollection);

  /**
   * Converting all parameters of an LCIO Object and attaching them to the
   * passed podio::Frame.
   */
  template<typename LCIOType>
  void convertObjectParameters(LCIOType* lcioobj, podio::Frame& event);

  inline edm4hep::Vector3f Vector3fFrom(const double* v) { return edm4hep::Vector3f(v[0], v[1], v[2]); }

  inline edm4hep::Vector3f Vector3fFrom(const EVENT::FloatVec& v) { return edm4hep::Vector3f(v[0], v[1], v[2]); }

  /**
   * Convert a TrackState
   */
  edm4hep::TrackState convertTrackState(const EVENT::TrackState* trackState);

  /**
   * Convert a ParticleID object.
   *
   * In LCIO ParticleIDs are persisted as part of the ReconstructedParticle they
   * are attached to. There are no ParticleID collections in the event. Hence,
   * this function converts single particle ID objects and the management of
   * putting them into a collection and of creating the LCIO to EDM4hep mapping
   * is done in the conversion of the ReconstructedParticles.
   */
  edm4hep::MutableParticleID convertPaticleID(const EVENT::ParticleID* pid);

  /**
   * Convert an MCParticle collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename MCParticleMapT = TypeMapT<lcio::MCParticle*, edm4hep::MutableMCParticle>>
  std::unique_ptr<edm4hep::MCParticleCollection>
  convertMCParticles(const std::string& name, EVENT::LCCollection* LCCollection, MCParticleMapT& mcparticlesMap);

  /**
   * Convert a ReconstructedParticle collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   *
   * NOTE: Also populates a ParticleID collection, as those are persisted as
   * part of the ReconstructedParticles in LCIO. The name of this collection is
   * <name>_particleIDs
   */
  template<
    typename RecoMapT = TypeMapT<lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle>,
    typename PIDMapT = TypeMapT<lcio::ParticleID*, edm4hep::MutableParticleID>>
  std::vector<CollNamePair> convertReconstructedParticles(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    RecoMapT& recoparticlesMap,
    PIDMapT& particleIDMap);

  /**
   * Convert a Vertex collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename VertexMapT = TypeMapT<lcio::Vertex*, edm4hep::MutableVertex>>
  std::unique_ptr<edm4hep::VertexCollection>
  convertVertices(const std::string& name, EVENT::LCCollection* LCCollection, VertexMapT& vertexMap);

  /**
   * Convert a SimTrackerHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename SimTrHitMapT = TypeMapT<lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit>>
  std::unique_ptr<edm4hep::SimTrackerHitCollection>
  convertSimTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, SimTrHitMapT& SimTrHitMap);

  /**
   * Convert a TPCHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT = TypeMapT<lcio::TPCHit*, edm4hep::MutableRawTimeSeries>>
  std::unique_ptr<edm4hep::RawTimeSeriesCollection>
  convertTPCHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TPCHitMap);

  /**
   * Convert a TrackerHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT = TypeMapT<lcio::TrackerHit*, edm4hep::MutableTrackerHit>>
  std::unique_ptr<edm4hep::TrackerHitCollection>
  convertTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitMap);

  /**
   * Convert a TrackerHitPlane collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT = TypeMapT<lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane>>
  std::unique_ptr<edm4hep::TrackerHitPlaneCollection>
  convertTrackerHitPlanes(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitPlaneMap);

  /**
   * Convert a Track collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename TrackMapT = TypeMapT<lcio::Track*, edm4hep::MutableTrack>>
  std::unique_ptr<edm4hep::TrackCollection>
  convertTracks(const std::string& name, EVENT::LCCollection* LCCollection, TrackMapT& TrackMap);

  /**
   * Convert a SimCalorimeterHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT = TypeMapT<lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit>>
  std::unique_ptr<edm4hep::SimCalorimeterHitCollection>
  convertSimCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& SimCaloHitMap);

  /**
   * Convert a RawCalorimeterHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT = TypeMapT<lcio::RawCalorimeterHit*, edm4hep::MutableRawCalorimeterHit>>
  std::unique_ptr<edm4hep::RawCalorimeterHitCollection>
  convertRawCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& rawCaloHitMap);

  /**
   * Convert a CalorimeterHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT = TypeMapT<lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit>>
  std::unique_ptr<edm4hep::CalorimeterHitCollection>
  convertCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& caloHitMap);

  /**
   * Convert a Cluster collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   *
   * NOTE: Also populates a ParticleID collection, as those are persisted as
   * part of the Cluster collection in LCIO. The name of this collection is
   * <name>_particleIDs
   */
  template<
    typename ClusterMapT = TypeMapT<lcio::Cluster*, edm4hep::MutableCluster>,
    typename PIDMapT = TypeMapT<lcio::ParticleID*, edm4hep::MutableParticleID>>
  std::vector<CollNamePair> convertClusters(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    ClusterMapT& clusterMap,
    PIDMapT& particleIDMap);

  /**
   * Create an EventHeaderCollection and fills it with the Metadata.
   */

  std::unique_ptr<edm4hep::EventHeaderCollection> createEventHeader(const EVENT::LCEvent* evt);

  /**
   * Helper function to create a subset collection from an existing (LCIO)
   * subset collection. Needs the object mapping as input to resolve to the
   * correct EDM4hep object for a given LCIO object.
   *
   * NOTE: Users responsibility to call this with the right inputs (i.e.
   * matching types)
   */
  template<
    typename CollT,
    typename ObjectMapT,
    typename LcioT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<ObjectMapT>>,
    typename Edm4hepT = k4EDM4hep2LcioConv::detail::mapped_t<ObjectMapT>>
  auto handleSubsetColl(EVENT::LCCollection* lcioColl, const ObjectMapT& elemMap)
  {
    auto edm4hepColl = std::make_unique<CollT>();
    edm4hepColl->setSubsetCollection();

    UTIL::LCIterator<LcioT> lcioIter(lcioColl);
    while (const auto lcioElem = lcioIter.next()) {
      if (auto edm4hepElem = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioElem, elemMap)) {
        edm4hepColl->push_back(edm4hepElem.value());
      }
      else {
        std::cerr << "Cannot find corresponding EDM4hep object for an LCIO object in a subset collection" << std::endl;
      }
    }

    return edm4hepColl;
  }

  /**
   * Create an Association collection from an LCRelations collection. Templated
   * on the From and To types as well as the direction of the relations in the
   * input LCRelations collection with respect to the order in which they are
   * mentioned in the Association collection of EDM4hep (since those are not
   * directed).
   *
   * Necessary inputs apart from the LCRelations collection are the correct LCIO
   * to EDM4hep object mappings to actually resolve the necessary relations.
   */
  template<
    typename CollT,
    bool Reverse,
    typename FromMapT,
    typename ToMapT,
    typename FromLCIOT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<FromMapT>>,
    typename ToLCIOT = std::remove_pointer_t<k4EDM4hep2LcioConv::detail::key_t<ToMapT>>,
    typename FromEDM4hepT = k4EDM4hep2LcioConv::detail::mapped_t<FromMapT>,
    typename ToEDM4hepT = k4EDM4hep2LcioConv::detail::mapped_t<ToMapT>>
  std::unique_ptr<CollT>
  createAssociationCollection(EVENT::LCCollection* relations, const FromMapT& fromMap, const ToMapT& toMap)
  {
    auto assocColl = std::make_unique<CollT>();
    auto relIter = UTIL::LCIterator<EVENT::LCRelation>(relations);

    while (const auto rel = relIter.next()) {
      auto assoc = assocColl->create();
      assoc.setWeight(rel->getWeight());
      const auto lcioTo = static_cast<ToLCIOT*>(rel->getTo());
      const auto lcioFrom = static_cast<FromLCIOT*>(rel->getFrom());
      const auto edm4hepTo = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioTo, toMap);
      const auto edm4hepFrom = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioFrom, fromMap);
      if (edm4hepTo.has_value() && edm4hepFrom.has_value()) {
        if constexpr (Reverse) {
          if constexpr (std::is_same_v<ToEDM4hepT, edm4hep::MutableVertex>) {
            assoc.setVertex(*edm4hepTo);
          }
          else {
            assoc.setSim(*edm4hepTo);
          }
          assoc.setRec(*edm4hepFrom);
        }
        else {
          if constexpr (std::is_same_v<FromEDM4hepT, edm4hep::MutableVertex>) {
            assoc.setVertex(*edm4hepFrom);
          }
          else {
            assoc.setSim(*edm4hepFrom);
          }
          assoc.setRec(*edm4hepTo);
        }
      }
    }

    return assocColl;
  }

  /**
   * Creates the CaloHitContributions for all SimCaloHits.
   * has to be done this way, since the converted McParticles are needeed.
   * The contributions are also attached to their corresponding SimCalorimeterHits.
   */
  template<
    typename HitMapT = TypeMapT<lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit>,
    typename MCParticleMapT = TypeMapT<lcio::MCParticle*, edm4hep::MutableMCParticle>>
  std::unique_ptr<edm4hep::CaloHitContributionCollection> createCaloHitContributions(
    HitMapT& SimCaloHitMap,
    const MCParticleMapT& mcparticlesMap);

  /**
   * Resolve the relations for the MCParticles.
   */
  template<typename MCParticleMapT = TypeMapT<lcio::MCParticle*, edm4hep::MutableMCParticle>>
  void resolveRelationsMCParticles(MCParticleMapT& mcparticlesMap);

  /**
   * Resolve the relations for SimTrackerHits
   */
  template<
    typename HitMapT = TypeMapT<lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit>,
    typename MCParticleMapT = TypeMapT<lcio::MCParticle*, edm4hep::MutableMCParticle>>
  void resolveRelationsSimTrackerHits(HitMapT& SimTrHitMap, const MCParticleMapT& mcparticlesMap);

  /**
   * Resolve the relations for ReconstructedParticles
   */
  template<
    typename RecoParticleMapT = TypeMapT<lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle>,
    typename VertexMapT = TypeMapT<lcio::Vertex*, edm4hep::MutableVertex>,
    typename ClusterMapT = TypeMapT<lcio::Cluster*, edm4hep::MutableCluster>,
    typename TrackMapT = TypeMapT<lcio::Track*, edm4hep::MutableTrack>>
  void resolveRelationsRecoParticles(
    RecoParticleMapT& recoparticlesMap,
    const VertexMapT& vertexMap,
    const ClusterMapT& clusterMap,
    const TrackMapT& tracksMap);

  /**
   * Resolve the relations for Clusters
   */
  template<
    typename ClusterMapT = TypeMapT<lcio::Cluster*, edm4hep::MutableCluster>,
    typename CaloHitMapT = TypeMapT<lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit>>
  void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap);

  /**
   * Resolve the relations for Tracks
   */
  template<
    typename TrackMapT = TypeMapT<lcio::Track*, edm4hep::MutableTrack>,
    typename TrackHitMapT = TypeMapT<lcio::TrackerHit*, edm4hep::MutableTrackerHit>,
    typename TPCHitMapT = TypeMapT<lcio::TPCHit*, edm4hep::MutableRawTimeSeries>,
    typename THPlaneHitMapT = TypeMapT<lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane>>
  void resolveRelationsTracks(
    TrackMapT& tracksMap,
    const TrackHitMapT& trackerHitMap,
    const TPCHitMapT&,
    const THPlaneHitMapT&);

  /**
   * Resolve the relations for Vertices
   */
  template<
    typename VertexMapT = TypeMapT<lcio::Vertex*, edm4hep::MutableVertex>,
    typename RecoParticleMapT = TypeMapT<lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle>>
  void resolveRelationsVertices(VertexMapT& vertexMap, const RecoParticleMapT& recoparticleMap);

  template<typename LCIOType>
  void convertObjectParameters(LCIOType* lcioobj, podio::Frame& event)
  {
    const auto& params = lcioobj->getParameters();
    // handle srting params
    EVENT::StringVec keys;
    const auto stringKeys = params.getStringKeys(keys);
    for (int i = 0; i < stringKeys.size(); i++) {
      EVENT::StringVec sValues;
      const auto stringVals = params.getStringVals(stringKeys[i], sValues);
      event.putParameter(stringKeys[i], stringVals);
    }
    // handle float params
    EVENT::StringVec fkeys;
    const auto floatKeys = params.getFloatKeys(fkeys);
    for (int i = 0; i < floatKeys.size(); i++) {
      EVENT::FloatVec fValues;
      const auto floatVals = params.getFloatVals(floatKeys[i], fValues);
      event.putParameter(floatKeys[i], floatVals);
    }
    // handle int params
    EVENT::StringVec ikeys;
    const auto intKeys = params.getIntKeys(ikeys);
    for (int i = 0; i < intKeys.size(); i++) {
      EVENT::IntVec iValues;
      const auto intVals = params.getIntVals(intKeys[i], iValues);
      event.putParameter(intKeys[i], intVals);
    }
    // handle double params
    EVENT::StringVec dkeys;
    const auto dKeys = params.getDoubleKeys(dkeys);
    for (int i = 0; i < dKeys.size(); i++) {
      EVENT::DoubleVec dValues;
      const auto dVals = params.getDoubleVals(dKeys[i], dValues);
      event.putParameter(dKeys[i], dVals);
    }
  }

  template<typename LCVecType>
  std::vector<CollNamePair> convertLCVec(const std::string& name, EVENT::LCCollection* LCCollection)
  {
    auto dest = std::make_unique<podio::UserDataCollection<typename LCVecType::value_type>>();
    auto vecSizes = std::make_unique<podio::UserDataCollection<uint32_t>>();
    if (LCCollection->getNumberOfElements() > 0) {
      vecSizes->push_back(0);
    }
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<LCVecType*>(LCCollection->getElementAt(i));
      for (unsigned j = 0; j < rval->size(); j++) {
        dest->push_back((*rval)[j]);
      }
      vecSizes->push_back(dest->size());
    }
    std::vector<CollNamePair> results;
    results.reserve(2);
    results.emplace_back(name, std::move(dest));
    results.emplace_back(name + "_VecLenghts", std::move(vecSizes));
    return results;
  }

} // namespace LCIO2EDM4hepConv

#endif // K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H
