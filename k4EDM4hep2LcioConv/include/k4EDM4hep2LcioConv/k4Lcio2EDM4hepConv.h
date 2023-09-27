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
  using ObjectMapT = k4EDM4hep2LcioConv::MapT<LcioT, EdmT>;

  /**
   * Maping holding all the original and converted objects in a 1:1 mapping in a
   * way that makes the lookup from LCIO to EDM4hep easy.
   */
  struct LcioEdmTypeMapping {
    ObjectMapT<lcio::Track*, edm4hep::MutableTrack> tracks {};
    ObjectMapT<lcio::TrackerHit*, edm4hep::MutableTrackerHit> trackerHits {};
    ObjectMapT<lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit> simTrackerHits {};
    ObjectMapT<lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit> caloHits {};
    ObjectMapT<lcio::RawCalorimeterHit*, edm4hep::MutableRawCalorimeterHit> rawCaloHits {};
    ObjectMapT<lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit> simCaloHits {};
    ObjectMapT<lcio::TPCHit*, edm4hep::MutableRawTimeSeries> tpcHits {};
    ObjectMapT<lcio::Cluster*, edm4hep::MutableCluster> clusters {};
    ObjectMapT<lcio::Vertex*, edm4hep::MutableVertex> vertices {};
    ObjectMapT<lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle> recoParticles {};
    ObjectMapT<lcio::MCParticle*, edm4hep::MutableMCParticle> mcParticles {};
    ObjectMapT<lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane> trackerHitPlanes {};
    ObjectMapT<lcio::ParticleID*, edm4hep::MutableParticleID> particleIDs {};
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
  template<typename ObjectMappingT>
  std::vector<CollNamePair>
  convertCollection(const std::string& name, EVENT::LCCollection* LCCollection, ObjectMappingT& typeMapping);

  /**
   * Resolve all relations in all converted objects that are held in the map.
   * Dispatch to the correpsonding implementation for all the types that have
   * relations
   */
  template<typename ObjectMappingT>
  void resolveRelations(ObjectMappingT& typeMapping);

  template<typename ObjectMappingT, typename ObjectMappingU>
  void resolveRelations(ObjectMappingT& updateMaps, const ObjectMappingU& lookupMaps);

  /**
   * Convert LCRelation collections into the corresponding Association collections in EDM4hep
   */
  template<typename ObjectMappingT>
  std::vector<CollNamePair> createAssociations(
    const ObjectMappingT& typeMapping,
    const std::vector<std::pair<std::string, EVENT::LCCollection*>>& LCRelation);

  /**
   * Convert a subset collection, dispatching to the correct function for the
   * type of the input collection
   */
  template<typename ObjectMappingT>
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
  edm4hep::MutableParticleID convertParticleID(const EVENT::ParticleID* pid);

  /**
   * Convert an MCParticle collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename MCParticleMapT>
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
  template<typename RecoMapT, typename PIDMapT>
  std::vector<CollNamePair> convertReconstructedParticles(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    RecoMapT& recoparticlesMap,
    PIDMapT& particleIDMap);

  /**
   * Convert a Vertex collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename VertexMapT>
  std::unique_ptr<edm4hep::VertexCollection>
  convertVertices(const std::string& name, EVENT::LCCollection* LCCollection, VertexMapT& vertexMap);

  /**
   * Convert a SimTrackerHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename SimTrHitMapT>
  std::unique_ptr<edm4hep::SimTrackerHitCollection>
  convertSimTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, SimTrHitMapT& SimTrHitMap);

  /**
   * Convert a TPCHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT>
  std::unique_ptr<edm4hep::RawTimeSeriesCollection>
  convertTPCHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TPCHitMap);

  /**
   * Convert a TrackerHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT>
  std::unique_ptr<edm4hep::TrackerHitCollection>
  convertTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitMap);

  /**
   * Convert a TrackerHitPlane collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT>
  std::unique_ptr<edm4hep::TrackerHitPlaneCollection>
  convertTrackerHitPlanes(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitPlaneMap);

  /**
   * Convert a Track collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename TrackMapT>
  std::unique_ptr<edm4hep::TrackCollection>
  convertTracks(const std::string& name, EVENT::LCCollection* LCCollection, TrackMapT& TrackMap);

  /**
   * Convert a SimCalorimeterHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT>
  std::unique_ptr<edm4hep::SimCalorimeterHitCollection>
  convertSimCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& SimCaloHitMap);

  /**
   * Convert a RawCalorimeterHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT>
  std::unique_ptr<edm4hep::RawCalorimeterHitCollection>
  convertRawCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& rawCaloHitMap);

  /**
   * Convert a CalorimeterHit collection and return the resulting collection.
   * Simultaneously populates the mapping from LCIO to EDM4hep objects.
   */
  template<typename HitMapT>
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
  template<typename ClusterMapT, typename PIDMapT>
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
  auto handleSubsetColl(EVENT::LCCollection* lcioColl, const ObjectMapT& elemMap);

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
  createAssociationCollection(EVENT::LCCollection* relations, const FromMapT& fromMap, const ToMapT& toMap);

  /**
   * Creates the CaloHitContributions for all SimCaloHits.
   * has to be done this way, since the converted McParticles are needeed.
   * The contributions are also attached to their corresponding SimCalorimeterHits.
   */
  template<typename HitMapT, typename MCParticleMapT>
  std::unique_ptr<edm4hep::CaloHitContributionCollection> createCaloHitContributions(
    HitMapT& SimCaloHitMap,
    const MCParticleMapT& mcparticlesMap);

  /**
   * Resolve the relations for the MCParticles.
   */
  template<typename MCParticleMapT>
  void resolveRelationsMCParticles(MCParticleMapT& mcparticlesMap);

  /**
   * Resolve the relations for SimTrackerHits
   */
  template<typename HitMapT, typename MCParticleMapT>
  void resolveRelationsSimTrackerHits(HitMapT& SimTrHitMap, const MCParticleMapT& mcparticlesMap);

  /**
   * Resolve the relations for ReconstructedParticles
   */
  template<typename RecoParticleMapT, typename VertexMapT, typename ClusterMapT, typename TrackMapT>
  void resolveRelationsRecoParticles(
    RecoParticleMapT& recoparticlesMap,
    const VertexMapT& vertexMap,
    const ClusterMapT& clusterMap,
    const TrackMapT& tracksMap);

  /**
   * Resolve the relations for Clusters
   */
  template<typename ClusterMapT, typename CaloHitMapT>
  void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap);

  /**
   * Resolve the relations for Tracks
   */
  template<typename TrackMapT, typename TrackHitMapT, typename TPCHitMapT, typename THPlaneHitMapT>
  void resolveRelationsTracks(
    TrackMapT& tracksMap,
    const TrackHitMapT& trackerHitMap,
    const TPCHitMapT&,
    const THPlaneHitMapT&);

  /**
   * Resolve the relations for Vertices
   */
  template<typename VertexMapT, typename RecoParticleMapT>
  void resolveRelationsVertices(VertexMapT& vertexMap, const RecoParticleMapT& recoparticleMap);

} // namespace LCIO2EDM4hepConv

#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.ipp"

#endif // K4EDM4HEP2LCIOCONV_K4LCIO2EDM4HEPCONV_H
