#ifndef K4EDM4HEP2LCIOCONV_H
#define K4EDM4HEP2LCIOCONV_H

#include "k4EDM4hep2LcioConv/MappingUtils.h"

// EDM4hep
#include <edm4hep/CaloHitContributionCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/ClusterCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/MCRecoCaloAssociationCollection.h>
#include <edm4hep/MCRecoCaloParticleAssociationCollection.h>
#include <edm4hep/MCRecoParticleAssociationCollection.h>
#include <edm4hep/MCRecoTrackParticleAssociationCollection.h>
#include <edm4hep/MCRecoTrackerAssociationCollection.h>
#include <edm4hep/ParticleIDCollection.h>
#include <edm4hep/RawCalorimeterHitCollection.h>
#include <edm4hep/RecoParticleVertexAssociationCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/RawTimeSeriesCollection.h>
#include <edm4hep/TrackCollection.h>
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
  using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
  using TrackerHit3D = edm4hep::TrackerHit;
} // namespace edm4hep
#endif
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/VertexCollection.h>

#include "podio/Frame.h"

// LCIO
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/ClusterImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/RawCalorimeterHitImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/SimCalorimeterHitImpl.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <IMPL/TPCHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/VertexImpl.h>
#include <LCIOSTLTypes.h>
#include <UTIL/CellIDEncoder.h>
#include <lcio.h>

#include <memory>

// Preprocessor symbol that can be used in downstream code to switch on the
// namespace for the conversion
#define EDM4HEP2LCIOCONV_NAMESPACE 1

namespace EDM4hep2LCIOConv {

  template<typename T1, typename T2>
  using ObjectMapT = k4EDM4hep2LcioConv::VecMapT<T1, T2>;

  template<typename T1, typename T2>
  using vec_pair [[deprecated("Use a more descriptive alias")]] = ObjectMapT<T1, T2>;

  struct CollectionsPairVectors {
    ObjectMapT<lcio::TrackImpl*, edm4hep::Track> tracks {};
    ObjectMapT<lcio::TrackerHitImpl*, edm4hep::TrackerHit3D> trackerHits {};
    ObjectMapT<lcio::TrackerHitPlaneImpl*, edm4hep::TrackerHitPlane> trackerHitPlanes {};
    ObjectMapT<lcio::SimTrackerHitImpl*, edm4hep::SimTrackerHit> simTrackerHits {};
    ObjectMapT<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit> caloHits {};
    ObjectMapT<lcio::RawCalorimeterHitImpl*, edm4hep::RawCalorimeterHit> rawCaloHits {};
    ObjectMapT<lcio::SimCalorimeterHitImpl*, edm4hep::SimCalorimeterHit> simCaloHits {};
    ObjectMapT<lcio::TPCHitImpl*, edm4hep::RawTimeSeries> tpcHits {};
    ObjectMapT<lcio::ClusterImpl*, edm4hep::Cluster> clusters {};
    ObjectMapT<lcio::VertexImpl*, edm4hep::Vertex> vertices {};
    ObjectMapT<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle> recoParticles {};
    ObjectMapT<lcio::MCParticleImpl*, edm4hep::MCParticle> mcParticles {};
  };

  /**
   * Convert EDM4hep Tracks to LCIO. Simultaneously populate the mapping from
   * EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename TrackMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertTracks(
    const edm4hep::TrackCollection* const edmCollection,
    TrackMapT& trackMap);

  template<typename TrackMapT, typename TrackerHitMapT>
  [[deprecated("Use convertTracks instead")]] lcio::LCCollectionVec*
  convTracks(const edm4hep::TrackCollection* const tracks_coll, TrackMapT& tracks_vec, const TrackerHitMapT&)
  {
    return convertTracks(tracks_coll, tracks_vec).release();
  }

  /**
   * Convert EDM4hep TrackerHits to LCIO. Simultaneously populate mapping from
   * EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename TrackerHitMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertTrackerHits(
    const edm4hep::TrackerHit3DCollection* const edmCollection,
    const std::string& cellIDStr,
    TrackerHitMapT& trackerHitMap);

  template<typename TrackerHitMapT>
  [[deprecated("Use convertTrackerHits instead")]] lcio::LCCollectionVec* convTrackerHits(
    const edm4hep::TrackerHit3DCollection* const trackerhits_coll,
    const std::string& cellIDstr,
    TrackerHitMapT& trackerhits_vec)
  {
    return convertTrackerHits(trackerhits_coll, cellIDstr, trackerhits_vec).release();
  }

  /**
   * Convert EDM4hep TrackerHitPlanes to LCIO. Simultaneously populate mapping
   * from EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename TrackerHitPlaneMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertTrackerHitPlanes(
    const edm4hep::TrackerHitPlaneCollection* const edmCollection,
    const std::string& cellIDstr,
    TrackerHitPlaneMapT& trackerHitsMap);

  template<typename TrackerHitPlaneMapT>
  [[deprecated("Use convertTrackerHitPlanes instead")]] lcio::LCCollectionVec* convTrackerHitPlanes(
    const edm4hep::TrackerHitPlaneCollection* const trackerhits_coll,
    const std::string& cellIDstr,
    TrackerHitPlaneMapT& trackerhits_vec)
  {
    return convertTrackerHitPlanes(trackerhits_coll, cellIDstr, trackerhits_vec).release();
  }

  /**
   * Convert EDM4hep SimTrackerHits to LCIO. Simultaneously populate mapping
   * from EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename SimTrHitMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertSimTrackerHits(
    const edm4hep::SimTrackerHitCollection* const edmCollection,
    const std::string& cellIDstr,
    SimTrHitMapT& simTrHitMap);

  template<typename SimTrHitMapT, typename MCParticleMapT>
  [[deprecated("Use convertSimTrackerHits instead")]] lcio::LCCollectionVec* convSimTrackerHits(
    const edm4hep::SimTrackerHitCollection* const simtrackerhits_coll,
    const std::string& cellIDstr,
    SimTrHitMapT& simtrackerhits_vec,
    const MCParticleMapT&)
  {
    return convertSimTrackerHits(simtrackerhits_coll, cellIDstr, simtrackerhits_vec).release();
  }

  /**
   * Convert EDM4hep CalorimeterHits to LCIO. Simultaneously populate mapping
   * from EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename CaloHitMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertCalorimeterHits(
    const edm4hep::CalorimeterHitCollection* const edmCollection,
    const std::string& cellIDstr,
    CaloHitMapT& caloHitMap);

  template<typename CaloHitMapT>
  [[deprecated("Use convertCalorimeterHits instead")]] lcio::LCCollectionVec* convCalorimeterHits(
    const edm4hep::CalorimeterHitCollection* const calohit_coll,
    const std::string& cellIDstr,
    CaloHitMapT& calo_hits_vec)
  {
    return convertCalorimeterHits(calohit_coll, cellIDstr, calo_hits_vec).release();
  }

  /**
   * Convert EDM4hep RawCalorimeterHits to LCIO. Simultaneously populate mapping
   * from EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename RawCaloHitMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertRawCalorimeterHits(
    const edm4hep::RawCalorimeterHitCollection* const edmCollection,
    RawCaloHitMapT& rawCaloHitMap);

  template<typename RawCaloHitMapT>
  [[deprecated("use convertRawCalorimeterHits instead")]] lcio::LCCollectionVec* convRawCalorimeterHits(
    const edm4hep::RawCalorimeterHitCollection* const rawcalohit_coll,
    RawCaloHitMapT& raw_calo_hits_vec)
  {
    return convertRawCalorimeterHits(rawcalohit_coll, raw_calo_hits_vec).release();
  }

  /**
   * Convert EDM4hep SimCalorimeterHits to LCIO. Simultaneously populate mapping
   * from EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename SimCaloHitMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertSimCalorimeterHits(
    const edm4hep::SimCalorimeterHitCollection* const edmCollection,
    const std::string& cellIDstr,
    SimCaloHitMapT& simCaloHitMap);

  template<typename SimCaloHitMapT>
  lcio::LCCollectionVec* convSimCalorimeterHits(
    const edm4hep::SimCalorimeterHitCollection* const simcalohit_coll,
    const std::string& cellIDstr,
    SimCaloHitMapT& sim_calo_hits_vec)
  {
    return convertSimCalorimeterHits(simcalohit_coll, cellIDstr, sim_calo_hits_vec).release();
  }

  template<typename SimCaloHitMapT, typename MCParticleMapT>
  [[deprecated("remove MCParticleMap argument since it is unused")]] lcio::LCCollectionVec* convSimCalorimeterHits(
    const edm4hep::SimCalorimeterHitCollection* const simcalohit_coll,
    const std::string& cellIDstr,
    SimCaloHitMapT& sim_calo_hits_vec,
    const MCParticleMapT&)
  {
    return convSimCalorimeterHits(simcalohit_coll, cellIDstr, sim_calo_hits_vec);
  }

  /**
   * Convert EDM4hep TPC Hits to LCIO. Simultaneously populate mapping from
   * EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename TPCHitMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertTPCHits(
    const edm4hep::RawTimeSeriesCollection* const edmCollection,
    TPCHitMapT& tpcHitMap);

  template<typename TPCHitMapT>
  [[deprecated("Use convertTPCHits instead")]] lcio::LCCollectionVec* convTPCHits(
    const edm4hep::RawTimeSeriesCollection* const tpchit_coll,
    TPCHitMapT& tpc_hits_vec)
  {
    return convertTPCHits(tpchit_coll, tpc_hits_vec).release();
  }

  /**
   * Convert EDM4hep Clusters to LCIO. Simultaneously populate mapping from
   * EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename ClusterMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertClusters(
    const edm4hep::ClusterCollection* const edmCollection,
    ClusterMapT& clusterMap);

  template<typename ClusterMapT>
  [[deprecated("Use convertClusters instead")]] lcio::LCCollectionVec* convClusters(
    const edm4hep::ClusterCollection* const cluster_coll,
    ClusterMapT& cluster_vec)
  {
    return convertClusters(cluster_coll, cluster_vec).release();
  }

  template<typename ClusterMapT, typename CaloHitMapT>
  [[deprecated("remove CaloHitMap argument since it is unused")]] lcio::LCCollectionVec*
  convClusters(const edm4hep::ClusterCollection* const cluster_coll, ClusterMapT& cluster_vec, const CaloHitMapT&)
  {
    return convClusters(cluster_coll, cluster_vec);
  }

  /**
   * Convert EDM4hep Vertices to LCIO. Simultaneously populate mapping from
   * EDM4hep to LCIO objects for relation resolving in a second step.
   */
  template<typename VertexMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertVertices(
    const edm4hep::VertexCollection* const edmCollection,
    VertexMapT& vertexMap);

  template<typename VertexMapT, typename RecoPartMapT>
  [[deprecated("Use convertVertices instead")]] lcio::LCCollectionVec*
  convVertices(const edm4hep::VertexCollection* const vertex_coll, VertexMapT& vertex_vec, const RecoPartMapT&)
  {
    return convertVertices(vertex_coll, vertex_vec).release();
  }

  /**
   * Convert EDM4hep ReconstructedParticles to LCIO. Simultaneously populate
   * mapping from EDM4hep to LCIO objects for relation resolving in a second
   * step.
   */
  template<typename RecoParticleMapT>
  std::unique_ptr<lcio::LCCollectionVec> convertReconstructedParticles(
    const edm4hep::ReconstructedParticleCollection* const edmCollection,
    RecoParticleMapT& recoParticleMap);

  template<typename RecoPartMapT, typename TrackMapT, typename VertexMapT, typename ClusterMapT>
  [[deprecated("Use convertReconstructedParticle instead")]] lcio::LCCollectionVec* convReconstructedParticles(
    const edm4hep::ReconstructedParticleCollection* const recos_coll,
    RecoPartMapT& recoparticles_vec,
    const TrackMapT&,
    const VertexMapT&,
    const ClusterMapT&)
  {
    return convertReconstructedParticles(recos_coll, recoparticles_vec).release();
  }

  template<typename MCPartMapT>
  lcio::LCCollectionVec* convMCParticles(
    const edm4hep::MCParticleCollection* const mcparticle_coll,
    MCPartMapT& mc_particles_vec);


  /**
   * Convert EDM4hep EventHeader to LCIO. This will directly populate the
   * corresponding information in the passed LCEvent. The input collection needs
   * to be of length 1!
   */
  void convertEventHeader(const edm4hep::EventHeaderCollection* const header_coll, lcio::LCEventImpl* const lcio_event);

  [[deprecated("Use convertEventheader instead")]] void convEventHeader(
    const edm4hep::EventHeaderCollection* const header_coll,
    lcio::LCEventImpl* const lcio_event);

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
   * Resolve the relations for SimTrackerHits
   */
  template<typename SimTrHitMapT, typename MCParticleMapT>
  void resolveRelationsSimTrackerHits(SimTrHitMapT& simTrHitMap, const MCParticleMapT& mcParticleMap);

  /**
   * Resolve the relations for Vertex
   */
  template<typename VertexMapT, typename RecoParticleMapT>
  void resolveRelationsVertices(VertexMapT& vertexMap, const RecoParticleMapT& recoParticleMap);

  /**
   * Resolve the relations for SimCalorimeterHit. This is also the step where
   * the LCIO SimCalorimeterHits get their contributions attached.
   */
  template<typename SimCaloHitMapT, typename MCParticleMapT>
  void resolveRelationsSimCaloHit(SimCaloHitMapT& simCaloHitMap, const MCParticleMapT& mcParticleMap);

  /**
   * Resolve the relations for ReconstructedParticles
   */
  template<
    typename RecoParticleMapT,
    typename RecoParticleLookupMapT,
    typename VertexMapT,
    typename ClusterMapT,
    typename TrackMapT>
  void resolveRelationsRecoParticles(
    RecoParticleMapT& recoParticleMap,
    const RecoParticleLookupMapT& recoLookupMap,
    const VertexMapT& vertexMap,
    const ClusterMapT& clusterMap,
    const TrackMapT& trackMap);

  /**
   * Resolve the relations for Clusters
   */
  template<typename ClusterMapT, typename CaloHitMapT>
  void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap);

  template<typename ObjectMappingT>
  void FillMissingCollections(ObjectMappingT& update_pairs);

  /// Update the relations of the objects in the update_pairs map, by linking
  /// them according to the contents of the lookup_pairs map
  template<typename ObjectMappingT, typename ObjectMappingU>
  void FillMissingCollections(ObjectMappingT& update_pairs, const ObjectMappingU& lookup_pairs);

  bool collectionExist(const std::string& collection_name, const lcio::LCEventImpl* lcio_event);

  /**
   * Convert an edm4hep event to an LCEvent
   */
  std::unique_ptr<lcio::LCEventImpl> convEvent(
    const podio::Frame& edmEvent,
    const podio::Frame& metadata = podio::Frame {});

} // namespace EDM4hep2LCIOConv

#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.ipp"

#endif
