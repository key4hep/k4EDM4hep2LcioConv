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
#if __has_include("edm4hep/RawTimeSeriesCollection.h")
#include <edm4hep/RawTimeSeriesCollection.h>
#else
#include <edm4hep/TPCHitCollection.h>
namespace edm4hep {
  using RawTimeSeries = TPCHit;
  using MutableRawTimeSeries = MutableTPCHit;
  using RawTimeSeriesCollection = TPCHitCollection;
} // namespace edm4hep
#endif

#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackerHitCollection.h>
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
  using TypeMapT = k4EDM4hep2LcioConv::MapT<T1, T2>;

  struct CollectionsPairVectors {
    TypeMapT<lcio::TrackImpl*, edm4hep::Track> tracks {};
    TypeMapT<lcio::TrackerHitImpl*, edm4hep::TrackerHit> trackerhits {};
    TypeMapT<lcio::SimTrackerHitImpl*, edm4hep::SimTrackerHit> simtrackerhits {};
    TypeMapT<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit> calohits {};
    TypeMapT<lcio::RawCalorimeterHitImpl*, edm4hep::RawCalorimeterHit> rawcalohits {};
    TypeMapT<lcio::SimCalorimeterHitImpl*, edm4hep::SimCalorimeterHit> simcalohits {};
    TypeMapT<lcio::TPCHitImpl*, edm4hep::RawTimeSeries> tpchits {};
    TypeMapT<lcio::ClusterImpl*, edm4hep::Cluster> clusters {};
    TypeMapT<lcio::VertexImpl*, edm4hep::Vertex> vertices {};
    TypeMapT<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle> recoparticles {};
    TypeMapT<lcio::MCParticleImpl*, edm4hep::MCParticle> mcparticles {};
  };

  template<
    typename TrackMapT = TypeMapT<lcio::TrackImpl*, edm4hep::Track>,
    typename TrackerHitMapT = TypeMapT<lcio::TrackerHitImpl*, edm4hep::TrackerHit>>
  lcio::LCCollectionVec* convTracks(
    const edm4hep::TrackCollection* const tracks_coll,
    TrackMapT& tracks_vec,
    const TrackerHitMapT& trackerhits_vec);

  template<typename TrackerHitMapT = TypeMapT<lcio::TrackerHitImpl*, edm4hep::TrackerHit>>
  lcio::LCCollectionVec* convTrackerHits(
    const edm4hep::TrackerHitCollection* const trackerhits_coll,
    const std::string& cellIDstr,
    TrackerHitMapT& trackerhits_vec);

  template<
    typename SimTrHitMapT = TypeMapT<lcio::SimTrackerHitImpl*, edm4hep::SimTrackerHit>,
    typename MCParticleMapT = TypeMapT<lcio::MCParticleImpl*, edm4hep::MCParticle>>
  lcio::LCCollectionVec* convSimTrackerHits(
    const edm4hep::SimTrackerHitCollection* const simtrackerhits_coll,
    const std::string& cellIDstr,
    SimTrHitMapT& simtrackerhits_vec,
    const MCParticleMapT& mcparticles_vec);

  template<typename CaloHitMapT = TypeMapT<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit>>
  lcio::LCCollectionVec* convCalorimeterHits(
    const edm4hep::CalorimeterHitCollection* const calohit_coll,
    const std::string& cellIDstr,
    CaloHitMapT& calo_hits_vec);

  template<typename RawCaloHitMapT = TypeMapT<lcio::RawCalorimeterHitImpl*, edm4hep::RawCalorimeterHit>>
  lcio::LCCollectionVec* convRawCalorimeterHits(
    const edm4hep::RawCalorimeterHitCollection* const rawcalohit_coll,
    RawCaloHitMapT& raw_calo_hits_vec);

  template<typename SimCaloHitMapT, typename MCParticleMapT = TypeMapT<lcio::MCParticleImpl*, edm4hep::MCParticle>>
  lcio::LCCollectionVec* convSimCalorimeterHits(
    const edm4hep::SimCalorimeterHitCollection* const simcalohit_coll,
    const std::string& cellIDstr,
    SimCaloHitMapT& sim_calo_hits_vec,
    const MCParticleMapT& mcparticles);

  template<typename TPCHitMapT = TypeMapT<lcio::TPCHitImpl*, edm4hep::RawTimeSeries>>
  lcio::LCCollectionVec* convTPCHits(
    const edm4hep::RawTimeSeriesCollection* const tpchit_coll,
    TPCHitMapT& tpc_hits_vec);

  template<
    typename ClusterMapT = TypeMapT<lcio::ClusterImpl*, edm4hep::Cluster>,
    typename CaloHitMapT = TypeMapT<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit>>
  lcio::LCCollectionVec* convClusters(
    const edm4hep::ClusterCollection* const cluster_coll,
    ClusterMapT& cluster_vec,
    const CaloHitMapT& calohits_vec);

  template<
    typename VertexMapT = TypeMapT<lcio::VertexImpl*, edm4hep::Vertex>,
    typename RecoPartMapT = TypeMapT<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle>>
  lcio::LCCollectionVec* convVertices(
    const edm4hep::VertexCollection* const vertex_coll,
    VertexMapT& vertex_vec,
    const RecoPartMapT& recoparticles_vec);

  template<
    typename RecoPartMapT = TypeMapT<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle>,
    typename TrackMapT = TypeMapT<lcio::TrackImpl*, edm4hep::Track>,
    typename VertexMapT = TypeMapT<lcio::VertexImpl*, edm4hep::Vertex>,
    typename ClusterMapT = TypeMapT<lcio::ClusterImpl*, edm4hep::Cluster>>
  lcio::LCCollectionVec* convReconstructedParticles(
    const edm4hep::ReconstructedParticleCollection* const recos_coll,
    RecoPartMapT& recoparticles_vec,
    const TrackMapT& tracks_vec,
    const VertexMapT& vertex_vec,
    const ClusterMapT& clusters_vec);

  template<typename MCPartMapT = TypeMapT<lcio::MCParticleImpl*, edm4hep::MCParticle>>
  lcio::LCCollectionVec* convMCParticles(
    const edm4hep::MCParticleCollection* const mcparticle_coll,
    MCPartMapT& mc_particles_vec);

  void convEventHeader(const edm4hep::EventHeaderCollection* const header_coll, lcio::LCEventImpl* const lcio_event);

  template<typename ObjectMappingT = CollectionsPairVectors>
  void FillMissingCollections(ObjectMappingT& collection_pairs);

  bool collectionExist(const std::string& collection_name, const lcio::LCEventImpl* lcio_event);

  /**
   * Convert an edm4hep event to an LCEvent
   */
  std::unique_ptr<lcio::LCEventImpl> convEvent(
    const podio::Frame& edmEvent,
    const podio::Frame& metadata = podio::Frame {});

} // namespace EDM4hep2LCIOConv

#endif
