#ifndef K4EDM4HEP2LCIOCONV_H
#define K4EDM4HEP2LCIOCONV_H

#include <cassert>

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
#include <edm4hep/TrackerHitCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/VertexCollection.h>

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

template<typename T1, typename T2>
using vec_pair = std::vector<std::pair<T1, T2>>;

struct CollectionsPairVectors {
  vec_pair<lcio::TrackImpl*, edm4hep::Track> tracks;
  vec_pair<lcio::TrackerHitImpl*, edm4hep::TrackerHit> trackerhits;
  vec_pair<lcio::SimTrackerHitImpl*, edm4hep::SimTrackerHit> simtrackerhits;
  vec_pair<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit> calohits;
  vec_pair<lcio::RawCalorimeterHitImpl*, edm4hep::RawCalorimeterHit> rawcalohits;
  vec_pair<lcio::SimCalorimeterHitImpl*, edm4hep::SimCalorimeterHit> simcalohits;
  vec_pair<lcio::TPCHitImpl*, edm4hep::RawTimeSeries> tpchits;
  vec_pair<lcio::ClusterImpl*, edm4hep::Cluster> clusters;
  vec_pair<lcio::VertexImpl*, edm4hep::Vertex> vertices;
  vec_pair<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle> recoparticles;
  vec_pair<lcio::MCParticleImpl*, edm4hep::MCParticle> mcparticles;
};

lcio::LCCollectionVec* convTracks(
  const edm4hep::TrackCollection* const tracks_coll,
  vec_pair<lcio::TrackImpl*, edm4hep::Track>& tracks_vec,
  const vec_pair<lcio::TrackerHitImpl*, edm4hep::TrackerHit>& trackerhits_vec);

lcio::LCCollectionVec* convTrackerHits(
  const edm4hep::TrackerHitCollection* const trackerhits_coll,
  const std::string cellIDstr,
  vec_pair<lcio::TrackerHitImpl*, edm4hep::TrackerHit>& trackerhits_vec);

lcio::LCCollectionVec* convSimTrackerHits(
  const edm4hep::SimTrackerHitCollection* const simtrackerhits_coll,
  const std::string cellIDstr,
  vec_pair<lcio::SimTrackerHitImpl*, edm4hep::SimTrackerHit>& simtrackerhits_vec,
  const vec_pair<lcio::MCParticleImpl*, edm4hep::MCParticle>& mcparticles_vec);

lcio::LCCollectionVec* convCalorimeterHits(
  const edm4hep::CalorimeterHitCollection* const calohit_coll,
  const std::string cellIDstr,
  vec_pair<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit>& calo_hits_vec);

lcio::LCCollectionVec* convRawCalorimeterHits(
  const edm4hep::RawCalorimeterHitCollection* const rawcalohit_coll,
  vec_pair<lcio::RawCalorimeterHitImpl*, edm4hep::RawCalorimeterHit>& raw_calo_hits_vec);

lcio::LCCollectionVec* convSimCalorimeterHits(
  const edm4hep::SimCalorimeterHitCollection* const simcalohit_coll,
  const std::string cellIDstr,
  vec_pair<lcio::SimCalorimeterHitImpl*, edm4hep::SimCalorimeterHit>& sim_calo_hits_vec,
  const vec_pair<lcio::MCParticleImpl*, edm4hep::MCParticle>& mcparticles);

lcio::LCCollectionVec* convTPCHits(
  const edm4hep::RawTimeSeriesCollection* const tpchit_coll,
  vec_pair<lcio::TPCHitImpl*, edm4hep::RawTimeSeries>& tpc_hits_vec);

lcio::LCCollectionVec* convClusters(
  const edm4hep::ClusterCollection* const cluster_coll,
  vec_pair<lcio::ClusterImpl*, edm4hep::Cluster>& cluster_vec,
  const vec_pair<lcio::CalorimeterHitImpl*, edm4hep::CalorimeterHit>& calohits_vec);

lcio::LCCollectionVec* convVertices(
  const edm4hep::VertexCollection* const vertex_coll,
  vec_pair<lcio::VertexImpl*, edm4hep::Vertex>& vertex_vec,
  const vec_pair<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle>& recoparticles_vec);

lcio::LCCollectionVec* convReconstructedParticles(
  const edm4hep::ReconstructedParticleCollection* const recos_coll,
  vec_pair<lcio::ReconstructedParticleImpl*, edm4hep::ReconstructedParticle>& recoparticles_vec,
  const vec_pair<lcio::TrackImpl*, edm4hep::Track>& tracks_vec,
  const vec_pair<lcio::VertexImpl*, edm4hep::Vertex>& vertex_vec,
  const vec_pair<lcio::ClusterImpl*, edm4hep::Cluster>& clusters_vec);

lcio::LCCollectionVec* convMCParticles(
  const edm4hep::MCParticleCollection* const mcparticle_coll,
  vec_pair<lcio::MCParticleImpl*, edm4hep::MCParticle>& mc_particles_vec);

void convEventHeader(const edm4hep::EventHeaderCollection* const header_coll, lcio::LCEventImpl* const lcio_event);

void FillMissingCollections(CollectionsPairVectors& collection_pairs);

bool collectionExist(const std::string& collection_name, const lcio::LCEventImpl* lcio_event);

#endif
