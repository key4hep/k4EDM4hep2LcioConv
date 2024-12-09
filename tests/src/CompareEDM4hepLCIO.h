#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H

#include "ComparisonUtils.h"
#include "ObjectMapping.h"

#include "edm4hep/CaloHitContributionCollection.h"
#include "edm4hep/CaloHitMCParticleLinkCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ClusterMCParticleLinkCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/VertexRecoParticleLinkCollection.h"
#include <edm4hep/VertexRecoParticleLink.h>
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/VertexCollection.h"
#include "edm4hep/utils/ParticleIDUtils.h"

#include "podio/Frame.h"

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
#include <UTIL/PIDHandler.h>
#include <lcio.h>

bool compare(const EVENT::CalorimeterHit* lcio, const edm4hep::CalorimeterHit& edm4hep,
             const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::CalorimeterHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::Cluster* lcio, const edm4hep::Cluster& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::ClusterCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::MCParticle* lcio, const edm4hep::MCParticle& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::MCParticleCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::RawCalorimeterHit* lcio, const edm4hep::RawCalorimeterHit& edm4hep,
             const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::RawCalorimeterHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::ReconstructedParticle* lcio, const edm4hep::ReconstructedParticle& edm4hep,
             const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection,
             const edm4hep::ReconstructedParticleCollection& edm4hepCollection, const ObjectMappings& objectMaps);

bool compare(const EVENT::SimCalorimeterHit* lcio, const edm4hep::SimCalorimeterHit& edm4hep,
             const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::SimCalorimeterHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::SimTrackerHit* lcio, const edm4hep::SimTrackerHit& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::SimTrackerHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::TPCHit* lcio, const edm4hep::RawTimeSeries& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::RawTimeSeriesCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::TrackerHit* lcio, const edm4hep::TrackerHit3D& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackerHit3DCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::TrackerHitPlane* lcio, const edm4hep::TrackerHitPlane& edm4hep,
             const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackerHitPlaneCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::TrackState* lcio, const edm4hep::TrackState& edm4hep);

bool compare(const EVENT::Track* lcio, const edm4hep::Track& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::Vertex* lcio, const edm4hep::Vertex& edm4hep, const ObjectMappings& objectMaps);
bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::VertexCollection& edm4hepCollection,
             const ObjectMappings& objectMaps);

bool compare(const EVENT::ParticleID* lcio, const edm4hep::ParticleID& edm4hep);

bool compareEventHeader(const EVENT::LCEvent* lcevt, const podio::Frame* edmEvent);

template <typename LinkCollT>
bool compare(const lcio::LCCollection* lcioCollection, const LinkCollT& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::LCRelation>(lcioCollection, edm4hepCollection, objectMaps);
}

namespace detail {
template <typename LinkT>
struct LcioFromToTypeHelper;

template <>
struct LcioFromToTypeHelper<edm4hep::RecoMCParticleLink> {
  using from_type = EVENT::ReconstructedParticle;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::CaloHitSimCaloHitLink> {
  using from_type = EVENT::CalorimeterHit;
  using to_type = EVENT::SimCalorimeterHit;
};

template <>
struct LcioFromToTypeHelper<edm4hep::TrackerHitSimTrackerHitLink> {
  using from_type = EVENT::TrackerHit;
  using to_type = EVENT::SimTrackerHit;
};

template <>
struct LcioFromToTypeHelper<edm4hep::CaloHitMCParticleLink> {
  using from_type = EVENT::CalorimeterHit;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::ClusterMCParticleLink> {
  using from_type = EVENT::Cluster;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::TrackMCParticleLink> {
  using from_type = EVENT::Track;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::VertexRecoParticleLink> {
  using from_type = EVENT::ReconstructedParticle;
  using to_type = EVENT::Vertex;
};

template <typename LinkT>
using getLcioFromType = typename LcioFromToTypeHelper<LinkT>::from_type;

template <typename LinkT>
using getLcioToType = typename LcioFromToTypeHelper<LinkT>::to_type;

template <typename T>
const auto& getObjectMap(const ObjectMappings& maps) {
  if constexpr (std::is_same_v<T, EVENT::MCParticle>) {
    return maps.mcParticles;
  } else if constexpr (std::is_same_v<T, EVENT::ReconstructedParticle>) {
    return maps.recoParticles;
  } else if constexpr (std::is_same_v<T, EVENT::SimCalorimeterHit>) {
    return maps.simCaloHits;
  } else if constexpr (std::is_same_v<T, EVENT::SimTrackerHit>) {
    return maps.simTrackerHits;
  } else if constexpr (std::is_same_v<T, EVENT::TrackerHit>) {
    return maps.trackerHits;
  } else if constexpr (std::is_same_v<T, EVENT::Cluster>) {
    return maps.clusters;
  } else if constexpr (std::is_same_v<T, EVENT::CalorimeterHit>) {
    return maps.caloHits;
  } else if constexpr (std::is_same_v<T, EVENT::Track>) {
    return maps.tracks;
  } else {
    return maps.vertices;
  }
}

} // namespace detail

template <typename LinkT>
bool compare(const EVENT::LCRelation* lcio, const LinkT& edm4hep, const ObjectMappings& objectMaps) {

  ASSERT_COMPARE(lcio, edm4hep, getWeight, "weight in relation / link");

  using LcioFromT = detail::getLcioFromType<LinkT>;
  using LcioToT = detail::getLcioToType<LinkT>;

  const auto lcioFrom = static_cast<LcioFromT*>(lcio->getFrom());
  const auto edm4hepFrom = edm4hep.getFrom();
  if (!compareRelation(lcioFrom, edm4hepFrom, detail::getObjectMap<LcioFromT>(objectMaps),
                       "from object in relation / link")) {
    return false;
  }

  const auto lcioTo = static_cast<LcioToT*>(lcio->getTo());
  const auto edm4hepTo = edm4hep.getTo();
  if (!compareRelation(lcioTo, edm4hepTo, detail::getObjectMap<LcioToT>(objectMaps), "to object in relation / link")) {
    return false;
  }

  return true;
}

/// Compare the information stored in startVertex in LCIO

bool compareStartVertexRelations(const EVENT::ReconstructedParticle* lcioReco,
                                 const edm4hep::VertexRecoParticleLink& link, const ObjectMappings& objectMaps);

/// Compare the information stored in associatedParticle in LCIO
bool compareVertexRecoLink(const EVENT::Vertex* lcioVtx, const edm4hep::VertexRecoParticleLink& link,
                           const ObjectMappings& objectMaps);

#define ASSERT_COMPARE_OR_EXIT(collType)                                                                               \
  if (type == #collType) {                                                                                             \
    auto& edmcoll = edmEvent.get<collType>(name);                                                                      \
    if (!compare(lcioColl, edmcoll, objectMapping)) {                                                                  \
      std::cerr << "in collection: " << name << std::endl;                                                             \
      return 1;                                                                                                        \
    }                                                                                                                  \
  }

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
