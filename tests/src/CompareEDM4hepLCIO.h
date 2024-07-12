#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H

#include "ComparisonUtils.h"
#include "ObjectMapping.h"

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
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/RecoParticleVertexAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
using TrackerHit3D = edm4hep::TrackerHit;
} // namespace edm4hep
#endif
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

template <typename AssocCollT>
bool compare(const lcio::LCCollection* lcioCollection, const AssocCollT& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::LCRelation>(lcioCollection, edm4hepCollection, objectMaps);
}

namespace detail {
template <typename AssocT>
struct LcioFromToTypeHelper;

template <>
struct LcioFromToTypeHelper<edm4hep::MCRecoParticleAssociation> {
  using from_type = EVENT::ReconstructedParticle;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::MCRecoCaloAssociation> {
  using from_type = EVENT::CalorimeterHit;
  using to_type = EVENT::SimCalorimeterHit;
};

template <>
struct LcioFromToTypeHelper<edm4hep::MCRecoTrackerAssociation> {
  using from_type = EVENT::TrackerHit;
  using to_type = EVENT::SimTrackerHit;
};

template <>
struct LcioFromToTypeHelper<edm4hep::MCRecoCaloParticleAssociation> {
  using from_type = EVENT::CalorimeterHit;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::MCRecoClusterParticleAssociation> {
  using from_type = EVENT::Cluster;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::MCRecoTrackParticleAssociation> {
  using from_type = EVENT::Track;
  using to_type = EVENT::MCParticle;
};

template <>
struct LcioFromToTypeHelper<edm4hep::RecoParticleVertexAssociation> {
  using from_type = EVENT::ReconstructedParticle;
  using to_type = EVENT::Vertex;
};

template <typename AssocT>
using getLcioFromType = typename LcioFromToTypeHelper<AssocT>::from_type;

template <typename AssocT>
using getLcioToType = typename LcioFromToTypeHelper<AssocT>::to_type;

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

template <typename AssocT>
bool compare(const EVENT::LCRelation* lcio, const AssocT& edm4hep, const ObjectMappings& objectMaps) {

  ASSERT_COMPARE(lcio, edm4hep, getWeight, "weight in relation / association");

  using LcioFromT = detail::getLcioFromType<AssocT>;
  using LcioToT = detail::getLcioToType<AssocT>;

  const auto lcioFrom = static_cast<LcioFromT*>(lcio->getFrom());
  const auto edm4hepFrom = edm4hep.getRec();
  if (!compareRelation(lcioFrom, edm4hepFrom, detail::getObjectMap<LcioFromT>(objectMaps),
                       "from / rec object in relation / association")) {
    return false;
  }

  const auto lcioTo = static_cast<LcioToT*>(lcio->getTo());
  if constexpr (std::is_same_v<AssocT, edm4hep::RecoParticleVertexAssociation>) {
    const auto edm4hepTo = edm4hep.getVertex();
    if (!compareRelation(lcioTo, edm4hepTo, detail::getObjectMap<LcioToT>(objectMaps),
                         "vertex object in relation / association")) {
      return false;
    }
  } else {
    const auto edm4hepTo = edm4hep.getSim();
    if (!compareRelation(lcioTo, edm4hepTo, detail::getObjectMap<LcioToT>(objectMaps),
                         "to / mc object in relation / association")) {
      return false;
    }
  }

  return true;
}

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPLCIO_H
