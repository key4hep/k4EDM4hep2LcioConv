#ifndef K4EDM4HEP2LCIOCONV_TEST_OBJECTMAPPINGS_H
#define K4EDM4HEP2LCIOCONV_TEST_OBJECTMAPPINGS_H

#include <edm4hep/ParticleID.h>
#include <edm4hep/utils/TrackUtils.h>

#include "podio/ObjectID.h"

#include <unordered_map>

namespace EVENT {
class Track;
class TrackerHit;
class TrackerHitPlane;
class SimTrackerHit;
class CalorimeterHit;
class RawCalorimeterHit;
class SimCalorimeterHit;
class TPCHit;
class Cluster;
class Vertex;
class ReconstructedParticle;
class MCParticle;
class LCEvent;
class ParticleID;
} // namespace EVENT

namespace podio {
class Frame;
} // namespace podio

struct ObjectMappings {
  template <typename K>
  using Map = std::unordered_map<K, podio::ObjectID>;

  Map<const EVENT::Track*> tracks{};
  Map<const EVENT::TrackerHit*> trackerHits{};
  Map<const EVENT::TrackerHitPlane*> trackerHitPlanes{};
  Map<const EVENT::SimTrackerHit*> simTrackerHits{};
  Map<const EVENT::Cluster*> clusters{};
  Map<const EVENT::CalorimeterHit*> caloHits{};
  Map<const EVENT::RawCalorimeterHit*> rawCaloHits{};
  Map<const EVENT::SimCalorimeterHit*> simCaloHits{};
  Map<const EVENT::MCParticle*> mcParticles{};
  Map<const EVENT::ReconstructedParticle*> recoParticles{};
  Map<const EVENT::TPCHit*> tpcHits{};
  Map<const EVENT::Vertex*> vertices{};

  std::unordered_map<const EVENT::ParticleID*, edm4hep::ParticleID> particleIDs{};

  edm4hep::utils::TrackPIDHandler trackPidHandler{};

  static ObjectMappings fromEvent(EVENT::LCEvent* lcEvt, const podio::Frame& edmEvt);
};

#endif
