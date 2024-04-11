#include "CompareEDM4hepEDM4hep.h"
#include "EDM4hep2LCIOUtilities.h"

#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.h"
#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

#include "podio/Frame.h"

#include <iostream>

#define ASSERT_SAME_OR_ABORT(type, name)                                     \
  if (!compare(origEvent.get<type>(name), roundtripEvent.get<type>(name))) { \
    std::cerr << "Comparison failure in " << name << std::endl;              \
    return 1;                                                                \
  }

int main()
{
  const auto& [origEvent, metadata] = createExampleEvent();
  const auto lcioEvent = EDM4hep2LCIOConv::convertEvent(origEvent, metadata);
  const auto roundtripEvent = LCIO2EDM4hepConv::convertEvent(lcioEvent.get());

  ASSERT_SAME_OR_ABORT(edm4hep::CalorimeterHitCollection, "caloHits");
  ASSERT_SAME_OR_ABORT(edm4hep::MCParticleCollection, "mcParticles");
  ASSERT_SAME_OR_ABORT(edm4hep::SimCalorimeterHitCollection, "simCaloHits");
  ASSERT_SAME_OR_ABORT(edm4hep::TrackCollection, "tracks");
  ASSERT_SAME_OR_ABORT(edm4hep::TrackerHit3DCollection, "trackerHits");
  ASSERT_SAME_OR_ABORT(edm4hep::TrackerHitPlaneCollection, "trackerHitPlanes");
  ASSERT_SAME_OR_ABORT(edm4hep::ClusterCollection, "clusters");
  ASSERT_SAME_OR_ABORT(edm4hep::ReconstructedParticleCollection, "recos");
  ASSERT_SAME_OR_ABORT(edm4hep::ParticleIDCollection, "ParticleID_coll_1");
  ASSERT_SAME_OR_ABORT(edm4hep::ParticleIDCollection, "ParticleID_coll_2");
  ASSERT_SAME_OR_ABORT(edm4hep::ParticleIDCollection, "ParticleID_coll_3");

  return 0;
}
