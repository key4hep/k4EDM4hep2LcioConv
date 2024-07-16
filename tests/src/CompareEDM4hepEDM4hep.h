#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPEDM4HEP_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPEDM4HEP_H

#include "ComparisonUtils.h"
#include "EDM4hep2LCIOUtilities.h"

#include <iostream>

#define REQUIRE_SAME(expected, actual, msg)                                                                            \
  {                                                                                                                    \
    if (!((expected) == (actual))) {                                                                                   \
      std::cerr << msg << " are not the same (expected: " << (expected) << ", actual: " << (actual) << ")"             \
                << std::endl;                                                                                          \
      return false;                                                                                                    \
    }                                                                                                                  \
  }

bool compare(const edm4hep::CalorimeterHitCollection& origColl, const edm4hep::CalorimeterHitCollection& roundtripColl);

bool compare(const edm4hep::MCParticleCollection& origColl, const edm4hep::MCParticleCollection& roundtripColl);

bool compare(const edm4hep::SimCalorimeterHitCollection& origColl,
             const edm4hep::SimCalorimeterHitCollection& roundtripColl);

bool compare(const edm4hep::TrackCollection& origColl, const edm4hep::TrackCollection& roundtripColl);

bool compare(const edm4hep::TrackerHit3DCollection& origColl, const edm4hep::TrackerHit3DCollection& roundtripColl);

bool compare(const edm4hep::TrackerHitPlaneCollection& origColl,
             const edm4hep::TrackerHitPlaneCollection& roundtripColl);

bool compare(const edm4hep::ClusterCollection& origColl, const edm4hep::ClusterCollection& roundtripColl);

bool compare(const edm4hep::ReconstructedParticleCollection& origColl,
             const edm4hep::ReconstructedParticleCollection& roundtripColl);

bool compare(const edm4hep::ParticleIDCollection& origColl, const edm4hep::ParticleIDCollection& roundtripColl);

bool compare(const edm4hep::RecoParticleVertexAssociationCollection& origColl,
             const edm4hep::RecoParticleVertexAssociationCollection& roundtripColl);

template <typename AssociationCollT>
bool compare(const AssociationCollT& origColl, const AssociationCollT& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    const auto origAssoc = origColl[i];
    const auto assoc = roundtripColl[i];

    REQUIRE_SAME(origAssoc.getWeight(), assoc.getWeight(), "weight in association " << i);
    REQUIRE_SAME(origAssoc.getSim().id(), assoc.getSim().id(), "MC part(icle) in association " << i);
    REQUIRE_SAME(origAssoc.getRec().id(), assoc.getRec().id(), "reco part(icle) in association " << i);
  }
  return true;
}

bool compare(const edm4hep::RecDqdxCollection& origColl, const edm4hep::RecDqdxCollection& roundtripColl);

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPEDM4HEP_H
