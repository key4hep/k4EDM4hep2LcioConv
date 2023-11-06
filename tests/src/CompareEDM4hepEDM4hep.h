#ifndef K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPEDM4HEP_H
#define K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPEDM4HEP_H

#include "EDM4hep2LCIOUtilities.h"

bool compare(const edm4hep::CalorimeterHitCollection& origColl, const edm4hep::CalorimeterHitCollection& roundtripColl);

bool compare(const edm4hep::MCParticleCollection& origColl, const edm4hep::MCParticleCollection& roundtripColl);

bool compare(
  const edm4hep::SimCalorimeterHitCollection& origColl,
  const edm4hep::SimCalorimeterHitCollection& roundtripColl);

bool compare(const edm4hep::TrackCollection& origColl, const edm4hep::TrackCollection& roundtripColl);

bool compare(const edm4hep::TrackerHitCollection& origColl, const edm4hep::TrackerHitCollection& roundtripColl);

bool compare(const edm4hep::ClusterCollection& origColl, const edm4hep::ClusterCollection& roundtripColl);

#endif // K4EDM4HEP2LCIOCONV_TEST_COMPAREEDM4HEPEDM4HEP_H
