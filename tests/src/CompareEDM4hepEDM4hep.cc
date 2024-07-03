#include "CompareEDM4hepEDM4hep.h"
#include "ComparisonUtils.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/RecDqdxCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include <edm4hep/TrackState.h>
#include <iostream>
#include <podio/RelationRange.h>

#define REQUIRE_SAME(expected, actual, msg)                                                                            \
  {                                                                                                                    \
    if (!((expected) == (actual))) {                                                                                   \
      std::cerr << msg << " are not the same (expected: " << (expected) << ", actual: " << (actual) << ")"             \
                << std::endl;                                                                                          \
      return false;                                                                                                    \
    }                                                                                                                  \
  }

bool compare(const edm4hep::CalorimeterHitCollection& origColl,
             const edm4hep::CalorimeterHitCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    const auto origHit = origColl[i];
    const auto hit = roundtripColl[i];

    REQUIRE_SAME(origHit.getCellID(), hit.getCellID(), "cellID in hit " << i);
    REQUIRE_SAME(origHit.getEnergy(), hit.getEnergy(), "energy in hit " << i);
    REQUIRE_SAME(origHit.getEnergyError(), hit.getEnergyError(), "energyError in hit " << i);
    REQUIRE_SAME(origHit.getTime(), hit.getTime(), "time in hit " << i);
    REQUIRE_SAME(origHit.getPosition(), hit.getPosition(), "position in hit " << i);
  }

  return true;
}

bool compare(const edm4hep::MCParticleCollection& origColl, const edm4hep::MCParticleCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");

  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origPart = origColl[i];
    auto part = roundtripColl[i];

    REQUIRE_SAME(origPart.getPDG(), part.getPDG(), "pdg in particle " << i);
    REQUIRE_SAME(origPart.getGeneratorStatus(), part.getGeneratorStatus(), "generatorStatus in particle " << i);
    REQUIRE_SAME(origPart.getVertex(), part.getVertex(), "vertex in particle " << i);
    REQUIRE_SAME(origPart.getTime(), part.getTime(), "time in particle " << i);
    REQUIRE_SAME(origPart.getEndpoint(), part.getEndpoint(), "endpoint in particle " << i);
    REQUIRE_SAME(origPart.getMomentum(), part.getMomentum(), "momentum in particle " << i);
    REQUIRE_SAME(origPart.getMomentumAtEndpoint(), part.getMomentumAtEndpoint(),
                 "momentumAtEndpoint in particle " << i);
    REQUIRE_SAME(origPart.getMass(), part.getMass(), "mass in particle " << i);
    REQUIRE_SAME(origPart.getCharge(), part.getCharge(), "charge in particle " << i);
    REQUIRE_SAME(origPart.getSpin(), part.getSpin(), "spin in particle " << i);
    REQUIRE_SAME(origPart.getColorFlow(), part.getColorFlow(), "colorFlow in particle " << i);

    REQUIRE_SAME(origPart.isCreatedInSimulation(), part.isCreatedInSimulation(), " in particle " << i);
    REQUIRE_SAME(origPart.isBackscatter(), part.isBackscatter(), " in particle " << i);
    REQUIRE_SAME(origPart.vertexIsNotEndpointOfParent(), part.vertexIsNotEndpointOfParent(), " in particle " << i);
    REQUIRE_SAME(origPart.isDecayedInTracker(), part.isDecayedInTracker(), " in particle " << i);
    REQUIRE_SAME(origPart.isDecayedInCalorimeter(), part.isDecayedInCalorimeter(), " in particle " << i);
    REQUIRE_SAME(origPart.hasLeftDetector(), part.hasLeftDetector(), " in particle " << i);
    REQUIRE_SAME(origPart.isStopped(), part.isStopped(), " in particle " << i);
    REQUIRE_SAME(origPart.isOverlay(), part.isOverlay(), " in particle " << i);

    // Check the mother / daughter relations. We use the ObjectIDs here
    // (assuming that the collection names are the same!)
    REQUIRE_SAME(origPart.getParents().size(), part.getParents().size(), "size of parents in particle " << i);
    for (size_t iP = 0; iP < origPart.getParents().size(); ++iP) {
      REQUIRE_SAME(origPart.getParents()[iP].getObjectID(), part.getParents()[iP].getObjectID(),
                   " parent " << iP << " in particle " << i);
    }
    REQUIRE_SAME(origPart.getDaughters().size(), part.getDaughters().size(), "size of daughters in particle " << i);
    for (size_t iD = 0; iD < origPart.getDaughters().size(); ++iD) {
      REQUIRE_SAME(origPart.getDaughters()[iD].getObjectID(), part.getDaughters()[iD].getObjectID(),
                   " daughter " << iD << " in particle " << i);
    }
  }

  return true;
}

bool compare(const edm4hep::SimCalorimeterHitCollection& origColl,
             const edm4hep::SimCalorimeterHitCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origHit = origColl[i];
    auto hit = roundtripColl[i];

    REQUIRE_SAME(origHit.getCellID(), hit.getCellID(), "cellID in hit " << i);
    REQUIRE_SAME(origHit.getEnergy(), hit.getEnergy(), "energy in hit " << i);
    REQUIRE_SAME(origHit.getPosition(), hit.getPosition(), "position in hit " << i);

    // Check the contributions
    const auto origContribs = origHit.getContributions();
    const auto contribs = hit.getContributions();
    REQUIRE_SAME(origContribs.size(), contribs.size(), "size of contributions in hit " << i);
    for (size_t iC = 0; iC < origHit.getContributions().size(); ++iC) {
      const auto origCont = origContribs[iC];
      const auto cont = contribs[iC];
      REQUIRE_SAME(origCont.getPDG(), cont.getPDG(), "pdg of contribution " << iC << " in hit " << i);
      REQUIRE_SAME(origCont.getEnergy(), cont.getEnergy(), "energy of contribution " << iC << " in hit " << i);
      REQUIRE_SAME(origCont.getTime(), cont.getTime(), "time of contribution " << iC << " in hit " << i);
      REQUIRE_SAME(origCont.getStepPosition(), cont.getStepPosition(),
                   "stepPosition of contribution " << iC << " in hit " << i);

      // Check the MCParticles via ObjectID (asssuming collection names remain
      // unchanged)
      REQUIRE_SAME(origCont.getParticle().getObjectID(), cont.getParticle().getObjectID(),
                   "particle of contribution " << iC << " in hit " << i);
    }
  }

  return true;
}

bool compare(const edm4hep::TrackState& orig, const edm4hep::TrackState& roundtrip) {
  REQUIRE_SAME(orig.location, roundtrip.location, "location in TrackState");
  REQUIRE_SAME(orig.D0, roundtrip.D0, "D0 in TrackState");
  REQUIRE_SAME(orig.Z0, roundtrip.Z0, "Z0 in TrackState");
  REQUIRE_SAME(orig.phi, roundtrip.phi, "phi in TrackState");
  REQUIRE_SAME(orig.omega, roundtrip.omega, "omega in TrackState");
  REQUIRE_SAME(orig.tanLambda, roundtrip.tanLambda, "tanLambda in TrackState");
  REQUIRE_SAME(orig.referencePoint, roundtrip.referencePoint, "referencePoint in TrackState");
  // The time related covariance matrix elements are lost in the EDM4hep -> LCIO
  // conversion and set to zero in the LCIO -> EDM4hep conversion
  const auto redOrigCov = std::vector(orig.covMatrix.begin(), orig.covMatrix.begin() + 15);
  const auto redCov = std::vector(roundtrip.covMatrix.begin(), roundtrip.covMatrix.begin() + 15);
  REQUIRE_SAME(redOrigCov, redCov, "covMatrix in TrackState");

  return true;
}

bool compare(const edm4hep::TrackCollection& origColl, const edm4hep::TrackCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");

  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origTrack = origColl[i];
    auto track = roundtripColl[i];

    REQUIRE_SAME(origTrack.getChi2(), track.getChi2(), "chi2 in track " << i);
    REQUIRE_SAME(origTrack.getNdf(), track.getNdf(), "chi2 in track " << i);

    const auto origSubDetHitNumbers = origTrack.getSubdetectorHitNumbers();
    const auto subDetHitNumbers = track.getSubdetectorHitNumbers();
    // Resizing in EDM4hep to LCIO conversion causes the roundtripped version to
    // have 50 hits. So here we just compare the ones that are available in the
    // original
    for (size_t iSN = 0; iSN < origSubDetHitNumbers.size(); ++iSN) {
      REQUIRE_SAME(origSubDetHitNumbers[iSN], subDetHitNumbers[iSN],
                   "subdetector hit numbers " << iSN << " in track " << i);
    }
    for (size_t iSN = origSubDetHitNumbers.size(); iSN < 50; ++iSN) {
      REQUIRE_SAME(0, subDetHitNumbers[iSN], "additional subdetector hit number in track " << i);
    }

    const auto origHits = origTrack.getTrackerHits();
    const auto hits = track.getTrackerHits();
    REQUIRE_SAME(origHits.size(), hits.size(), "number of hits in track " << i);
    for (size_t iH = 0; iH < origHits.size(); ++iH) {
      REQUIRE_SAME(origHits[iH].getObjectID(), hits[iH].getObjectID(), "tracker hit " << iH << " in track " << i);
    }

    const auto origStates = origTrack.getTrackStates();
    const auto states = track.getTrackStates();
    REQUIRE_SAME(origStates.size(), states.size(), "number of track states in track " << i);
    for (size_t iS = 0; iS < origStates.size(); ++iS) {
      if (!compare(origStates[iS], states[iS])) {
        std::cerr << "in track state " << iS << " in track " << i << std::endl;
        return false;
      }
    }

    const auto origTracks = origTrack.getTracks();
    const auto tracks = track.getTracks();
    REQUIRE_SAME(origTracks.size(), tracks.size(), "number of tracks in track " << i);
    for (size_t iT = 0; iT < origTracks.size(); ++iT) {
      REQUIRE_SAME(origTracks[iT].getObjectID(), tracks[iT].getObjectID(), "track " << iT << " in track " << i);
    }
  }
  return true;
}

bool compare(const edm4hep::TrackerHit3DCollection& origColl, const edm4hep::TrackerHit3DCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origHit = origColl[i];
    auto hit = roundtripColl[i];

    REQUIRE_SAME(origHit.getCellID(), hit.getCellID(), "cellID in hit " << i);
    REQUIRE_SAME(origHit.getType(), hit.getType(), "type in hit " << i);
    REQUIRE_SAME(origHit.getQuality(), hit.getQuality(), "quality in hit " << i);
    REQUIRE_SAME(origHit.getTime(), hit.getTime(), "time in hit " << i);
    REQUIRE_SAME(origHit.getEDep(), hit.getEDep(), "EDep in hit " << i);
    REQUIRE_SAME(origHit.getEDepError(), hit.getEDepError(), "EDepError in hit " << i);
    REQUIRE_SAME(origHit.getPosition(), hit.getPosition(), "Position in hit " << i);
    REQUIRE_SAME(origHit.getCovMatrix(), hit.getCovMatrix(), "CovMatrix in hit " << i);
  }

  return true;
}

bool compare(const edm4hep::TrackerHitPlaneCollection& origColl,
             const edm4hep::TrackerHitPlaneCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origHit = origColl[i];
    auto hit = roundtripColl[i];

    REQUIRE_SAME(origHit.getCellID(), hit.getCellID(), "cellID in hit " << i);
    REQUIRE_SAME(origHit.getType(), hit.getType(), "type in hit " << i);
    REQUIRE_SAME(origHit.getQuality(), hit.getQuality(), "quality in hit " << i);
    REQUIRE_SAME(origHit.getTime(), hit.getTime(), "time in hit " << i);
    REQUIRE_SAME(origHit.getEDep(), hit.getEDep(), "EDep in hit " << i);
    REQUIRE_SAME(origHit.getEDepError(), hit.getEDepError(), "EDepError in hit " << i);
    REQUIRE_SAME(origHit.getPosition(), hit.getPosition(), "Position in hit " << i);
    REQUIRE_SAME(origHit.getCovMatrix(), hit.getCovMatrix(), "CovMatrix in hit " << i);
    REQUIRE_SAME(origHit.getU(), hit.getU(), "U in hit " << i);
    REQUIRE_SAME(origHit.getV(), hit.getV(), "V in hit " << i);
    REQUIRE_SAME(origHit.getDu(), hit.getDu(), "Du in hit " << i);
    REQUIRE_SAME(origHit.getDv(), hit.getDv(), "Dv in hit " << i);
  }

  return true;
}

bool compare(const edm4hep::ClusterCollection& origColl, const edm4hep::ClusterCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origCluster = origColl[i];
    auto cluster = roundtripColl[i];

    const auto origRelClusters = origCluster.getClusters();
    const auto relClusters = cluster.getClusters();
    REQUIRE_SAME(origRelClusters.size(), relClusters.size(), "number of related clusters in cluster " << i);
    for (size_t iC = 0; iC < origRelClusters.size(); ++iC) {
      REQUIRE_SAME(origRelClusters[iC].getObjectID(), relClusters[iC].getObjectID(),
                   "related cluster " << iC << " in cluster " << i);
    }

    const auto origHits = origCluster.getHits();
    const auto hits = cluster.getHits();
    REQUIRE_SAME(origHits.size(), hits.size(), "number of calorimeter hits in cluster " << i);
    for (size_t iH = 0; iH < origHits.size(); ++iH) {
      REQUIRE_SAME(origHits[iH].getObjectID(), hits[iH].getObjectID(), "calorimeter hit " << iH << " in cluster " << i);
    }

    const auto& origSubdetE = origCluster.getSubdetectorEnergies();
    const auto& subdetE = cluster.getSubdetectorEnergies();
    REQUIRE_SAME(origSubdetE.size(), subdetE.size(), "sizes of subdetector energies in cluster " << i);
    for (size_t iSE = 0; iSE < origSubdetE.size(); ++iSE) {
      REQUIRE_SAME(origSubdetE[iSE], subdetE[iSE], "subdetector energy " << iSE << " in cluster " << i);
    }
  }

  return true;
}

bool compare(const edm4hep::ReconstructedParticleCollection& origColl,
             const edm4hep::ReconstructedParticleCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    auto origReco = origColl[i];
    auto reco = roundtripColl[i];

    REQUIRE_SAME(origReco.getCharge(), reco.getCharge(), "charge in reco particle " << i);
    REQUIRE_SAME(origReco.getMomentum(), reco.getMomentum(), "momentum in reco particle " << i);

    const auto origRelClusters = origReco.getClusters();
    const auto relClusters = reco.getClusters();
    REQUIRE_SAME(origRelClusters.size(), relClusters.size(), "number of related clusters in reco particle " << i);
    for (size_t iC = 0; iC < relClusters.size(); ++iC) {
      REQUIRE_SAME(origRelClusters[iC].getObjectID(), relClusters[iC].getObjectID(),
                   "related cluster " << iC << " in reco particle " << i);
    }

    const auto origRelTracks = origReco.getTracks();
    const auto relTracks = reco.getTracks();
    REQUIRE_SAME(origRelTracks.size(), relTracks.size(), "number of related tracks in reco particle " << i);
    for (size_t iT = 0; iT < relTracks.size(); ++iT) {
      REQUIRE_SAME(origRelTracks[iT].getObjectID(), relTracks[iT].getObjectID(),
                   "related track " << iT << " in reco particle " << i);
    }

    const auto origRelParticles = origReco.getParticles();
    const auto relParticles = reco.getParticles();
    REQUIRE_SAME(origRelParticles.size(), relParticles.size(), "number of related particles in reco particle " << i);
    for (size_t iP = 0; iP < relParticles.size(); ++iP) {
      REQUIRE_SAME(origRelParticles[iP].getObjectID(), relParticles[iP].getObjectID(),
                   "related particle " << iP << " in reco particle " << i);
    }
  }

  return true;
}

bool compare(const edm4hep::ParticleIDCollection& origColl, const edm4hep::ParticleIDCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    const auto origPid = origColl[i];
    const auto pid = roundtripColl[i];

    // This might not be preserved in roundtripping
    REQUIRE_SAME(origPid.getAlgorithmType(), pid.getAlgorithmType(), "algorithm type in ParticleID " << i);

    REQUIRE_SAME(origPid.getType(), pid.getType(), "type in ParticleID " << i);

    const auto origParams = origPid.getParameters();
    const auto params = pid.getParameters();
    REQUIRE_SAME(origParams.size(), params.size(), "parameter sizes in ParticleID " << i);
    for (size_t iP = 0; iP < params.size(); ++iP) {
      REQUIRE_SAME(origParams[iP], params[iP], "parameter " << iP << " in ParticleID " << i);
    }

    REQUIRE_SAME(origPid.getParticle().getObjectID(), pid.getParticle().getObjectID(),
                 "related particle in ParticleID " << i);
  }
  return true;
}

bool compare(const edm4hep::RecoParticleVertexAssociationCollection& origColl,
             const edm4hep::RecoParticleVertexAssociationCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    const auto origAssoc = origColl[i];
    const auto assoc = roundtripColl[i];

    REQUIRE_SAME(origAssoc.getWeight(), assoc.getWeight(), "weight in association " << i);
    REQUIRE_SAME(origAssoc.getVertex().id(), assoc.getVertex().id(), "vertex in association " << i);
    REQUIRE_SAME(origAssoc.getRec().id(), assoc.getRec().id(), "reco particle in association " << i);
  }

  return true;
}

bool compare(const edm4hep::RecDqdxCollection& origColl, const edm4hep::RecDqdxCollection& roundtripColl) {
  REQUIRE_SAME(origColl.size(), roundtripColl.size(), "collection sizes");
  for (size_t i = 0; i < origColl.size(); ++i) {
    const auto origDqdx = origColl[i];
    const auto dQdx = roundtripColl[i];

    const auto origQuant = origDqdx.getDQdx();
    const auto quant = dQdx.getDQdx();

    REQUIRE_SAME(origQuant.value, quant.value, "value of quantity in dQ/dx " << i);
    REQUIRE_SAME(origQuant.error, quant.error, "error of quantity in dQ/dx " << i);
    REQUIRE_SAME(origQuant.type, quant.value, "value of quantity in dQ/dx " << i);

    REQUIRE_SAME(origDqdx.getTrack().id(), dQdx.getTrack().id(), "related track in dQ/dx " << i);
  }

  return true;
}
