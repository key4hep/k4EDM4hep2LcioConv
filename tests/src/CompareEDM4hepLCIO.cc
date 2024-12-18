#include "CompareEDM4hepLCIO.h"
#include "ComparisonUtils.h"
#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.h"

#include "IMPL/TrackerHitImpl.h"

#include <cmath>
#include <cstdint>
#include <edm4hep/VertexRecoParticleLinkCollection.h>

#include "TMath.h"

/**
 * The basic implementation of the functionality has been generated via modified
 * podio templates, employing some handwritten macros to facilitate the task.
 * These macros and the code generation cover the data members of all involved
 * data types but no relations.
 *
 * The basic working principle is to implement an overload that compares an LCIO
 * and EDM4hep object for each corresponding type. A second overload for
 * collections of the involved types simply forwards to the templated
 * compareCollection function, which then does the looping and eventually calls
 * the single object compare. All comparison functions return true if all
 * comparisons succeeded or false in case any comparison failed returning as
 * early as possible.
 */

/// Convert the two 32 bit cellIDs into one 64 bit value
template <typename LcioT>
auto to64BitCellID(LcioT* obj) {
  const auto cellID0 = obj->getCellID0();
  const auto cellID1 = obj->getCellID1();
  uint64_t cellID = cellID1;
  cellID = (cellID << 32) | cellID0;
  return cellID;
}

// ================= CalorimeterHit ================

bool compare(const EVENT::CalorimeterHit* lcioElem, const edm4hep::CalorimeterHit& edm4hepElem, const ObjectMappings&) {
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in CalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergyError, "energyError in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in CalorimeterHit");
  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::CalorimeterHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::CalorimeterHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= Cluster ================

bool compare(const EVENT::Cluster* lcioElem, const edm4hep::Cluster& edm4hepElem, const ObjectMappings& objectMaps) {
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in Cluster");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in Cluster");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergyError, "energyError in Cluster");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in Cluster");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPositionError, "positionError in Cluster");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getITheta, "iTheta in Cluster");
  // LCIO has getIPhi and EDM4hep has getPhi
  ASSERT_COMPARE_VALS(lcioElem->getIPhi(), edm4hepElem.getPhi(), "phi in Cluster");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getDirectionError, "directionError in Cluster");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getSubdetectorEnergies, "subdetectorEnergies in Cluster");
  ASSERT_COMPARE_VALS(lcioElem->getShape(), edm4hepElem.getShapeParameters(), "shape / shapeParameters in Cluster");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getClusters, objectMaps.clusters, "related clusters in Cluster");

  // Different names of related calorimeter hits in interfaces
  if (!compareRelation(lcioElem->getCalorimeterHits(), edm4hepElem.getHits(), objectMaps.caloHits,
                       "calorimeter hits in Cluster")) {
    return false;
  }

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::ClusterCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::Cluster>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= MCParticle ================

bool compare(const EVENT::MCParticle* lcioElem, const edm4hep::MCParticle& edm4hepElem,
             const ObjectMappings& objectMaps) {
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPDG, "PDG in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getGeneratorStatus, "generatorStatus in MCParticle");
  // LCIO changes the SimulatorStatus during I/O, so here we have to check the
  // individual bits  which are untouched instead of just doing one comparison for the SimulatorStatus
  ASSERT_COMPARE(lcioElem, edm4hepElem, isCreatedInSimulation, "Created in Simulation");
  ASSERT_COMPARE(lcioElem, edm4hepElem, isBackscatter, "particle is from backscatter of calorimeter shower");
  ASSERT_COMPARE(lcioElem, edm4hepElem, vertexIsNotEndpointOfParent, "checks vertex, enpoint of parent");
  ASSERT_COMPARE(lcioElem, edm4hepElem, isDecayedInTracker, "isDecayedInTracker");
  ASSERT_COMPARE(lcioElem, edm4hepElem, isDecayedInCalorimeter, "isDecayedInCalorimeter,");
  ASSERT_COMPARE(lcioElem, edm4hepElem, hasLeftDetector, "hasLeftDetector,");
  ASSERT_COMPARE(lcioElem, edm4hepElem, isStopped, "isStopped,");
  ASSERT_COMPARE(lcioElem, edm4hepElem, isOverlay, "isOverlay,");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getCharge, "charge in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMass, "mass in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getVertex, "vertex in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEndpoint, "endpoint in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMomentum, "momentum in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMomentumAtEndpoint, "momentumAtEndpoint in MCParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getSpin, "spin in MCParticle");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getDaughters, objectMaps.mcParticles, "daughters in MCParticle");
  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getParents, objectMaps.mcParticles, "parents in MCParticle");

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::MCParticleCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::MCParticle>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= ParticleID ================

bool compare(const EVENT::ParticleID* lcioElem, const edm4hep::ParticleID& edm4hepElem) {
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in ParticleID");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPDG, "PDG in ParticleID");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getAlgorithmType, "algorithmType in ParticleID");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getLikelihood, "likelihood in ParticleID");
  return true;
}

// ================= RawCalorimeterHit ================

bool compare(const EVENT::RawCalorimeterHit* lcioElem, const edm4hep::RawCalorimeterHit& edm4hepElem,
             const ObjectMappings&) {
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in RawCalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getAmplitude, "amplitude in RawCalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTimeStamp, "timeStamp in RawCalorimeterHit");
  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::RawCalorimeterHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::RawCalorimeterHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= ReconstructedParticle ================

bool compare(const EVENT::ReconstructedParticle* lcioElem, const edm4hep::ReconstructedParticle& edm4hepElem,
             const ObjectMappings& objectMaps) {
  ASSERT_COMPARE_VALS(lcioElem->getType(), edm4hepElem.getPDG(), "type/PDG in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMomentum, "momentum in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getReferencePoint, "referencePoint in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCharge, "charge in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMass, "mass in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getGoodnessOfPID, "goodnessOfPID in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCovMatrix, "covMatrix in ReconstructedParticle");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getClusters, objectMaps.clusters, "clusters in ReonstructedParticle");
  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getTracks, objectMaps.tracks, "tracks in ReonstructedParticle");

  // The attached particles need special treatment because there is a difference
  // in EDM4hep and LCIO regarding where the decay particles of a Vertex go
  const auto edm4hepDecayVtx = edm4hepElem.getDecayVertex();
  if (edm4hepDecayVtx.isAvailable()) {
    if (lcioElem->getParticles().size() == edm4hepDecayVtx.getParticles().size()) {
      ASSERT_COMPARE_RELATION(lcioElem, edm4hepDecayVtx, getParticles, objectMaps.recoParticles,
                              "decay particles in ReconstructedParticle (decayVertex)");
    } else {
      // This doesn't seem to happen, in case it ever does at least make the test fail
      std::cerr << "Potentially inconsistent setting of (decay) particles in ReconstructedParticle" << std::endl;
      return false;
    }
  }

  // TODO: check for start vertex. Might not be possible here but needs to be done externally

  // ParticleIDs need special treatment because they live in different
  // collections in EDM4hep. Here we make sure that all ParticleIDs have been
  // converted and mapped correctly by checking the ParticleID content and
  // making sure that it points back to the converted reco particle
  const auto& lcioPIDs = lcioElem->getParticleIDs();
  for (size_t i = 0; i < lcioPIDs.size(); ++i) {
    if (auto it = objectMaps.particleIDs.find(lcioPIDs[i]); it != objectMaps.particleIDs.end()) {
      const auto& [lcioPid, edm4hepPid] = *it;
      if (!compare(lcioPid, edm4hepPid) || edm4hepPid.getParticle() != edm4hepElem) {
        std::cerr << "particle ID " << i << " is not mapped to the correct EDM4hep particle ID (LCIO: " << lcioPid
                  << ", EDM4hep: " << edm4hepPid << ")" << std::endl;
        return false;
      }
    } else {
      std::cerr << "Cannot find a converted ParticleID object for particle ID " << i << std::endl;
      return false;
    }
  }

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection,
             const edm4hep::ReconstructedParticleCollection& edm4hepCollection, const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::ReconstructedParticle>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= SimCalorimeterHit ================

bool compare(const EVENT::SimCalorimeterHit* lcioElem, const edm4hep::SimCalorimeterHit& edm4hepElem,
             const ObjectMappings& objectMaps) {
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in SimCalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in SimCalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in SimCalorimeterHit");

  // Contributions are not part of the "proper LCIO"
  const auto edmContributions = edm4hepElem.getContributions();
  ASSERT_COMPARE_VALS((unsigned)lcioElem->getNMCContributions(), edmContributions.size(),
                      "number of CaloHitContributions");

  for (int iCont = 0; iCont < lcioElem->getNMCContributions(); ++iCont) {
    const auto& edmContrib = edmContributions[iCont];
    ASSERT_COMPARE_VALS(lcioElem->getEnergyCont(iCont), edmContrib.getEnergy(),
                        "energy in CaloHitContribution " + std::to_string(iCont));
    ASSERT_COMPARE_VALS(lcioElem->getStepPosition(iCont), edmContrib.getStepPosition(),
                        "stepPosition in CaloHitContribution " + std::to_string(iCont));
    ASSERT_COMPARE_VALS(lcioElem->getTimeCont(iCont), edmContrib.getTime(),
                        "time in CaloHitContribution " + std::to_string(iCont));

    if (!compareRelation(lcioElem->getParticleCont(iCont), edmContrib.getParticle(), objectMaps.mcParticles,
                         " MCParticle in CaloHitContribution " + std::to_string(iCont))) {
      return false;
    }
  }

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::SimCalorimeterHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::SimCalorimeterHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= SimTrackerHit ================

bool compare(const EVENT::SimTrackerHit* lcioElem, const edm4hep::SimTrackerHit& edm4hepElem,
             const ObjectMappings& objectMaps) {
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEDep, "EDep in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPathLength, "pathLength in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getQuality, "quality in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMomentum, "momentum in SimTrackerHit");

  if (!compareRelation(lcioElem->getMCParticle(), edm4hepElem.getParticle(), objectMaps.mcParticles,
                       "MC particle in SimTrackerHit")) {
    return false;
  }

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::SimTrackerHitCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::SimTrackerHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TPCHit ================

bool compare(const EVENT::TPCHit* lcioElem, const edm4hep::RawTimeSeries& edm4hepElem, const ObjectMappings&) {
  const uint64_t lcioCellID = lcioElem->getCellID();
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in TPCHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getQuality, "quality in TPCHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in TPCHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCharge, "charge in TPCHit");
  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::RawTimeSeriesCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::TPCHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TrackState ================

bool compare(const EVENT::TrackState* lcio, const edm4hep::TrackState& edm4hep) {
  ASSERT_COMPARE_VALS(lcio->getLocation(), edm4hep.location, "location in TrackState");
  ASSERT_COMPARE_VALS(lcio->getD0(), edm4hep.D0, "D0 in TrackState");
  ASSERT_COMPARE_VALS(lcio->getZ0(), edm4hep.Z0, "Z0 in TrackState");
  ASSERT_COMPARE_VALS(lcio->getPhi(), edm4hep.phi, "phi in TrackState");
  ASSERT_COMPARE_VALS(lcio->getOmega(), edm4hep.omega, "omega in TrackState");
  ASSERT_COMPARE_VALS(lcio->getTanLambda(), edm4hep.tanLambda, "tanLambda in TrackState");
  ASSERT_COMPARE_VALS(lcio->getReferencePoint(), edm4hep.referencePoint, "referencePoint in TrackState");

  // Need to make sure to only compare the part of the covariance matrix that is
  // actually available in LCIO
  const auto& lcioCov = lcio->getCovMatrix();
  const auto& edmCov = edm4hep.covMatrix;
  ASSERT_COMPARE_VALS(lcioCov, std::vector(edmCov.begin(), edmCov.begin() + 15), "covMatrix in TrackState");

  return true;
}

// ================= Track ================

bool compare(const EVENT::Track* lcioElem, const edm4hep::Track& edm4hepElem, const ObjectMappings& objectMaps) {
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in Track");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getChi2, "chi2 in Track");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getNdf, "ndf in Track");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getNholes, "Nholes in Track");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getSubdetectorHitNumbers, "subdetectorHitNumbers");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getSubdetectorHoleNumbers, "subdetectorHitNumbers");

  // EDM4hep does not have the dEdx information inside the track, but rather
  // inside a RecDqdx object
  const auto& dQdxInfos = objectMaps.trackPidHandler.getDqdxValues(edm4hepElem);
  if (dQdxInfos.size() != 1) {
    std::cerr << "Could not find a unique RecDqdx object for track" << std::endl;
    return false;
  }
  ASSERT_COMPARE_VALS(lcioElem->getdEdx(), dQdxInfos[0].getDQdx().value, "dEdx in Track");
  ASSERT_COMPARE_VALS(lcioElem->getdEdxError(), dQdxInfos[0].getDQdx().error, "dEdxError in Track");

  double radius = EDM4hep2LCIOConv::getRadiusOfStateAtFirstHit(edm4hepElem).value_or(-1.0);
  double radius3D = EDM4hep2LCIOConv::getRadiusOfStateAtFirstHit(edm4hepElem, true).value_or(-1.0);
  const double radiusLCIO = lcioElem->getRadiusOfInnermostHit();
  if (std::abs(radius - radiusLCIO) > radiusLCIO / 1e6 && std::abs(radius3D - radiusLCIO) > radiusLCIO / 1e6) {
    std::cerr << "radiusOfInnermostHit in Track (LCIO: " << radiusLCIO << "), EDM4hep: 2d: " << radius
              << ", 3d: " << radius3D << ")" << std::endl;
  }

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getTracks, objectMaps.tracks, "Tracks in Track");

  const auto& lcioTrackStates = lcioElem->getTrackStates();
  const auto& edm4hepTrackStates = edm4hepElem.getTrackStates();
  ASSERT_COMPARE_VALS(lcioTrackStates.size(), edm4hepTrackStates.size(), "number of TrackStates in Track");
  for (size_t i = 0; i < lcioTrackStates.size(); ++i) {
    if (!compare(lcioTrackStates[i], edm4hepTrackStates[i])) {
      std::cerr << "TrackState " << i << " in Track" << std::endl;
      return false;
    }
  }

  const auto& lcioHits = lcioElem->getTrackerHits();
  const auto edmHits = edm4hepElem.getTrackerHits();
  ASSERT_COMPARE_VALS(lcioHits.size(), edmHits.size(), "number of tracker hits in Track");
  int iHit = 0;
  for (const auto* lcioHit : lcioElem->getTrackerHits()) {
    if (const auto typedHit = dynamic_cast<const EVENT::TrackerHitPlane*>(lcioHit)) {
      if (!compareRelation(typedHit, edmHits[iHit], objectMaps.trackerHitPlanes,
                           "TrackerHit " + std::to_string(iHit) + " in Track")) {
        return false;
      }
      iHit++;
    } else if (dynamic_cast<const IMPL::TrackerHitImpl*>(lcioHit)) {
      if (!compareRelation(lcioHit, edmHits[iHit], objectMaps.trackerHits,
                           "TrackerHit " + std::to_string(iHit) + " in Track")) {
        return false;
      }
      iHit++;
    }
  }

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::Track>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TrackerHit ================

bool compare(const EVENT::TrackerHit* lcioElem, const edm4hep::TrackerHit3D& edm4hepElem, const ObjectMappings&) {
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in TrackerHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in TrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getQuality, "quality in TrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in TrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEDep, "eDep in TrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEDepError, "eDepError in TrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in TrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCovMatrix, "covMatrix in TrackerHit");
  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackerHit3DCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::TrackerHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TrackerHitPlane ================

bool compare(const EVENT::TrackerHitPlane* lcioElem, const edm4hep::TrackerHitPlane& edm4hepElem,
             const ObjectMappings&) {
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in TrackerHitPlane");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getQuality, "quality in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEDep, "eDep in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEDepError, "eDepError in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getU, "u in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getV, "v in TrackerHitPlane");
  // LCIO has getdU and getdV instead of getDu and getDv
  ASSERT_COMPARE_VALS(lcioElem->getdV(), edm4hepElem.getDv(), "dv in TrackerHitPlane");
  ASSERT_COMPARE_VALS(lcioElem->getdU(), edm4hepElem.getDu(), "du in TrackerHitPlane");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in TrackerHitPlane");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCovMatrix, "covMatrix in TrackerHitPlane");
  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::TrackerHitPlaneCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::TrackerHitPlane>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= Vertex ================

bool compare(const EVENT::Vertex* lcioElem, const edm4hep::Vertex& edm4hepElem, const ObjectMappings& objectMaps) {
  ASSERT_COMPARE(lcioElem, edm4hepElem, isPrimary, "primary in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getChi2, "chi2 in Vertex");
  ASSERT_COMPARE_VALS_FLOAT(lcioElem->getProbability(), TMath::Prob(edm4hepElem.getChi2(), edm4hepElem.getNdf()), 1e-6,
                            "probability in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCovMatrix, "covMatrix in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getParameters, "parameters in Vertex");
  // TODO: LCIO with std::string vs. EDM4hep with int
  // ASSERT_COMPARE(lcioElem, edm4hepElem, getAlgorithmType,
  //                "algorithmType in Vertex");

  // Vertices work conceptually different in EDM4hep and LCIO. Here we compare
  // the decay particles only, since that is straight forward
  const auto lcioParticles = lcioElem->getAssociatedParticle()->getParticles();
  const auto edmParticles = edm4hepElem.getParticles();
  ASSERT_COMPARE_VALS(lcioParticles.size(), edmParticles.size(), "number of associated / decay particles in Vertex");
  for (size_t iV = 0; iV < lcioParticles.size(); ++iV) {
    if (!compareRelation(lcioParticles[iV], edmParticles[iV], objectMaps.recoParticles,
                         "particle " + std::to_string(iV) + " in Vertex")) {
      return false;
    }
  }

  return true;
}

bool compare(const lcio::LCCollection* lcioCollection, const edm4hep::VertexCollection& edm4hepCollection,
             const ObjectMappings& objectMaps) {
  return compareCollection<EVENT::Vertex>(lcioCollection, edm4hepCollection, objectMaps);
}

bool compareEventHeader(const EVENT::LCEvent* lcevt, const podio::Frame* edmEvent) {
  const auto& edmEventHeader = edmEvent->get<edm4hep::EventHeaderCollection>("EventHeader")[0];

  ASSERT_COMPARE(lcevt, edmEventHeader, getEventNumber, "Event Number is not the same");
  ASSERT_COMPARE(lcevt, edmEventHeader, getRunNumber, "Run Number is not the same");
  const uint64_t lcioTimeStamp = lcevt->getTimeStamp();
  ASSERT_COMPARE_VALS(lcioTimeStamp, edmEventHeader.getTimeStamp(), "TimeStamp in EventHeader is not the same");
  ASSERT_COMPARE(lcevt, edmEventHeader, getWeight, "Weight in EventHeader is not the same");

  return true;
}

// bool compareStartVertexRelations(const edm4hep::VertexRecoParticleLinkCollection& startVtxLinks, const
// ObjectMappings& objectMaps, const podio::Frame& event) {
//   for (const auto& link : startVtxLinks) {
//     const auto edmVtx = link.getVertex();
//     const auto edmReco = link.getRec();

//   }
// }

bool compareStartVertexRelations(const EVENT::ReconstructedParticle* lcioReco,
                                 const edm4hep::VertexRecoParticleLink& link, const ObjectMappings& objectMaps) {
  const auto lcioVertex = lcioReco->getStartVertex();
  const auto edm4hepVertex = link.getFrom();
  if (!compareRelation(lcioVertex, edm4hepVertex, objectMaps.vertices, "")) {
    return false;
  }

  return true;
}

bool compareVertexRecoLink(const EVENT::Vertex*, const edm4hep::VertexRecoParticleLink&, const ObjectMappings&) {
  // TODO: Actually implement the checks
  // TODO: Figure out if this is even the right interface here
  return false;
}
