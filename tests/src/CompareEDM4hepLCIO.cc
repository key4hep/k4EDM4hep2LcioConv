#include "CompareEDM4hepLCIO.h"
#include "ComparisonUtils.h"

#include "IMPL/TrackerHitImpl.h"

#include <cstdint>

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
 *
 * TODO: Also compare relations
 */

/// Convert the two 32 bit cellIDs into one 64 bit value
template<typename LcioT>
auto to64BitCellID(LcioT* obj)
{
  const auto cellID0 = obj->getCellID0();
  const auto cellID1 = obj->getCellID1();
  uint64_t cellID = cellID1;
  cellID = (cellID << 32) | cellID0;
  return cellID;
}

// ================= CalorimeterHit ================

bool compare(
  const EVENT::CalorimeterHit* lcioElem,
  const edm4hep::CalorimeterHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in CalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergyError, "energyError in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in CalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in CalorimeterHit");
  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::CalorimeterHitCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::CalorimeterHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= Cluster ================

bool compare(const EVENT::Cluster* lcioElem, const edm4hep::Cluster& edm4hepElem, const ObjectMappings& objectMaps)
{
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
  if (!compareRelation(
        lcioElem->getCalorimeterHits(), edm4hepElem.getHits(), objectMaps.caloHits, "calorimeter hits in Cluster")) {
    return false;
  }

  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::ClusterCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::Cluster>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= MCParticle ================

bool compare(
  const EVENT::MCParticle* lcioElem,
  const edm4hep::MCParticle& edm4hepElem,
  const ObjectMappings& objectMaps)
{
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
  ASSERT_COMPARE(lcioElem, edm4hepElem, getColorFlow, "colorFlow in MCParticle");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getDaughters, objectMaps.mcParticles, "daughters in MCParticle");
  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getParents, objectMaps.mcParticles, "parents in MCParticle");

  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::MCParticleCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::MCParticle>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= ParticleID ================

bool compare(const EVENT::ParticleID* lcioElem, const edm4hep::ParticleID& edm4hepElem)
{
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in ParticleID");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPDG, "PDG in ParticleID");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getAlgorithmType, "algorithmType in ParticleID");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getLikelihood, "likelihood in ParticleID");
  return true;
}

// ================= RawCalorimeterHit ================

bool compare(
  const EVENT::RawCalorimeterHit* lcioElem,
  const edm4hep::RawCalorimeterHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in RawCalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getAmplitude, "amplitude in RawCalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTimeStamp, "timeStamp in RawCalorimeterHit");
  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::RawCalorimeterHitCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::RawCalorimeterHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= ReconstructedParticle ================

bool compare(
  const EVENT::ReconstructedParticle* lcioElem,
  const edm4hep::ReconstructedParticle& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMomentum, "momentum in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getReferencePoint, "referencePoint in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCharge, "charge in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMass, "mass in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getGoodnessOfPID, "goodnessOfPID in ReconstructedParticle");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCovMatrix, "covMatrix in ReconstructedParticle");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getClusters, objectMaps.clusters, "clusters in ReonstructedParticle");
  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getTracks, objectMaps.tracks, "tracks in ReonstructedParticle");
  ASSERT_COMPARE_RELATION(
    lcioElem, edm4hepElem, getParticles, objectMaps.recoParticles, "particles in ReonstructedParticle");
  ASSERT_COMPARE_RELATION(
    lcioElem, edm4hepElem, getStartVertex, objectMaps.vertices, "startVertex in ReconstructedParticle");

  const auto& lcioPIDs = lcioElem->getParticleIDs();
  const auto edmPIDs = edm4hepElem.getParticleIDs();
  ASSERT_COMPARE_VALS(lcioPIDs.size(), edmPIDs.size(), "particleIDs with different sizes in ReconstructedParticle");

  for (size_t i = 0; i < lcioPIDs.size(); ++i) {
    if (!compare(lcioPIDs[i], edmPIDs[i])) {
      std::cerr << "particle ID " << i << " differs in ReconstructedParticle (LCIO: " << lcioPIDs[i]
                << ", EDM4hep: " << edmPIDs[i] << ")" << std::endl;
      return false;
    }
  }

  const auto lcioPIDUsed = lcioElem->getParticleIDUsed();
  const auto edmPIDUsed = edm4hepElem.getParticleIDUsed();
  if (lcioPIDUsed == nullptr) {
    if (edmPIDUsed.isAvailable()) {
      std::cerr << "particleIDUsed is not available in LCIO, but points to " << edmPIDUsed.getObjectID()
                << " in EDM4hep for ReconstructedParticle" << std::endl;
      return false;
    }
  }
  else {
    if (!compare(lcioPIDUsed, edmPIDUsed)) {
      std::cerr << "particleIDUsed differs in ReconstructedParticle (LCIO: " << lcioPIDUsed
                << ", EDM4hep: " << edmPIDUsed << ")" << std::endl;
      return false;
    }
  }
  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::ReconstructedParticleCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::ReconstructedParticle>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= SimCalorimeterHit ================

bool compare(
  const EVENT::SimCalorimeterHit* lcioElem,
  const edm4hep::SimCalorimeterHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  const auto lcioCellID = to64BitCellID(lcioElem);
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in SimCalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in SimCalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in SimCalorimeterHit");

  // Contributions are not part of the "proper LCIO"
  const auto edmContributions = edm4hepElem.getContributions();
  ASSERT_COMPARE_VALS(lcioElem->getNMCContributions(), edmContributions.size(), "number of CaloHitContributions");

  for (int iCont = 0; iCont < lcioElem->getNMCContributions(); ++iCont) {
    const auto& edmContrib = edmContributions[iCont];
    ASSERT_COMPARE_VALS(
      lcioElem->getEnergyCont(iCont), edmContrib.getEnergy(), "energy in CaloHitContribution " + std::to_string(iCont));
    ASSERT_COMPARE_VALS(
      lcioElem->getStepPosition(iCont),
      edmContrib.getStepPosition(),
      "stepPosition in CaloHitContribution " + std::to_string(iCont));
    ASSERT_COMPARE_VALS(
      lcioElem->getTimeCont(iCont), edmContrib.getTime(), "time in CaloHitContribution " + std::to_string(iCont));

    if (!compareRelation(
          lcioElem->getParticleCont(iCont),
          edmContrib.getParticle(),
          objectMaps.mcParticles,
          " MCParticle in CaloHitContribution " + std::to_string(iCont))) {
      return false;
    }
  }

  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::SimCalorimeterHitCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::SimCalorimeterHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= SimTrackerHit ================

bool compare(
  const EVENT::SimTrackerHit* lcioElem,
  const edm4hep::SimTrackerHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCellID, "cellID in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getEDep, "EDep in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPathLength, "pathLength in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getQuality, "quality in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in SimTrackerHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getMomentum, "momentum in SimTrackerHit");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getMCParticle, objectMaps.mcParticles, "MCParticle in SimTrackerHit");

  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::SimTrackerHitCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::SimTrackerHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TPCHit ================

bool compare(const EVENT::TPCHit* lcioElem, const edm4hep::RawTimeSeries& edm4hepElem, const ObjectMappings& objectMaps)
{
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCellID, "cellID in TPCHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getQuality, "quality in TPCHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getTime, "time in TPCHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCharge, "charge in TPCHit");
  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::RawTimeSeriesCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::TPCHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TrackState ================

bool compare(const EVENT::TrackState* lcio, const edm4hep::TrackState& edm4hep)
{
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

bool compare(const EVENT::Track* lcioElem, const edm4hep::Track& edm4hepElem, const ObjectMappings& objectMaps)
{
  ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in Track");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getChi2, "chi2 in Track");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getNdf, "ndf in Track");
  // LCIO has getdEdx instead of getDEdx
  ASSERT_COMPARE_VALS(lcioElem->getdEdx(), edm4hepElem.getDEdx(), "dEdx in Track");
  ASSERT_COMPARE_VALS(lcioElem->getdEdxError(), edm4hepElem.getDEdxError(), "dEdxError in Track");
  // Also check whether these have been corretly put into the dQQuantities
  const auto dxQuantities = edm4hepElem.getDxQuantities();
  if (dxQuantities.size() != 1) {
    std::cerr << "DxQuantities have not been filled correctly, expected exactly 1, got " << dxQuantities.size()
              << " in Track" << std::endl;
    return false;
  }
  ASSERT_COMPARE_VALS(lcioElem->getdEdx(), dxQuantities[0].value, "dEdx in DxQuantities in Track");
  ASSERT_COMPARE_VALS(lcioElem->getdEdxError(), dxQuantities[0].error, "dEdxError in DxQuantities in Track");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getRadiusOfInnermostHit, "radiusOfInnermostHit in Track");

  ASSERT_COMPARE_RELATION(lcioElem, edm4hepElem, getTracks, objectMaps.tracks, "Tracks in Track");

  const auto& lcioTrackStates = lcioElem->getTrackStates();
  const auto& edm4hepTrackStates = edm4hepElem.getTrackStates();
  ASSERT_COMPARE_VALS(lcioTrackStates.size(), edm4hepTrackStates.size(), "number of TrackStates in Track");
  for (size_t i = 0; i < lcioTrackStates.size(); ++i) {
    if (!compare(lcioTrackStates[i], edm4hepTrackStates[i])) {
      std::cerr << " " << i << " in Track" << std::endl;
      return false;
    }
  }

  const auto& lcioHits = lcioElem->getTrackerHits();
  const auto edmHits = edm4hepElem.getTrackerHits();
  ASSERT_COMPARE_VALS(lcioHits.size(), edmHits.size(), "number of tracker hits in Track");
  int iHit = 0;
  for (const auto* lcioHit : lcioElem->getTrackerHits()) {
    if (const auto typedHit = dynamic_cast<const EVENT::TrackerHitPlane*>(lcioHit)) {
      if (!compareRelation(
            typedHit, edmHits[iHit], objectMaps.trackerHitPlanes, "TrackerHit " + std::to_string(iHit) + " in Track")) {
        return false;
      }
      iHit++;
    }
    else if (dynamic_cast<const IMPL::TrackerHitImpl*>(lcioHit)) {
      if (!compareRelation(
            lcioHit, edmHits[iHit], objectMaps.trackerHits, "TrackerHit " + std::to_string(iHit) + " in Track")) {
        return false;
      }
      iHit++;
    }
  }

  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::TrackCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::Track>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TrackerHit ================

bool compare(
  const EVENT::TrackerHit* lcioElem,
  const edm4hep::TrackerHit3D& edm4hepElem,
  const ObjectMappings& objectMaps)
{
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

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::TrackerHit3DCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::TrackerHit>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= TrackerHitPlane ================

bool compare(
  const EVENT::TrackerHitPlane* lcioElem,
  const edm4hep::TrackerHitPlane& edm4hepElem,
  const ObjectMappings& objectMaps)
{
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

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::TrackerHitPlaneCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::TrackerHitPlane>(lcioCollection, edm4hepCollection, objectMaps);
}

// ================= Vertex ================

bool compare(const EVENT::Vertex* lcioElem, const edm4hep::Vertex& edm4hepElem, const ObjectMappings& objectMaps)
{
  // LCIO has isPrimary (bool), EDM4hep has getPrimary (int32_t)
  ASSERT_COMPARE_VALS(lcioElem->isPrimary(), edm4hepElem.getPrimary(), "primary in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getChi2, "chi2 in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getProbability, "probability in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getCovMatrix, "covMatrix in Vertex");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getParameters, "parameters in Vertex");
  // TODO: LCIO with std::string vs. EDM4hep with int
  // ASSERT_COMPARE(lcioElem, edm4hepElem, getAlgorithmType,
  //                "algorithmType in Vertex");

  ASSERT_COMPARE_RELATION(
    lcioElem, edm4hepElem, getAssociatedParticle, objectMaps.recoParticles, "associatedParticle in Vertex");

  return true;
}

bool compare(
  const lcio::LCCollection* lcioCollection,
  const edm4hep::VertexCollection& edm4hepCollection,
  const ObjectMappings& objectMaps)
{
  return compareCollection<EVENT::Vertex>(lcioCollection, edm4hepCollection, objectMaps);
}

bool compareEventHeader(const EVENT::LCEvent* lcevt, const podio::Frame* edmEvent)
{

  const auto& edmEventHeader = edmEvent->get<edm4hep::EventHeaderCollection>("EventHeader")[0];

  ASSERT_COMPARE(lcevt, edmEventHeader, getEventNumber, "Event Number is not the same");
  ASSERT_COMPARE(lcevt, edmEventHeader, getRunNumber, "Run Number is not the same");
  ASSERT_COMPARE(lcevt, edmEventHeader, getTimeStamp, "TimeStamp in EventHeader is not the same");
  ASSERT_COMPARE(lcevt, edmEventHeader, getWeight, "Weight in EventHeader is not the same");

  return true;
}
