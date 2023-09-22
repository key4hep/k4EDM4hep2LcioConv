#include "CompareEDM4hepLCIO.h"
#include "ComparisonUtils.h"

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

// ================= CalorimeterHit ================

bool compare(
  const EVENT::CalorimeterHit* lcioElem,
  const edm4hep::CalorimeterHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  // LCIO has CellID0 and CellID1
  const auto lcioCellID = [&]() {
    const auto cellID0 = lcioElem->getCellID0();
    const auto cellID1 = lcioElem->getCellID1();

    uint64_t cellID = cellID1;
    cellID = (cellID << 32) | cellID0;
    return cellID;
  }();
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

// bool compare(const EVENT::ParticleID * lcioElem, const edm4hep::ParticleID &
// edm4hepElem) {
//     ASSERT_COMPARE(lcioElem, edm4hepElem, getType, "type in ParticleID");
//     ASSERT_COMPARE(lcioElem, edm4hepElem, getPDG, "PDG in ParticleID");
//     ASSERT_COMPARE(lcioElem, edm4hepElem, getAlgorithmType, "algorithmType in
//     ParticleID"); ASSERT_COMPARE(lcioElem, edm4hepElem, getLikelihood,
//     "likelihood in ParticleID"); return true;
// }

// ================= RawCalorimeterHit ================

bool compare(
  const EVENT::RawCalorimeterHit* lcioElem,
  const edm4hep::RawCalorimeterHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  // LCIO has CellID0 and CellID1
  const auto lcioCellID = [&]() {
    const auto cellID0 = lcioElem->getCellID0();
    const auto cellID1 = lcioElem->getCellID1();

    uint64_t cellID = cellID1;
    cellID = (cellID << 32) | cellID0;
    return cellID;
  }();
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
  // LCIO has CellID0 and CellID1
  const auto lcioCellID = [&]() {
    const auto cellID0 = lcioElem->getCellID0();
    const auto cellID1 = lcioElem->getCellID1();
    uint64_t cellID = cellID1;
    cellID = (cellID << 32) | cellID0;
    return cellID;
  }();
  ASSERT_COMPARE_VALS(lcioCellID, edm4hepElem.getCellID(), "cellID in SimCalorimeterHit");

  ASSERT_COMPARE(lcioElem, edm4hepElem, getEnergy, "energy in SimCalorimeterHit");
  ASSERT_COMPARE(lcioElem, edm4hepElem, getPosition, "position in SimCalorimeterHit");
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
  const edm4hep::TrackerHit& edm4hepElem,
  const ObjectMappings& objectMaps)
{
  // LCIO has CellID0 and CellID1
  const auto lcioCellID = [&]() {
    const auto cellID0 = lcioElem->getCellID0();
    const auto cellID1 = lcioElem->getCellID1();
    uint64_t cellID = cellID1;
    cellID = (cellID << 32) | cellID0;
    return cellID;
  }();
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
  const edm4hep::TrackerHitCollection& edm4hepCollection,
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
  // LCIO has CellID0 and CellID1
  const auto lcioCellID = [&]() {
    const auto cellID0 = lcioElem->getCellID0();
    const auto cellID1 = lcioElem->getCellID1();
    uint64_t cellID = cellID1;
    cellID = (cellID << 32) | cellID0;
    return cellID;
  }();
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
  // TODO: LCIO with std::string vs. EDM4hep with int
  // ASSERT_COMPARE(lcioElem, edm4hepElem, getAlgorithmType,
  //                "algorithmType in Vertex");
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
