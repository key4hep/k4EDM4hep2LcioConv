#include "k4EDM4hep2LcioConv/MappingUtils.h"

#include <UTIL/PIDHandler.h>

#include <edm4hep/ParticleIDCollection.h>
#include <edm4hep/RecDqdxCollection.h>

namespace LCIO2EDM4hepConv {
template <typename LCIOType>
void convertObjectParameters(LCIOType* lcioobj, podio::Frame& event) {
  const auto& params = lcioobj->getParameters();
  // handle srting params
  EVENT::StringVec keys;
  const auto stringKeys = params.getStringKeys(keys);
  for (auto i = 0u; i < stringKeys.size(); i++) {
    EVENT::StringVec sValues;
    const auto stringVals = params.getStringVals(stringKeys[i], sValues);
    event.putParameter(stringKeys[i], stringVals);
  }
  // handle float params
  EVENT::StringVec fkeys;
  const auto floatKeys = params.getFloatKeys(fkeys);
  for (auto i = 0u; i < floatKeys.size(); i++) {
    EVENT::FloatVec fValues;
    const auto floatVals = params.getFloatVals(floatKeys[i], fValues);
    event.putParameter(floatKeys[i], floatVals);
  }
  // handle int params
  EVENT::StringVec ikeys;
  const auto intKeys = params.getIntKeys(ikeys);
  for (auto i = 0u; i < intKeys.size(); i++) {
    EVENT::IntVec iValues;
    const auto intVals = params.getIntVals(intKeys[i], iValues);
    event.putParameter(intKeys[i], intVals);
  }
  // handle double params
  EVENT::StringVec dkeys;
  const auto dKeys = params.getDoubleKeys(dkeys);
  for (auto i = 0u; i < dKeys.size(); i++) {
    EVENT::DoubleVec dValues;
    const auto dVals = params.getDoubleVals(dKeys[i], dValues);
    event.putParameter(dKeys[i], dVals);
  }
}

template <typename LCVecType>
std::vector<CollNamePair> convertLCVec(const std::string& name, EVENT::LCCollection* LCCollection) {
  auto dest = std::make_unique<podio::UserDataCollection<typename LCVecType::value_type>>();
  auto vecSizes = std::make_unique<podio::UserDataCollection<uint32_t>>();
  if (LCCollection->getNumberOfElements() > 0) {
    vecSizes->push_back(0);
  }
  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    const auto* rval = static_cast<LCVecType*>(LCCollection->getElementAt(i));
    for (unsigned j = 0; j < rval->size(); j++) {
      dest->push_back((*rval)[j]);
    }
    vecSizes->push_back(dest->size());
  }
  std::vector<CollNamePair> results;
  results.reserve(2);
  results.emplace_back(name, std::move(dest));
  results.emplace_back(name + "_VecLenghts", std::move(vecSizes));
  return results;
}

template <typename CollT, typename ObjectMapT, typename LcioT, typename Edm4hepT>
auto handleSubsetColl(EVENT::LCCollection* lcioColl, const ObjectMapT& elemMap) {
  auto edm4hepColl = std::make_unique<CollT>();
  edm4hepColl->setSubsetCollection();

  UTIL::LCIterator<LcioT> lcioIter(lcioColl);
  while (const auto lcioElem = lcioIter.next()) {
    if (auto edm4hepElem = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioElem, elemMap)) {
      edm4hepColl->push_back(edm4hepElem.value());
    } else {
      std::cerr << "Cannot find corresponding EDM4hep object for an LCIO "
                   "object in a subset collection of type "
                << Edm4hepT::collection_type::valueTypeName << std::endl;
    }
  }

  return edm4hepColl;
}

template <typename MCParticleMapT>
std::unique_ptr<edm4hep::MCParticleCollection>
convertMCParticles(const std::string& name, EVENT::LCCollection* LCCollection, MCParticleMapT& mcparticlesMap) {
  auto dest = std::make_unique<edm4hep::MCParticleCollection>();
  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::MCParticle*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    lval.setPDG(rval->getPDG());
    lval.setGeneratorStatus(rval->getGeneratorStatus());

    // Convert LCIO simulator status to EDM4hep simulator status bit by
    // bit to avoid issues with different ordering
    lval.setCreatedInSimulation(rval->isCreatedInSimulation());
    lval.setBackscatter(rval->isBackscatter());
    lval.setVertexIsNotEndpointOfParent(rval->vertexIsNotEndpointOfParent());
    lval.setDecayedInTracker(rval->isDecayedInTracker());
    lval.setDecayedInCalorimeter(rval->isDecayedInCalorimeter());
    lval.setHasLeftDetector(rval->hasLeftDetector());
    lval.setStopped(rval->isStopped());
    lval.setOverlay(rval->isOverlay());

    lval.setCharge(rval->getCharge());
    lval.setTime(rval->getTime());
    lval.setMass(rval->getMass());
#ifdef EDM4HEP_MCPARTICLE_HAS_HELICITY
    lval.setHelicity(rval->getSpin()[2]);
#else
    lval.setSpin(edm4hep::Vector3f(rval->getSpin()));
#endif
    lval.setVertex(edm4hep::Vector3d(rval->getVertex()));
    lval.setEndpoint(edm4hep::Vector3d(rval->getEndpoint()));
    lval.setMomentum(rval->getMomentum());
    lval.setMomentumAtEndpoint(rval->getMomentumAtEndpoint());

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, mcparticlesMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element" << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }
  return dest;
}

template <typename RecoMapT>
std::vector<CollNamePair> convertReconstructedParticles(const std::string& name, EVENT::LCCollection* LCCollection,
                                                        RecoMapT& recoparticlesMap) {
  auto dest = std::make_unique<edm4hep::ReconstructedParticleCollection>();
  auto startVertexLinks = std::make_unique<edm4hep::VertexRecoParticleLinkCollection>();

  // Set up a PIDHandler to split off the ParticlID objects stored in the
  // reconstructed particles into separate collections. Each algorithm id /
  // name get's a separate collection
  auto pidHandler = UTIL::PIDHandler(LCCollection);
  // TODO: parameter names
  std::map<int, std::unique_ptr<edm4hep::ParticleIDCollection>> particleIDs;
  for (const auto id : pidHandler.getAlgorithmIDs()) {
    particleIDs.emplace(id, std::make_unique<edm4hep::ParticleIDCollection>());
  }

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::ReconstructedParticle*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    lval.setCharge(rval->getCharge());
    auto& m = rval->getCovMatrix(); // 10 parameters
    lval.setCovMatrix({m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9]});
    lval.setEnergy(rval->getEnergy());
    lval.setGoodnessOfPID(rval->getGoodnessOfPID());
    lval.setMass(rval->getMass());
    lval.setMomentum(Vector3fFrom(rval->getMomentum()));
    lval.setReferencePoint(rval->getReferencePoint());
    lval.setPDG(rval->getType());

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, recoparticlesMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }

    // Need to convert the particle IDs here, since they are part of the reco
    // particle collection in LCIO.
    for (const auto lcioPid : rval->getParticleIDs()) {
      auto pid = convertParticleID(lcioPid);
      pid.setParticle(lval);
      if (auto pidIt = particleIDs.find(pid.getAlgorithmType()); pidIt != particleIDs.end()) {
        pidIt->second->push_back(pid);
      } else {
        std::cerr << "ERROR: Found a PID object with an algorithm ID that is "
                     "not known to the PIDHandler (id = "
                  << pid.getAlgorithmType() << ")" << std::endl;
      }
    }

    // Keep track of the startVertex links
    if (rval->getStartVertex() != nullptr) {
      auto link = startVertexLinks->create();
      link.setTo(lval);
    }
  }

  std::vector<CollNamePair> results;
  results.reserve(particleIDs.size() + 2);
  results.emplace_back(name, std::move(dest));
  for (auto& [id, coll] : particleIDs) {
    results.emplace_back(getPIDCollName(name, pidHandler.getAlgorithmName(id)), std::move(coll));
  }
  results.emplace_back(name + "_startVertices", std::move(startVertexLinks));
  return results;
}

template <typename VertexMapT>
std::vector<CollNamePair> convertVertices(const std::string& name, EVENT::LCCollection* LCCollection,
                                          VertexMapT& vertexMap) {
  auto dest = std::make_unique<edm4hep::VertexCollection>();
  auto assocParticles = std::make_unique<edm4hep::VertexRecoParticleLinkCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::Vertex*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    lval.setPrimary(rval->isPrimary() ? 1 : 0); // 1 for primary and 0 for not primary
    lval.setChi2(rval->getChi2());
    lval.setNdf(find_ndf(rval->getChi2(), rval->getProbability()));
    lval.setPosition(rval->getPosition());
    auto& m = rval->getCovMatrix(); // 6 parameters
    lval.setCovMatrix({m[0], m[1], m[2], m[3], m[4], m[5]});
    // NOTE: the algorithm type in LCIO is a string, but an integer is expected
    // lval.setAlgorithmType(rval->getAlgorithmType());

    for (auto v : rval->getParameters()) {
      lval.addToParameters(v);
    }

    auto assoc = assocParticles->create();
    assoc.setFrom(lval);

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, vertexMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  std::vector<CollNamePair> results;
  results.reserve(2);
  results.emplace_back(name, std::move(dest));
  results.emplace_back(name + "_associatedParticles", std::move(assocParticles));
  return results;
}

template <typename SimTrHitMapT>
std::unique_ptr<edm4hep::SimTrackerHitCollection>
convertSimTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, SimTrHitMapT& SimTrHitMap) {
  auto dest = std::make_unique<edm4hep::SimTrackerHitCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::SimTrackerHit*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    uint64_t cellID = rval->getCellID1();
    cellID = (cellID << 32) | rval->getCellID0();
    lval.setCellID(cellID);
    lval.setEDep(rval->getEDep());
    lval.setTime(rval->getTime());
    lval.setPathLength(rval->getPathLength());
    lval.setQuality(rval->getQuality());
    lval.setPosition(rval->getPosition());
    lval.setMomentum(rval->getMomentum());

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, SimTrHitMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }
  return dest;
}

template <typename HitMapT>
std::unique_ptr<edm4hep::RawTimeSeriesCollection>
convertTPCHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TPCHitMap) {
  auto dest = std::make_unique<edm4hep::RawTimeSeriesCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::TPCHit*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    lval.setCellID(rval->getCellID());
    lval.setTime(rval->getTime());
    lval.setCharge(rval->getCharge());
    lval.setQuality(rval->getQuality());
    for (unsigned j = 0, M = rval->getNRawDataWords(); j < M; j++) {
      lval.addToAdcCounts(rval->getRawDataWord(j));
    }
    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, TPCHitMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  return dest;
}

template <typename HitMapT>
std::unique_ptr<edm4hep::TrackerHit3DCollection>
convertTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitMap) {
  auto dest = std::make_unique<edm4hep::TrackerHit3DCollection>();
  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::TrackerHit*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    uint64_t cellID = rval->getCellID1();
    cellID = (cellID << 32) | rval->getCellID0();
    lval.setCellID(cellID);
    lval.setType(rval->getType());
    lval.setQuality(rval->getQuality());
    lval.setTime(rval->getTime());
    lval.setEDep(rval->getEDep());
    lval.setEDepError(rval->getEDepError());
    lval.setPosition(rval->getPosition());
    auto& m = rval->getCovMatrix();
    lval.setCovMatrix({m[0], m[1], m[2], m[3], m[4], m[5]});

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, TrackerHitMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }
  return dest;
}

template <typename HitMapT>
std::unique_ptr<edm4hep::TrackerHitPlaneCollection>
convertTrackerHitPlanes(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitPlaneMap) {
  auto dest = std::make_unique<edm4hep::TrackerHitPlaneCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::TrackerHitPlane*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    uint64_t cellID = rval->getCellID1();
    cellID = (cellID << 32) | rval->getCellID0();
    lval.setCellID(cellID);
    lval.setType(rval->getType());
    lval.setQuality(rval->getQuality());
    lval.setTime(rval->getTime());
    lval.setEDep(rval->getEDep());
    lval.setEDepError(rval->getEDepError());
    lval.setPosition(rval->getPosition());
    lval.setU({rval->getU()[0], rval->getU()[1]});
    lval.setV({rval->getV()[0], rval->getV()[1]});
    lval.setDu(rval->getdU());
    lval.setDv(rval->getdV());
    auto& m = rval->getCovMatrix();
    lval.setCovMatrix({m[0], m[1], m[2], m[3], m[4], m[5]});

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, TrackerHitPlaneMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  return dest;
}

template <typename TrackMapT>
std::vector<CollNamePair> convertTracks(const std::string& name, EVENT::LCCollection* LCCollection,
                                        TrackMapT& TrackMap) {
  auto dest = std::make_unique<edm4hep::TrackCollection>();
  auto dQdxColl = std::make_unique<edm4hep::RecDqdxCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::Track*>(LCCollection->getElementAt(i));
    auto lval = dest->create();
    auto trackDqdx = dQdxColl->create();
    trackDqdx.setTrack(lval);

    lval.setType(rval->getType());
    lval.setChi2(rval->getChi2());
    lval.setNdf(rval->getNdf());
    lval.setNholes(rval->getNholes());

    auto& dqdx = trackDqdx.getDQdx();
    dqdx.value = rval->getdEdx();
    dqdx.error = rval->getdEdxError();

    const auto& subdetectorHitNum = rval->getSubdetectorHitNumbers();
    for (auto hitNum : subdetectorHitNum) {
      lval.addToSubdetectorHitNumbers(hitNum);
    }
    const auto& subdetectorHoleNum = rval->getSubdetectorHoleNumbers();
    for (auto holeNum : subdetectorHoleNum) {
      lval.addToSubdetectorHoleNumbers(holeNum);
    }

    const auto& trackStates = rval->getTrackStates();
    for (auto& trackState : trackStates) {
      lval.addToTrackStates(convertTrackState(trackState));
    }

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, TrackMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  std::vector<CollNamePair> results{};
  results.emplace_back(name, std::move(dest));
  results.emplace_back(name + "_dQdx", std::move(dQdxColl));
  return results;
}

template <typename HitMapT>
std::unique_ptr<edm4hep::SimCalorimeterHitCollection>
convertSimCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& SimCaloHitMap) {
  auto dest = std::make_unique<edm4hep::SimCalorimeterHitCollection>();
  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::SimCalorimeterHit*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    uint64_t cellID = rval->getCellID1();
    cellID = (cellID << 32) | rval->getCellID0();
    lval.setCellID(cellID);
    lval.setEnergy(rval->getEnergy());
    lval.setPosition(rval->getPosition());

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, SimCaloHitMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  return dest;
}

template <typename HitMapT>
std::unique_ptr<edm4hep::RawCalorimeterHitCollection>
convertRawCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& rawCaloHitMap) {
  auto dest = std::make_unique<edm4hep::RawCalorimeterHitCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::RawCalorimeterHit*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    uint64_t cellID = rval->getCellID1();
    cellID = (cellID << 32) | rval->getCellID0();
    lval.setCellID(cellID);
    lval.setAmplitude(rval->getAmplitude());
    lval.setTimeStamp(rval->getTimeStamp());

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, rawCaloHitMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  return dest;
}

template <typename HitMapT>
std::unique_ptr<edm4hep::CalorimeterHitCollection>
convertCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& caloHitMap) {
  auto dest = std::make_unique<edm4hep::CalorimeterHitCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::CalorimeterHit*>(LCCollection->getElementAt(i));
    auto lval = dest->create();
    uint64_t cellID = rval->getCellID1();
    cellID = (cellID << 32) | rval->getCellID0();
    lval.setCellID(cellID);
    lval.setEnergy(rval->getEnergy());
    lval.setEnergyError(rval->getEnergyError());
    lval.setPosition(rval->getPosition());
    lval.setTime(rval->getTime());
    lval.setType(rval->getType());

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, caloHitMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
      const auto existingId = existing.id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }

  return dest;
}

template <typename ClusterMapT>
std::unique_ptr<edm4hep::ClusterCollection> convertClusters(const std::string& name, EVENT::LCCollection* LCCollection,
                                                            ClusterMapT& clusterMap) {
  auto dest = std::make_unique<edm4hep::ClusterCollection>();

  for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
    auto* rval = static_cast<EVENT::Cluster*>(LCCollection->getElementAt(i));
    auto lval = dest->create();

    lval.setEnergy(rval->getEnergy());
    lval.setEnergyError(rval->getEnergyError());
    lval.setITheta(rval->getITheta());
    lval.setPhi(rval->getIPhi());
    lval.setPosition(rval->getPosition());
    auto& m = rval->getPositionError();
    lval.setPositionError({m[0], m[1], m[2], m[3], m[4], m[5]});
    lval.setType(rval->getType());
    lval.setDirectionError(Vector3fFrom(rval->getDirectionError()));

    const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, clusterMap);
    if (!inserted) {
      auto existing = k4EDM4hep2LcioConv::detail::getKey(iterator);
      const auto existingId = existing->id();
      std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                << " collection" << std::endl;
    }
  }
  return dest;
}

template <typename ObjectMappingT>
std::vector<CollNamePair> convertCollection(const std::string& name, EVENT::LCCollection* LCCollection,
                                            ObjectMappingT& typeMapping) {
  const auto& type = LCCollection->getTypeName();
  std::vector<CollNamePair> retColls;
  if (type == "MCParticle") {
    retColls.emplace_back(name, convertMCParticles(name, LCCollection, typeMapping.mcParticles));
  } else if (type == "ReconstructedParticle") {
    return convertReconstructedParticles(name, LCCollection, typeMapping.recoParticles);
  } else if (type == "Vertex") {
    return convertVertices(name, LCCollection, typeMapping.vertices);
  } else if (type == "Track") {
    return convertTracks(name, LCCollection, typeMapping.tracks);
  } else if (type == "Cluster") {
    retColls.emplace_back(name, convertClusters(name, LCCollection, typeMapping.clusters));
  } else if (type == "SimCalorimeterHit") {
    retColls.emplace_back(name, convertSimCalorimeterHits(name, LCCollection, typeMapping.simCaloHits));
  } else if (type == "RawCalorimeterHit") {
    retColls.emplace_back(name, convertRawCalorimeterHits(name, LCCollection, typeMapping.rawCaloHits));
  } else if (type == "CalorimeterHit") {
    retColls.emplace_back(name, convertCalorimeterHits(name, LCCollection, typeMapping.caloHits));
  } else if (type == "SimTrackerHit") {
    retColls.emplace_back(name, convertSimTrackerHits(name, LCCollection, typeMapping.simTrackerHits));
  } else if (type == "TPCHit") {
    retColls.emplace_back(name, convertTPCHits(name, LCCollection, typeMapping.tpcHits));
  } else if (type == "TrackerHit") {
    retColls.emplace_back(name, convertTrackerHits(name, LCCollection, typeMapping.trackerHits));
  } else if (type == "TrackerHitPlane") {
    retColls.emplace_back(name, convertTrackerHitPlanes(name, LCCollection, typeMapping.trackerHitPlanes));
  } else if (type == "LCIntVec") {
    return convertLCVec<EVENT::LCIntVec>(name, LCCollection);
  } else if (type == "LCFloatVec") {
    return convertLCVec<EVENT::LCFloatVec>(name, LCCollection);
  } else if (type != "LCRelation") {
    std::cerr << type << " is a collection type for which no known conversion exists." << std::endl;
  }
  return retColls;
}

template <typename HitMapT, typename MCParticleMapT>
std::unique_ptr<edm4hep::CaloHitContributionCollection>
createCaloHitContributions(HitMapT& SimCaloHitMap, const MCParticleMapT& mcparticlesMap) {
  auto contrCollection = std::make_unique<edm4hep::CaloHitContributionCollection>();
  for (auto& [lcioHit, edmHit] : SimCaloHitMap) {
    auto NMCParticle = lcioHit->getNMCParticles();
    for (int j = 0; j < NMCParticle; j++) {
      auto edm_contr = contrCollection->create();
      edmHit.addToContributions(edm_contr);

      edm_contr.setPDG(lcioHit->getPDGCont(j));
      edm_contr.setTime(lcioHit->getTimeCont(j));
      edm_contr.setEnergy(lcioHit->getEnergyCont(j));
      edm_contr.setStepPosition(lcioHit->getStepPosition(j));
      auto lcioParticle = (lcioHit->getParticleCont(j));
      if (lcioParticle != nullptr) {
        if (const auto edm4hepParticle = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioParticle, mcparticlesMap)) {
          edm_contr.setParticle(edm4hepParticle.value());
        } else {
          std::cerr << "Cannot find corresponding EDM4hep MCParticle for a "
                       "LCIO MCParticle, "
                    << "while trying to build CaloHitContributions " << std::endl;
        }
      }
    }
  }
  return contrCollection;
}

template <typename MCParticleMapT, typename MCParticleLookupMapT>
void resolveRelationsMCParticles(MCParticleMapT& mcparticlesMap, const MCParticleLookupMapT& lookupMap) {
  for (auto& [lcio, edm] : mcparticlesMap) {
    auto daughters = lcio->getDaughters();
    auto parents = lcio->getParents();

    for (auto d : daughters) {
      if (d == nullptr) {
        continue;
      }
      if (const auto edmD = k4EDM4hep2LcioConv::detail::mapLookupTo(d, lookupMap)) {
        edm.addToDaughters(edmD.value());
      } else {
        std::cerr << "Cannot find corresponding EDM4hep MCParticle for an LCIO "
                     "MCParticle, "
                     "while trying to resolve the daughters of MCParticles"
                  << std::endl;
      }
    }
    for (auto p : parents) {
      if (p == nullptr) {
        continue;
      }
      if (const auto edmP = k4EDM4hep2LcioConv::detail::mapLookupTo(p, lookupMap)) {
        edm.addToParents(edmP.value());
      } else {
        std::cerr << "Cannot find corresponding EDM4hep MCParticle for the LCIO "
                     "MCParticle, "
                     "while trying to resolve the parents of MCParticles Collections"
                  << std::endl;
      }
    }
  }
}

template <typename HitMapT, typename MCParticleMapT>
void resolveRelationsSimTrackerHits(HitMapT& SimTrHitMap, const MCParticleMapT& mcparticlesMap) {
  for (auto& [lcio, edm] : SimTrHitMap) {
    auto mcps = lcio->getMCParticle();
    if (mcps == nullptr) {
      continue;
    }
    if (const auto edmP = k4EDM4hep2LcioConv::detail::mapLookupTo(mcps, mcparticlesMap)) {
      edm.setParticle(edmP.value());
    } else {
      std::cerr << "Cannot find corresponding EDM4hep MCParticle for the LCIO "
                   "MCParticle, "
                   "while trying to resolve the SimTrackHit Relations"
                << std::endl;
    }
  }
}

template <typename RecoParticleMapT, typename RecoParticleLookupMapT, typename ClusterMapT, typename TrackMapT>
void resolveRelationsRecoParticles(RecoParticleMapT& recoparticlesMap, const RecoParticleLookupMapT& recoLookupMap,
                                   const ClusterMapT& clusterMap, const TrackMapT& tracksMap) {
  for (auto& [lcio, edm] : recoparticlesMap) {
    auto clusters = lcio->getClusters();
    for (auto c : clusters) {
      if (c == nullptr) {
        continue;
      }
      if (const auto edmC = k4EDM4hep2LcioConv::detail::mapLookupTo(c, clusterMap)) {
        edm.addToClusters(edmC.value());
      } else {
        std::cerr << "Cannot find corresponding EDM4hep Cluster for a LCIO Cluster, "
                     "while trying to resolve the ReconstructedParticle Relations"
                  << std::endl;
      }
    }

    auto tracks = lcio->getTracks();
    for (auto t : tracks) {
      if (t == nullptr) {
        continue;
      }
      if (const auto edmT = k4EDM4hep2LcioConv::detail::mapLookupTo(t, tracksMap)) {
        edm.addToTracks(edmT.value());
      } else {
        std::cerr << "Cannot find corresponding EDM4hep Tracks for a LCIO Tracks, "
                     "while trying to resolve the ReconstructedParticle Relations"
                  << std::endl;
      }
    }

    auto parents = lcio->getParticles();
    for (auto p : parents) {
      if (p == nullptr) {
        continue;
      }
      if (const auto edmReco = k4EDM4hep2LcioConv::detail::mapLookupTo(p, recoLookupMap)) {
        edm.addToParticles(edmReco.value());
      } else {
        std::cerr << "Cannot find corresponding EDM4hep RecoParticle for a "
                     "LCIO RecoParticle, "
                     "while trying to resolve the ReconstructedParticles "
                     "parents Relations"
                  << std::endl;
      }
    }
  }
}

template <typename ClusterMapT, typename CaloHitMapT>
void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap) {
  for (auto& [lcio, edm] : clustersMap) {
    auto clusters = lcio->getClusters();
    auto calohits = lcio->getCalorimeterHits();
    auto shape = lcio->getShape();
    auto subdetectorEnergies = lcio->getSubdetectorEnergies();
    for (auto c : clusters) {
      if (c == nullptr) {
        continue;
      }
      if (const auto edmC = k4EDM4hep2LcioConv::detail::mapLookupTo(c, clustersMap)) {
        edm.addToClusters(edmC.value());
      } else {
        std::cerr << "Couldn't find cluster to add to Relations in edm" << std::endl;
      }
    }
    for (auto cal : calohits) {
      if (cal == nullptr) {
        continue;
      }
      if (const auto edmCaloHit = k4EDM4hep2LcioConv::detail::mapLookupTo(cal, caloHitMap)) {
        edm.addToHits(edmCaloHit.value());
      } else {
        std::cerr << "Couldn't find CaloHit to add to Relations for Clusters in edm" << std::endl;
      }
    }
    for (auto s : shape) {
      edm.addToShapeParameters(s);
    }
    for (auto subE : subdetectorEnergies) {
      edm.addToSubdetectorEnergies(subE);
    }
  }
}

template <typename TrackMapT, typename TrackHitMapT, typename THPlaneHitMapT, typename TPCHitMapT>
void resolveRelationsTracks(TrackMapT& tracksMap, const TrackHitMapT& trackerHitMap,
                            const THPlaneHitMapT& trackerHitPlaneMap, const TPCHitMapT&) {
  for (auto& [lcio, edm] : tracksMap) {
    auto tracks = lcio->getTracks();
    for (auto t : tracks) {
      if (t == nullptr) {
        continue;
      }
      if (const auto track = k4EDM4hep2LcioConv::detail::mapLookupTo(t, tracksMap)) {
        edm.addToTracks(track.value());
      } else {
        // std::cerr << "Couldn't find tracks to add to Tracks Relations in edm"
        // << std::endl;
      }
    }
    const auto trackerHits = lcio->getTrackerHits();
    for (const auto th : trackerHits) {
      bool found = false;
      if (th == nullptr) {
        continue;
      }
      if (const auto typedTHPlane = dynamic_cast<EVENT::TrackerHitPlane*>(th)) {
        if (const auto trHit = k4EDM4hep2LcioConv::detail::mapLookupTo(typedTHPlane, trackerHitPlaneMap)) {
          edm.addToTrackerHits(trHit.value());
          found = true;
        }
      } else if (auto typedTH = dynamic_cast<EVENT::TrackerHit*>(th)) {
        if (const auto trHit = k4EDM4hep2LcioConv::detail::mapLookupTo(typedTH, trackerHitMap)) {
          edm.addToTrackerHits(trHit.value());
          found = true;
        }
      }

      if (!found) {
        std::cerr << "Couldn't find a edm4hep TrackerHit for an LCIO "
                     "TrackerHit when resolving "
                  << "relations for a Track" << std::endl;
      }
    }
  }
}

template <typename VertexMapT, typename URecoParticleMapT, typename LURecoParticleMapT>
void resolveRelationsVertices(VertexMapT& vertexMap, URecoParticleMapT& updateRPMap,
                              const LURecoParticleMapT& lookupRPMap) {
  for (auto& [lcioVtx, edmVtx] : vertexMap) {
    const auto recoparticle = lcioVtx->getAssociatedParticle();
    if (recoparticle == nullptr) {
      continue;
    }
    if (auto edmReco = k4EDM4hep2LcioConv::detail::mapLookupTo(recoparticle, updateRPMap)) {
      edmReco->setDecayVertex(edmVtx);
    } else {
      std::cerr << "Could not find a reco particle to attach a Vertex to" << std::endl;
    }

    // Attach the decay particles
    for (const auto& p : recoparticle->getParticles()) {
      if (const auto& edm_p = k4EDM4hep2LcioConv::detail::mapLookupTo(p, lookupRPMap)) {
        edmVtx.addToParticles(edm_p.value());
      } else {
        std::cerr << "Could not find a (decay) particle to add to a Vertex" << std::endl;
      }
    }
  }
}

template <typename VertexMapT, typename RecoParticleMapT>
void finalizeVertexRecoParticleLinks(edm4hep::VertexRecoParticleLinkCollection& links, const VertexMapT& vertexMap,
                                     const RecoParticleMapT& recoParticleMap) {
  for (auto link : links) {
    auto linkRec = link.getTo();
    auto linkVtx = link.getFrom();

    if (linkRec.isAvailable()) {
      // This is the link that points from the particles to their start vertex
      if (const auto lcioRec = k4EDM4hep2LcioConv::detail::mapLookupFrom(linkRec, recoParticleMap)) {
        const auto lcioStartVtx = lcioRec.value()->getStartVertex();
        if (const auto startVtx = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioStartVtx, vertexMap)) {
          link.setFrom(startVtx.value());
        } else {
          std::cerr << "Could not find start vertex while finalizing the RecoParticle - Vertex links" << std::endl;
        }
      } else {
        std::cerr << "Could not find a corresponding LCIO reco particle for finalizing the RecoParticle - Vertex links"
                  << std::endl;
      }
    } else {
      // This is the link that points from the vertex to the linkiated particle
      if (const auto lcioVtx = k4EDM4hep2LcioConv::detail::mapLookupFrom(linkVtx, vertexMap)) {
        const auto lcioLinkParticle = lcioVtx.value()->getAssociatedParticle();
        if (const auto linkParticle = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioLinkParticle, recoParticleMap)) {
          link.setTo(linkParticle.value());
        } else {
          std::cerr << "Could not find an linkiated particle while finalizing the RecoParticle - Vertex links"
                    << std::endl;
        }
      } else {
        std::cerr << "Could not find a corresponding LCIO vertex for finalizing the RecoParticle - Vertex links"
                  << std::endl;
      }
    }
  }
}

template <typename ObjectMappingT>
void resolveRelations(ObjectMappingT& typeMapping) {
  resolveRelations(typeMapping, typeMapping);
}

template <typename ObjectMappingT, typename ObjectMappingU>
void resolveRelations(ObjectMappingT& updateMaps, const ObjectMappingU& lookupMaps) {
  resolveRelationsMCParticles(updateMaps.mcParticles, lookupMaps.mcParticles);
  resolveRelationsRecoParticles(updateMaps.recoParticles, lookupMaps.recoParticles, lookupMaps.clusters,
                                lookupMaps.tracks);
  resolveRelationsSimTrackerHits(updateMaps.simTrackerHits, lookupMaps.mcParticles);
  resolveRelationsClusters(updateMaps.clusters, lookupMaps.caloHits);
  resolveRelationsTracks(updateMaps.tracks, lookupMaps.trackerHits, lookupMaps.trackerHitPlanes, lookupMaps.tpcHits);
  resolveRelationsVertices(updateMaps.vertices, updateMaps.recoParticles, lookupMaps.recoParticles);
}

template <typename CollT, bool Reverse, typename FromMapT, typename ToMapT, typename FromLCIOT, typename ToLCIOT,
          typename FromEDM4hepT, typename ToEDM4hepT>
std::unique_ptr<CollT> createLinkCollection(EVENT::LCCollection* relations, const FromMapT& fromMap,
                                            const ToMapT& toMap) {
  auto linkColl = std::make_unique<CollT>();
  auto relIter = UTIL::LCIterator<EVENT::LCRelation>(relations);

  while (const auto rel = relIter.next()) {
    auto link = linkColl->create();
    link.setWeight(rel->getWeight());
    const auto lcioTo = static_cast<ToLCIOT*>(rel->getTo());
    const auto lcioFrom = static_cast<FromLCIOT*>(rel->getFrom());
    const auto edm4hepTo = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioTo, toMap);
    const auto edm4hepFrom = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioFrom, fromMap);
    if (edm4hepTo.has_value() && edm4hepFrom.has_value()) {
      if constexpr (Reverse) {
        if constexpr (std::is_same_v<k4EDM4hep2LcioConv::detail::mutable_t<ToEDM4hepT>, edm4hep::MutableVertex>) {
          link.setFrom(*edm4hepTo);
          link.setTo(*edm4hepFrom);
        } else {
          link.setTo(*edm4hepTo);
          link.setFrom(*edm4hepFrom);
        }
      } else {
        if constexpr (std::is_same_v<k4EDM4hep2LcioConv::detail::mutable_t<FromEDM4hepT>, edm4hep::MutableVertex>) {
          link.setFrom(*edm4hepFrom);
          link.setTo(*edm4hepTo);
        } else {
          link.setTo(*edm4hepFrom);
          link.setFrom(*edm4hepTo);
        }
      }
    }
  }

  return linkColl;
}

template <typename SimHitMap, typename Hit3DMap, typename HitPlaneMap>
std::unique_ptr<edm4hep::TrackerHitSimTrackerHitLinkCollection>
createLinkCollection(EVENT::LCCollection* relations, const SimHitMap& simHitMap, const Hit3DMap& hit3DMap,
                     const HitPlaneMap& hitPlaneMap, bool reverse) {
  auto linkColl = std::make_unique<edm4hep::TrackerHitSimTrackerHitLinkCollection>();
  auto relIter = UTIL::LCIterator<EVENT::LCRelation>(relations);

  while (const auto rel = relIter.next()) {
    auto link = linkColl->create();
    link.setWeight(rel->getWeight());

    const auto lcioSimHit = [&]() {
      if (reverse) {
        return static_cast<EVENT::SimTrackerHit*>(rel->getFrom());
      } else {
        return static_cast<EVENT::SimTrackerHit*>(rel->getTo());
      }
    }();

    const auto lcioHit = [&]() {
      if (reverse) {
        return static_cast<EVENT::TrackerHit*>(rel->getTo());
      } else {
        return static_cast<EVENT::TrackerHit*>(rel->getFrom());
      }
    }();

    const auto edm4hepSimHit = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioSimHit, simHitMap);
    if (edm4hepSimHit) {
      link.setTo(edm4hepSimHit.value());
    }

    // Have to look in both maps here to be sure
    const auto hit3D = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioHit, hit3DMap);
    if (hit3D) {
      link.setFrom(hit3D.value());
      continue; // skip a lookup that will fail in any case
    }
    const auto hitPlane =
        k4EDM4hep2LcioConv::detail::mapLookupTo(static_cast<EVENT::TrackerHitPlane*>(lcioHit), hitPlaneMap);
    if (hitPlane) {
      link.setFrom(hitPlane.value());
    }
  }
  return linkColl;
}

template <typename ObjectMappingT>
std::vector<CollNamePair> createLinks(const ObjectMappingT& typeMapping,
                                      const std::vector<std::pair<std::string, EVENT::LCCollection*>>& LCRelation) {
  std::vector<CollNamePair> linksCollVec;
  for (const auto& [name, relations] : LCRelation) {
    const auto& params = relations->getParameters();

    const auto& fromType = params.getStringVal("FromType");
    const auto& toType = params.getStringVal("ToType");
    if (fromType.empty() || toType.empty()) {
      std::cerr << "LCRelation collection " << name << " has missing FromType or ToType parameters. "
                << "Cannot convert it without this information." << std::endl;
      continue;
    }

    if (fromType == "MCParticle" && toType == "ReconstructedParticle") {
      auto mc_a = createLinkCollection<edm4hep::RecoMCParticleLinkCollection, false>(relations, typeMapping.mcParticles,
                                                                                     typeMapping.recoParticles);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "ReconstructedParticle" && toType == "MCParticle") {
      auto mc_a = createLinkCollection<edm4hep::RecoMCParticleLinkCollection, true>(
          relations, typeMapping.recoParticles, typeMapping.mcParticles);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "CalorimeterHit" && toType == "SimCalorimeterHit") {
      auto mc_a = createLinkCollection<edm4hep::CaloHitSimCaloHitLinkCollection, true>(relations, typeMapping.caloHits,
                                                                                       typeMapping.simCaloHits);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "SimCalorimeterHit" && toType == "CalorimeterHit") {
      auto mc_a = createLinkCollection<edm4hep::CaloHitSimCaloHitLinkCollection, false>(
          relations, typeMapping.simCaloHits, typeMapping.caloHits);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "Cluster" && toType == "MCParticle") {
      auto mc_a = createLinkCollection<edm4hep::ClusterMCParticleLinkCollection, true>(relations, typeMapping.clusters,
                                                                                       typeMapping.mcParticles);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "MCParticle" && toType == "Cluster") {
      auto mc_a = createLinkCollection<edm4hep::ClusterMCParticleLinkCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.clusters);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "MCParticle" && toType == "Track") {
      auto mc_a = createLinkCollection<edm4hep::TrackMCParticleLinkCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.tracks);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "Track" && toType == "MCParticle") {
      auto mc_a = createLinkCollection<edm4hep::TrackMCParticleLinkCollection, true>(relations, typeMapping.tracks,
                                                                                     typeMapping.mcParticles);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if ((fromType == "TrackerHit" || fromType == "TrackerHitPlane") && toType == "SimTrackerHit") {
      auto mc_a = createLinkCollection(relations, typeMapping.simTrackerHits, typeMapping.trackerHits,
                                       typeMapping.trackerHitPlanes, false);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "SimTrackerHit" && (toType == "TrackerHit" || toType == "TrackerHitPlane")) {
      auto mc_a = createLinkCollection(relations, typeMapping.simTrackerHits, typeMapping.trackerHits,
                                       typeMapping.trackerHitPlanes, true);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "ReconstructedParticle" && toType == "Vertex") {
      auto mc_a = createLinkCollection<edm4hep::VertexRecoParticleLinkCollection, true>(
          relations, typeMapping.recoParticles, typeMapping.vertices);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "Vertex" && toType == "ReconstructedParticle") {
      auto mc_a = createLinkCollection<edm4hep::VertexRecoParticleLinkCollection, false>(
          relations, typeMapping.vertices, typeMapping.recoParticles);
      linksCollVec.emplace_back(name, std::move(mc_a));
    } else if (fromType == "CalorimeterHit" && toType == "MCParticle") {
      auto assoc = createLinkCollection<edm4hep::CaloHitMCParticleLinkCollection, true>(relations, typeMapping.caloHits,
                                                                                        typeMapping.mcParticles);
      linksCollVec.emplace_back(name, std::move(assoc));
    } else if (fromType == "MCParticle" && toType == "CalorimeterHit") {
      auto assoc = createLinkCollection<edm4hep::CaloHitMCParticleLinkCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.caloHits);
      linksCollVec.emplace_back(name, std::move(assoc));
    } else {
      std::cout << "Relation from: " << fromType << " to: " << toType << " (" << name
                << ") is not beeing handled during creation of associations" << std::endl;
    }
  }

  return linksCollVec;
}

template <typename ObjectMappingT>
std::unique_ptr<podio::CollectionBase> fillSubset(EVENT::LCCollection* LCCollection, const ObjectMappingT& typeMapping,
                                                  const std::string& type) {
  if (type == "MCParticle") {
    return handleSubsetColl<edm4hep::MCParticleCollection>(LCCollection, typeMapping.mcParticles);
  } else if (type == "ReconstructedParticle") {
    return handleSubsetColl<edm4hep::ReconstructedParticleCollection>(LCCollection, typeMapping.recoParticles);
  } else if (type == "Vertex") {
    return handleSubsetColl<edm4hep::VertexCollection>(LCCollection, typeMapping.vertices);
  } else if (type == "Track") {
    return handleSubsetColl<edm4hep::TrackCollection>(LCCollection, typeMapping.tracks);
  } else if (type == "Cluster") {
    return handleSubsetColl<edm4hep::ClusterCollection>(LCCollection, typeMapping.clusters);
  } else if (type == "SimCalorimeterHit") {
    return handleSubsetColl<edm4hep::SimCalorimeterHitCollection>(LCCollection, typeMapping.simCaloHits);
  } else if (type == "RawCalorimeterHit") {
    return handleSubsetColl<edm4hep::RawCalorimeterHitCollection>(LCCollection, typeMapping.rawCaloHits);
  } else if (type == "CalorimeterHit") {
    return handleSubsetColl<edm4hep::CalorimeterHitCollection>(LCCollection, typeMapping.caloHits);
  } else if (type == "SimTrackerHit") {
    return handleSubsetColl<edm4hep::SimTrackerHitCollection>(LCCollection, typeMapping.simTrackerHits);
  } else if (type == "TPCHit") {
    return handleSubsetColl<edm4hep::RawTimeSeriesCollection>(LCCollection, typeMapping.tpcHits);
  } else if (type == "TrackerHit") {
    return handleSubsetColl<edm4hep::TrackerHit3DCollection>(LCCollection, typeMapping.trackerHits);
  } else if (type == "TrackerHitPlane") {
    return handleSubsetColl<edm4hep::TrackerHitPlaneCollection>(LCCollection, typeMapping.trackerHitPlanes);
  } else {
    return nullptr;
  }
}

} // namespace LCIO2EDM4hepConv
