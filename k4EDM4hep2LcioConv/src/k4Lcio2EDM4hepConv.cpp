#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

#include <iostream>

namespace LCIO2EDM4hepConv {
  edm4hep::TrackState convertTrackState(const EVENT::TrackState* trackState)
  {
    auto edmtrackState = edm4hep::TrackState {};
    edmtrackState.location = trackState->getLocation();
    edmtrackState.D0 = trackState->getD0();
    edmtrackState.phi = trackState->getPhi();
    edmtrackState.omega = trackState->getOmega();
    edmtrackState.Z0 = trackState->getZ0();
    edmtrackState.tanLambda = trackState->getTanLambda();
    // not available in lcio
    edmtrackState.time = -1;
    const auto refPoint = trackState->getReferencePoint();
    edmtrackState.referencePoint = Vector3fFrom({refPoint[0], refPoint[1], refPoint[2]});
    const auto& covMatrix = trackState->getCovMatrix();
    edmtrackState.covMatrix = {
      covMatrix[0],
      covMatrix[1],
      covMatrix[2],
      covMatrix[3],
      covMatrix[4],
      covMatrix[5],
      covMatrix[6],
      covMatrix[7],
      covMatrix[8],
      covMatrix[9],
      covMatrix[10],
      covMatrix[11],
      covMatrix[12],
      covMatrix[13],
      covMatrix[15],
      0,
      0,
      0,
      0,
      0,
      0};

    return edmtrackState;
  }

  std::unique_ptr<edm4hep::MCParticleCollection> convertMCParticle(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::MCParticle*, edm4hep::MutableMCParticle>& mcparticlesMap)
  {
    auto dest = std::make_unique<edm4hep::MCParticleCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::MCParticle*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setPDG(rval->getPDG());
      lval.setGeneratorStatus(rval->getGeneratorStatus());
      lval.setSimulatorStatus(rval->getSimulatorStatus());
      lval.setCharge(rval->getCharge());
      lval.setTime(rval->getTime());
      lval.setMass(rval->getMass());
      lval.setSpin(edm4hep::Vector3f(rval->getSpin()));
      lval.setColorFlow(edm4hep::Vector2i(rval->getColorFlow()));
      lval.setVertex(edm4hep::Vector3d(rval->getVertex()));
      lval.setEndpoint(edm4hep::Vector3d(rval->getEndpoint()));
      lval.setMomentum(Vector3fFrom(rval->getMomentum()));
      lval.setMomentumAtEndpoint(Vector3fFrom(rval->getMomentumAtEndpoint()));

      const auto [iterator, inserted] = mcparticlesMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element" << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }
    return dest;
  }

  std::unique_ptr<edm4hep::ReconstructedParticleCollection> convertReconstructedParticle(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle>& recoparticlesMap)
  {
    auto dest = std::make_unique<edm4hep::ReconstructedParticleCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::ReconstructedParticle*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setCharge(rval->getCharge());
      auto& m = rval->getCovMatrix(); // 10 parameters
      lval.setCovMatrix({m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9]});
      lval.setEnergy(rval->getEnergy());
      lval.setGoodnessOfPID(rval->getGoodnessOfPID());
      lval.setMass(rval->getMass());
      lval.setMomentum(Vector3fFrom(rval->getMomentum()));
      lval.setReferencePoint(rval->getReferencePoint());
      lval.setType(rval->getType());

      const auto [iterator, inserted] = recoparticlesMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }
    return dest;
  }

  std::unique_ptr<edm4hep::VertexCollection> convertVertex(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::Vertex*, edm4hep::MutableVertex>& vertexMap)
  {
    auto dest = std::make_unique<edm4hep::VertexCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::Vertex*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setPrimary(rval->isPrimary() ? 1 : 0); // 1 for primary and 0 for not primary
      lval.setChi2(rval->getChi2());
      lval.setProbability(rval->getProbability());
      lval.setPosition(rval->getPosition());
      auto& m = rval->getCovMatrix(); // 6 parameters
      lval.setCovMatrix({m[0], m[1], m[2], m[3], m[4], m[5]});
      // FIXME: the algorithm type in LCIO is a string, but an integer is expected
      // lval.setAlgorithmType(rval->getAlgorithmType());
      // lval.setAssociatedParticle();  //set it when convert ReconstructedParticle
      //
      for (auto v : rval->getParameters()) {
        lval.addToParameters(v);
      }

      const auto [iterator, inserted] = vertexMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }
    return dest;
  }

  std::unique_ptr<edm4hep::SimTrackerHitCollection> convertSimTrackerHit(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit>& SimTrHitMap)
  {
    auto dest = std::make_unique<edm4hep::SimTrackerHitCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::SimTrackerHit*>(LCCollection->getElementAt(i));
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

      const auto [iterator, inserted] = SimTrHitMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }
    return dest;
  }

  std::unique_ptr<edm4hep::TPCHitCollection> convertTPCHit(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::TPCHit*, edm4hep::MutableTPCHit>& TPCHitMap)
  {
    auto dest = std::make_unique<edm4hep::TPCHitCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::TPCHit*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setCellID(rval->getCellID());
      lval.setTime(rval->getTime());
      lval.setCharge(rval->getCharge());
      lval.setQuality(rval->getQuality());
      for (unsigned j = 0, M = rval->getNRawDataWords(); j < M; j++) {
        lval.addToRawDataWords(rval->getRawDataWord(j));
      }
      const auto [iterator, inserted] = TPCHitMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  std::unique_ptr<edm4hep::TrackerHitCollection> convertTrackerHit(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::TrackerHit*, edm4hep::MutableTrackerHit>& TrackerHitMap)
  {
    auto dest = std::make_unique<edm4hep::TrackerHitCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::TrackerHit*>(LCCollection->getElementAt(i));
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

      const auto [iterator, inserted] = TrackerHitMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }
    return dest;
  }

  std::unique_ptr<edm4hep::TrackerHitPlaneCollection> convertTrackerHitPlane(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane>& TrackerHitPlaneMap)
  {
    auto dest = std::make_unique<edm4hep::TrackerHitPlaneCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::TrackerHitPlane*>(LCCollection->getElementAt(i));
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

      const auto [iterator, inserted] = TrackerHitPlaneMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  std::unique_ptr<edm4hep::TrackCollection> convertTrack(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::Track*, edm4hep::MutableTrack>& TrackMap)
  {
    auto dest = std::make_unique<edm4hep::TrackCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::Track*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setType(rval->getType());
      lval.setChi2(rval->getChi2());
      lval.setNdf(rval->getNdf());
      lval.setDEdx(rval->getdEdx());
      lval.setDEdxError(rval->getdEdxError());
      lval.setRadiusOfInnermostHit(rval->getRadiusOfInnermostHit());

      auto subdetectorHitNum = rval->getSubdetectorHitNumbers();
      for (auto hitNum : subdetectorHitNum) {
        lval.addToSubDetectorHitNumbers(hitNum);
      }
      auto& trackStates = rval->getTrackStates();
      for (auto& trackState : trackStates) {
        lval.addToTrackStates(convertTrackState(trackState));
      }
      auto quantities = edm4hep::Quantity {};
      quantities.value = rval->getdEdx();
      quantities.error = rval->getdEdxError();
      lval.addToDxQuantities(quantities);
      const auto [iterator, inserted] = TrackMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  std::vector<CollNamePair> convertSimCalorimeterHit(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit>& SimCaloHitMap)
  {
    auto dest = std::make_unique<edm4hep::SimCalorimeterHitCollection>();
    auto contr = std::make_unique<edm4hep::CaloHitContributionCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::SimCalorimeterHit*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      uint64_t cellID = rval->getCellID1();
      cellID = (cellID << 32) | rval->getCellID0();
      lval.setCellID(cellID);
      lval.setEnergy(rval->getEnergy());
      lval.setPosition(rval->getPosition());

      auto NMCParticle = rval->getNMCParticles();
      for (unsigned j = 0; j < NMCParticle; j++) {
        auto edm_contr = contr->create();

        edm_contr.setPDG(rval->getPDGCont(j));
        edm_contr.setTime(rval->getTimeCont(j));
        edm_contr.setEnergy(rval->getEnergyCont(j));
        edm_contr.setStepPosition(rval->getStepPosition(j));

        lval.addToContributions(edm_contr);
      }
      const auto [iterator, inserted] = SimCaloHitMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    std::vector<CollNamePair> results;
    results.emplace_back(name, std::move(dest));
    results.emplace_back(name + "_contribution", std::move(contr));

    return results;
  }

  std::unique_ptr<edm4hep::RawCalorimeterHitCollection> convertRawCalorimeterHit(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::RawCalorimeterHit*, edm4hep::MutableRawCalorimeterHit>& rawCaloHitMap)
  {
    auto dest = std::make_unique<edm4hep::RawCalorimeterHitCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::RawCalorimeterHit*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      uint64_t cellID = rval->getCellID1();
      cellID = (cellID << 32) | rval->getCellID0();
      lval.setCellID(cellID);
      lval.setAmplitude(rval->getAmplitude());
      lval.setTimeStamp(rval->getTimeStamp());

      const auto [iterator, inserted] = rawCaloHitMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  std::unique_ptr<edm4hep::CalorimeterHitCollection> convertCalorimeterHit(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit>& caloHitMap)
  {
    auto dest = std::make_unique<edm4hep::CalorimeterHitCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::CalorimeterHit*>(LCCollection->getElementAt(i));
      auto lval = dest->create();
      uint64_t cellID = rval->getCellID1();
      cellID = (cellID << 32) | rval->getCellID0();
      lval.setCellID(cellID);
      lval.setEnergy(rval->getEnergy());
      lval.setEnergyError(rval->getEnergyError());
      lval.setPosition(rval->getPosition());
      lval.setTime(rval->getTime());
      lval.setType(rval->getType());

      const auto [iterator, inserted] = caloHitMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->second;
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  std::unique_ptr<edm4hep::ParticleIDCollection> convertParticleID(
    const std::string& name,
    EVENT::LCCollection* LCCollection)
  {
    auto dest = std::make_unique<edm4hep::ParticleIDCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::ParticleID*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setType(rval->getType());
      lval.setPDG(rval->getPDG());
      lval.setAlgorithmType(rval->getAlgorithmType());
      lval.setLikelihood(rval->getLikelihood());

      for (auto v : rval->getParameters()) {
        lval.addToParameters(v);
      }
    }

    return dest;
  }

  std::unique_ptr<edm4hep::ClusterCollection> convertCluster(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    TypeMapT<const lcio::Cluster*, edm4hep::MutableCluster>& clusterMap)
  {
    auto dest = std::make_unique<edm4hep::ClusterCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      const auto* rval = static_cast<EVENT::Cluster*>(LCCollection->getElementAt(i));
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

      const auto [iterator, inserted] = clusterMap.emplace(rval, lval);
      if (!inserted) {
        auto existing = iterator->first;
        const auto existingId = existing->id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  std::vector<CollNamePair>
  convertCollection(const std::string& name, EVENT::LCCollection* LCCollection, LcioEdmTypeMapping& typeMapping)
  {
    const auto& type = LCCollection->getTypeName();
    std::vector<CollNamePair> retColls;
    if (type == "MCParticle") {
      retColls.emplace_back(name, convertMCParticle(name, LCCollection, typeMapping.mcParticles));
    }
    else if (type == "ReconstructedParticle") {
      retColls.emplace_back(name, convertReconstructedParticle(name, LCCollection, typeMapping.recoParticles));
    }
    else if (type == "Vertex") {
      retColls.emplace_back(name, convertVertex(name, LCCollection, typeMapping.vertices));
    }
    else if (type == "Track") {
      retColls.emplace_back(name, convertTrack(name, LCCollection, typeMapping.tracks));
    }
    else if (type == "Cluster") {
      retColls.emplace_back(name, convertCluster(name, LCCollection, typeMapping.clusters));
    }
    else if (type == "SimCalorimeterHit") {
      return convertSimCalorimeterHit(name, LCCollection, typeMapping.simCaloHits);
    }
    else if (type == "RawCalorimeterHit") {
      retColls.emplace_back(name, convertRawCalorimeterHit(name, LCCollection, typeMapping.rawCaloHits));
    }
    else if (type == "CalorimeterHit") {
      retColls.emplace_back(name, convertCalorimeterHit(name, LCCollection, typeMapping.caloHits));
    }
    else if (type == "SimTrackerHit") {
      retColls.emplace_back(name, convertSimTrackerHit(name, LCCollection, typeMapping.simTrackerHits));
    }
    else if (type == "TPCHit") {
      retColls.emplace_back(name, convertTPCHit(name, LCCollection, typeMapping.tpcHits));
    }
    else if (type == "TrackerHit") {
      retColls.emplace_back(name, convertTrackerHit(name, LCCollection, typeMapping.trackerHits));
    }
    else if (type == "TrackerHitPlane") {
      retColls.emplace_back(name, convertTrackerHitPlane(name, LCCollection, typeMapping.trackerHitPlanes));
    }
    else if (type == "ParticleID") {
      retColls.emplace_back(name, convertParticleID(name, LCCollection));
    }
    return retColls;
  }

  podio::Frame convertEvent(EVENT::LCEvent* evt)
  {
    auto typeMapping = LcioEdmTypeMapping {};
    std::vector<CollNamePair> edmevent;
    std::vector<std::pair<std::string, EVENT::LCCollection*>> LCRelations;
    const auto& lcnames = evt->getCollectionNames();
    // In this loop the data gets converted.
    for (const auto& lcioname : *lcnames) {
      const auto& lcioColl = evt->getCollection(lcioname);
      const auto& lciotype = lcioColl->getTypeName();
      if (lciotype == "LCRelations") {
        LCRelations.push_back(std::make_pair(lcioname, lcioColl));
        // We handle Relations (aka Associations) once we have converted all the
        // data parts.
        continue;
      }

      if (!lcioColl->isSubset()) {
        for (auto&& [name, edmColl] : convertCollection(lcioname, lcioColl, typeMapping)) {
          if (edmColl != nullptr) {
            edmevent.emplace_back(std::move(name), std::move(edmColl));
          }
        }
      }
    }

    // Filling of the Subset Colections
    for (const auto& lcioname : *lcnames) {

      auto lcioColl = evt->getCollection(lcioname);
      if (lcioColl->isSubset()) {
        const auto& lciotype = lcioColl->getTypeName();
        auto edmColl = fillSubSet(lcioColl, typeMapping, lciotype);
        if (edmColl != nullptr) {
          edmevent.emplace_back(lcioname, std::move(edmColl));
        }
      }
    }
    // Filling all the OneToMany and OnToOne Relations and creating the AssociationCollections.
    resolveRelations(typeMapping);
    auto assoCollVec = createAssociations(typeMapping, LCRelations);

    podio::Frame event;
    // Now everything is done and we simply populate a Frame
    for (auto& [name, coll] : edmevent) {
      event.put(std::move(coll), name);
    }
    for (auto& [name, coll] : assoCollVec) {
      event.put(std::move(coll), name);
    }
    return event;
  }

  void resolveRelationsMCParticle(TypeMapT<const lcio::MCParticle*, edm4hep::MutableMCParticle>& mcparticlesMap)
  {
    int edmnum = 1;
    for (auto& [lcio, edm] : mcparticlesMap) {
      edmnum++;
      auto daughters = lcio->getDaughters();
      auto parents = lcio->getParents();

      int dnum = 1;
      for (auto d : daughters) {
        const auto it = mcparticlesMap.find(d);
        dnum++;
        if (it != mcparticlesMap.end()) {
          edm.addToDaughters(it->second);
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep MCParticle for an LCIO MCParticle, "
                       "while trying to resolve the daughters of MCParticles"
                    << std::endl;
        }
      }
      for (auto p : parents) {
        const auto it = mcparticlesMap.find(p);
        if (it != mcparticlesMap.end()) {
          edm.addToParents(it->second);
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep MCParticle for the LCIO MCParticle, "
                       "while trying to resolve the parents of MCParticles Collections"
                    << std::endl;
        }
      }
    }
  }

  void resolveRelationsSimTrackerHit(
    TypeMapT<const lcio::SimTrackerHit*, edm4hep::MutableSimTrackerHit>& SimTrHitMap,
    TypeMapT<const lcio::MCParticle*, edm4hep::MutableMCParticle>& mcparticlesMap)
  {
    for (auto& [lcio, edm] : SimTrHitMap) {
      auto mcps = lcio->getMCParticle();
      const auto it = mcparticlesMap.find(mcps);
      if (it != mcparticlesMap.end()) {
        edm.setMCParticle(it->second);
      }
      else {
        std::cerr << "Cannot find corresponding EDM4hep MCParticle for the LCIO MCParticle, "
                     "while trying to resolve the SimTrackHit Relations"
                  << std::endl;
      }
    }
  }

  void resolveRelationsRecocParticle(
    TypeMapT<const lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle>& recoparticlesMap,
    const TypeMapT<const lcio::Vertex*, edm4hep::MutableVertex>& vertexMap,
    const TypeMapT<const lcio::Cluster*, edm4hep::MutableCluster>& clusterMap,
    const TypeMapT<const lcio::Track*, edm4hep::MutableTrack>& tracksMap,
    const TypeMapT<const lcio::ParticleID*, edm4hep::MutableParticleID>& particleIDMap)
  {
    int edmnum = 1;
    for (auto& [lcio, edm] : recoparticlesMap) {
      edmnum++;

      auto vertex = lcio->getStartVertex();
      if (const auto it = vertexMap.find(vertex); it != vertexMap.end()) {
        edm.setStartVertex(it->second);
      }
      else {
        std::cerr << "Cannot find corresponding EDM4hep Vertex for a LCIO Vertex, "
                     "while trying to resolve the ReconstructedParticle Relations"
                  << std::endl;
      }

      auto clusters = lcio->getClusters();
      for (auto c : clusters) {
        const auto it = clusterMap.find(c);
        if (it != clusterMap.end()) {
          edm.addToClusters(it->second);
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep Cluster for a LCIO Cluster, "
                       "while trying to resolve the ReconstructedParticle Relations"
                    << std::endl;
        }
      }

      auto tracks = lcio->getTracks();
      for (auto t : tracks) {
        const auto it = tracksMap.find(t);
        if (it != tracksMap.end()) {
          edm.addToTracks(it->second);
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep Tracks for a LCIO Tracks, "
                       "while trying to resolve the ReconstructedParticle Relations"
                    << std::endl;
        }
      }

      auto parents = lcio->getParticles();
      for (auto p : parents) {
        const auto it = recoparticlesMap.find(p);
        if (it != recoparticlesMap.end()) {
          edm.addToParticles(it->second);
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep RecoParticle for a LCIO RecoParticle, "
                       "while trying to resolve the ReconstructedParticles parents Relations"
                    << std::endl;
        }
      }

      auto particleIDUsed = lcio->getParticleIDUsed();
      if (const auto it = particleIDMap.find(particleIDUsed); it != particleIDMap.end()) {
        edm.setParticleIDUsed(it->second);
      }
      else {
        std::cerr << "Cannot find corresponding ParticleIDUsed for a LCIO RecoParticle" << std::endl;
      }

      for (auto pid : lcio->getParticleIDs()) {
        if (const auto it = particleIDMap.find(pid); it != particleIDMap.end()) {
          edm.addToParticleIDs(it->second);
        }
        else {
          std::cerr << "Cannot find corresponding ParticleID for a LCIO Recoparticle" << std::endl;
        }
      }
    }
  }

  void resolveRelationsSimCalorimeterHit(
    TypeMapT<const lcio::SimCalorimeterHit*, edm4hep::MutableSimCalorimeterHit>& SimCaloHitMap,
    const TypeMapT<const lcio::MCParticle*, edm4hep::MutableMCParticle>& mcparticlesMap)
  {
    // TODO: Currently not doing anything here because we cannot get a mutable
    // CaloHitContribution from the SimCalorimeterHit via the podio generated
    // interface. The underlying issue is https://github.com/AIDASoft/podio/issues/347

    // for (auto& [lcio, edm] : SimCaloHitMap) {
    //   auto contributionLen = lcio->getNMCParticles();
    //   for (auto i = 0; i < contributionLen; i++) {
    //     auto mcp = lcio->getParticleCont(i);
    //     const auto it = mcparticlesMap.find(mcp);
    //     if (it != mcparticlesMap.end()) {
    //       edm.getContributions(i).setParticle(it->second); // This breaks
    //     }
    //   }
    // }
  }

  void resolveRelationsCluster(
    TypeMapT<const lcio::Cluster*, edm4hep::MutableCluster>& clustersMap,
    const TypeMapT<const lcio::CalorimeterHit*, edm4hep::MutableCalorimeterHit>& caloHitMap)
  {
    for (auto& [lcio, edm] : clustersMap) {
      auto clusters = lcio->getClusters();
      auto calohits = lcio->getCalorimeterHits();
      auto shape = lcio->getShape();
      auto subdetectorEnergies = lcio->getSubdetectorEnergies();
      for (auto c : clusters) {
        const auto it = clustersMap.find(c);
        if (it != clustersMap.end()) {
          edm.addToClusters(it->second);
        }
        else {
          std::cerr << "Couldn't find cluster to add to Relations in edm" << std::endl;
        }
      }
      for (auto cal : calohits) {
        const auto it = caloHitMap.find(cal);
        if (it != caloHitMap.end()) {
          edm.addToHits(it->second);
        }
        else {
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

  void resolveRelationsTrack(
    TypeMapT<const lcio::Track*, edm4hep::MutableTrack>& tracksMap,
    const TypeMapT<const lcio::TrackerHit*, edm4hep::MutableTrackerHit>& trackerHitMap,
    const TypeMapT<const lcio::TPCHit*, edm4hep::MutableTPCHit>& TPCHitMap,
    const TypeMapT<const lcio::TrackerHitPlane*, edm4hep::MutableTrackerHitPlane>& trackerhitplaneMap)
  {
    for (auto& [lcio, edm] : tracksMap) {
      auto tracks = lcio->getTracks();
      auto trackerHits = lcio->getTrackerHits();
      for (auto t : tracks) {
        const auto it = tracksMap.find(t);
        if (it != tracksMap.end()) {
          edm.addToTracks(it->second);
        }
        else {
          std::cerr << "Couldn't find tracks to add to Tracks Relations in edm" << std::endl;
        }
      }
      for (auto th : trackerHits) {
        const auto it = trackerHitMap.find(th);
        if (it != trackerHitMap.end()) {
          edm.addToTrackerHits(it->second);
        }
        else {
          //         std::cerr << "Couldn't find trackerHit to add to Relations for tracks in edm\n"
          //         <<" This is due to it being a TrackerHitPlane or TPCHit" << std::endl;

          /*
          // This Code looks for the trackerHit in the TPCHit Map aswell as the trackerHitPlane Map. Those relations can
          not be set for a track in edm4HEP. In all tests the missing trackerHits were located in either of these maps.
                const auto tpchit = dynamic_cast<lcio::TPCHit*>(th);
                const auto trackerhitplane = dynamic_cast<lcio::TrackerHitPlane*>(th);
                  if (tpchit != nullptr){
                    const auto it = TPCHitMap.find(tpchit);
                    if (it != TPCHitMap.end()) {
                      std::cout << "trackerHit found in TPCHit map !" << std::endl;
                      }
                      else {
                        std::cerr << "TRACKERHIT also could not be found in TPCHit Map" << std::endl;
                      }
                    }
                  else if(trackerhitplane != nullptr){
                    const auto it = trackerhitplaneMap.find(trackerhitplane);
                    if (it != trackerhitplaneMap.end()) {
                      std::cout << "trackerHit found in TrackerHitPlane map !" << std::endl;
                      }
                      else {
                        std::cerr << "TRACKERHIT also could not be found in TrackerHitPlane Map" << std::endl;
                        }
                     }
          */
        }
      }
    }
  }

  void resolveRelationsVertex(
    TypeMapT<const lcio::Vertex*, edm4hep::MutableVertex>& vertexMap,
    const TypeMapT<const lcio::ReconstructedParticle*, edm4hep::MutableReconstructedParticle>& recoparticleMap)
  {
    for (auto& [lcio, edm] : vertexMap) {
      auto recoparticle = lcio->getAssociatedParticle();

      const auto it = recoparticleMap.find(recoparticle);
      if (it != recoparticleMap.end()) {
        edm.setAssociatedParticle(it->second);
      }
      else {
        std::cerr << "Couldn't find associated Particle to add to Vertex "
                  << "Relations in edm" << std::endl;
      }
    }
  }

  void resolveRelations(LcioEdmTypeMapping& typeMapping)
  {
    resolveRelationsMCParticle(typeMapping.mcParticles);
    resolveRelationsRecoParticle(
      typeMapping.recoParticles,
      typeMapping.vertices,
      typeMapping.clusters,
      typeMapping.tracks,
      typeMapping.particleIDs);
    resolveRelationsSimTrackerHit(typeMapping.simTrackerHits, typeMapping.mcParticles);
    resolveRelationsSimCalorimeterHit(typeMapping.simCaloHits, typeMapping.mcParticles);
    resolveRelationsCluster(typeMapping.clusters, typeMapping.caloHits);
    resolveRelationsTrack(
      typeMapping.tracks, typeMapping.trackerHits, typeMapping.tpcHits, typeMapping.trackerHitPlanes);
    resolveRelationsVertex(typeMapping.vertices, typeMapping.recoParticles);
  }

  std::vector<CollNamePair> createAssociations(
    const LcioEdmTypeMapping& typeMapping,
    const std::vector<std::pair<std::string, EVENT::LCCollection*>>& LCRelation)
  {
    std::vector<CollNamePair> assoCollVec;
    for (const auto& [name, relations] : LCRelation) {
      const auto& params = relations->getParameters();

      const auto& fromType = params.getStringVal("FromType");
      const auto& toType = params.getStringVal("ToType");

      if (fromType == "MCParticle" && toType == "ReconstructedParticle") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoParticleAssociationCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.recoParticles);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "ReconstructedParticle" && toType == "MCParticle") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoParticleAssociationCollection, true>(
          relations, typeMapping.recoParticles, typeMapping.mcParticles);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "CalorimeterHit" && toType == "SimCalorimeterHit") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoCaloAssociationCollection, true>(
          relations, typeMapping.caloHits, typeMapping.simCaloHits);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "SimCalorimeterHit" && toType == "CalorimeterHit") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoCaloAssociationCollection, false>(
          relations, typeMapping.simCaloHits, typeMapping.caloHits);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "Cluster" && toType == "MCParticle") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoClusterParticleAssociationCollection, true>(
          relations, typeMapping.clusters, typeMapping.mcParticles);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "MCParticle" && toType == "Cluster") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoClusterParticleAssociationCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.clusters);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "MCParticle" && toType == "Track") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoTrackParticleAssociationCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.tracks);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "Track" && toType == "MCParticle") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoTrackParticleAssociationCollection, true>(
          relations, typeMapping.tracks, typeMapping.mcParticles);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "TrackerHit" && toType == "SimTrackerHit") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoTrackerAssociationCollection, true>(
          relations, typeMapping.trackerHits, typeMapping.simTrackerHits);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "SimTrackerHit" && toType == "TrackerHit") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoTrackerAssociationCollection, false>(
          relations, typeMapping.simTrackerHits, typeMapping.trackerHits);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "SimTrackerHit" && toType == "TrackerHitPlane") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoTrackerHitPlaneAssociationCollection, false>(
          relations, typeMapping.simTrackerHits, typeMapping.trackerHitPlanes);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "TrackerHitPlane" && toType == "SimTrackerHit") {
        auto mc_a = createAssociationCollection<edm4hep::MCRecoTrackerHitPlaneAssociationCollection, true>(
          relations, typeMapping.trackerHitPlanes, typeMapping.simTrackerHits);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "ReconstructedParticle" && toType == "Vertex") {
        auto mc_a = createAssociationCollection<edm4hep::RecoParticleVertexAssociationCollection, true>(
          relations, typeMapping.recoParticles, typeMapping.vertices);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "Vertex" && toType == "reconstructedparticle") {
        auto mc_a = createAssociationCollection<edm4hep::RecoParticleVertexAssociationCollection, false>(
          relations, typeMapping.vertices, typeMapping.recoParticles);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
    }

    return assoCollVec;
  }

  std::unique_ptr<podio::CollectionBase>
  fillSubSet(EVENT::LCCollection* LCCollection, const LcioEdmTypeMapping& typeMapping, const std::string& type)
  {
    if (type == "MCParticle") {
      return handleSubsetColl<edm4hep::MCParticleCollection>(LCCollection, typeMapping.mcParticles);
    }
    else if (type == "ReconstructedParticle") {
      return handleSubsetColl<edm4hep::ReconstructedParticleCollection>(LCCollection, typeMapping.recoParticles);
    }
    else if (type == "Vertex") {
      return handleSubsetColl<edm4hep::VertexCollection>(LCCollection, typeMapping.vertices);
    }
    else if (type == "Track") {
      return handleSubsetColl<edm4hep::TrackCollection>(LCCollection, typeMapping.tracks);
    }
    else if (type == "Cluster") {
      return handleSubsetColl<edm4hep::ClusterCollection>(LCCollection, typeMapping.clusters);
    }
    else if (type == "SimCalorimeterHit") {
      return handleSubsetColl<edm4hep::SimCalorimeterHitCollection>(LCCollection, typeMapping.simCaloHits);
    }
    else if (type == "RawCalorimeterHit") {
      return handleSubsetColl<edm4hep::RawCalorimeterHitCollection>(LCCollection, typeMapping.rawCaloHits);
    }
    else if (type == "CalorimeterHit") {
      return handleSubsetColl<edm4hep::CalorimeterHitCollection>(LCCollection, typeMapping.caloHits);
    }
    else if (type == "SimTrackerHit") {
      return handleSubsetColl<edm4hep::SimTrackerHitCollection>(LCCollection, typeMapping.simTrackerHits);
    }
    else if (type == "TPCHit") {
      return handleSubsetColl<edm4hep::TPCHitCollection>(LCCollection, typeMapping.tpcHits);
    }
    else if (type == "TrackerHit") {
      return handleSubsetColl<edm4hep::TrackerHitCollection>(LCCollection, typeMapping.trackerHits);
    }
    else if (type == "TrackerHitPlane") {
      return handleSubsetColl<edm4hep::TrackerHitPlaneCollection>(LCCollection, typeMapping.trackerHitPlanes);
    }
    else if (type == "ParticleID") {
      return handleSubsetColl<edm4hep::ParticleIDCollection>(LCCollection, typeMapping.particleIDs);
    }
    else {
      return nullptr;
    }
  }

} // namespace LCIO2EDM4hepConv
