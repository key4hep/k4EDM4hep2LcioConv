#include "k4EDM4hep2LcioConv/MappingUtils.h"

namespace LCIO2EDM4hepConv {
  template<typename LCIOType>
  void convertObjectParameters(LCIOType* lcioobj, podio::Frame& event)
  {
    const auto& params = lcioobj->getParameters();
    // handle srting params
    EVENT::StringVec keys;
    const auto stringKeys = params.getStringKeys(keys);
    for (int i = 0; i < stringKeys.size(); i++) {
      EVENT::StringVec sValues;
      const auto stringVals = params.getStringVals(stringKeys[i], sValues);
      event.putParameter(stringKeys[i], stringVals);
    }
    // handle float params
    EVENT::StringVec fkeys;
    const auto floatKeys = params.getFloatKeys(fkeys);
    for (int i = 0; i < floatKeys.size(); i++) {
      EVENT::FloatVec fValues;
      const auto floatVals = params.getFloatVals(floatKeys[i], fValues);
      event.putParameter(floatKeys[i], floatVals);
    }
    // handle int params
    EVENT::StringVec ikeys;
    const auto intKeys = params.getIntKeys(ikeys);
    for (int i = 0; i < intKeys.size(); i++) {
      EVENT::IntVec iValues;
      const auto intVals = params.getIntVals(intKeys[i], iValues);
      event.putParameter(intKeys[i], intVals);
    }
    // handle double params
    EVENT::StringVec dkeys;
    const auto dKeys = params.getDoubleKeys(dkeys);
    for (int i = 0; i < dKeys.size(); i++) {
      EVENT::DoubleVec dValues;
      const auto dVals = params.getDoubleVals(dKeys[i], dValues);
      event.putParameter(dKeys[i], dVals);
    }
  }

  template<typename LCVecType>
  std::vector<CollNamePair> convertLCVec(const std::string& name, EVENT::LCCollection* LCCollection)
  {
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

  template<typename CollT, typename ObjectMapT, typename LcioT, typename Edm4hepT>
  auto handleSubsetColl(EVENT::LCCollection* lcioColl, const ObjectMapT& elemMap)
  {
    auto edm4hepColl = std::make_unique<CollT>();
    edm4hepColl->setSubsetCollection();

    UTIL::LCIterator<LcioT> lcioIter(lcioColl);
    while (const auto lcioElem = lcioIter.next()) {
      if (auto edm4hepElem = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioElem, elemMap)) {
        edm4hepColl->push_back(edm4hepElem.value());
      }
      else {
        std::cerr << "Cannot find corresponding EDM4hep object for an LCIO object in a subset collection of type "
                  << Edm4hepT::collection_type::valueTypeName << std::endl;
      }
    }

    return edm4hepColl;
  }

  template<typename MCParticleMapT>
  std::unique_ptr<edm4hep::MCParticleCollection>
  convertMCParticles(const std::string& name, EVENT::LCCollection* LCCollection, MCParticleMapT& mcparticlesMap)
  {
    auto dest = std::make_unique<edm4hep::MCParticleCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      auto* rval = static_cast<EVENT::MCParticle*>(LCCollection->getElementAt(i));
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

  template<typename RecoMapT, typename PIDMapT>
  std::vector<CollNamePair> convertReconstructedParticles(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    RecoMapT& recoparticlesMap,
    PIDMapT& particleIDMap)
  {
    auto dest = std::make_unique<edm4hep::ReconstructedParticleCollection>();
    auto particleIDs = std::make_unique<edm4hep::ParticleIDCollection>();
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
      lval.setType(rval->getType());

      const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, recoparticlesMap);
      if (!inserted) {
        auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }

      // Need to convert the particle IDs here, since they are part of the reco
      // particle collection in LCIO
      for (const auto lcioPid : rval->getParticleIDs()) {
        auto pid = convertParticleID(lcioPid);
        const auto [pidIt, pidInserted] = k4EDM4hep2LcioConv::detail::mapInsert(
          lcioPid, pid, particleIDMap, k4EDM4hep2LcioConv::detail::InsertMode::Checked);
        if (pidInserted) {
          lval.addToParticleIDs(pid);
          particleIDs->push_back(pid);
        }
        else {
          lval.addToParticleIDs(k4EDM4hep2LcioConv::detail::getMapped(pidIt));
        }
      }

      const auto lcioPidUsed = rval->getParticleIDUsed();
      if (lcioPidUsed == nullptr) {
        continue;
      }
      if (const auto edm4hepPid = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioPidUsed, particleIDMap)) {
        lval.setParticleIDUsed(edm4hepPid.value());
      }
      else {
        auto pid = convertParticleID(lcioPidUsed);
        particleIDs->push_back(pid);
        k4EDM4hep2LcioConv::detail::mapInsert(lcioPidUsed, pid, particleIDMap);
        lval.setParticleIDUsed(pid);
      }
    }

    std::vector<CollNamePair> results;
    results.reserve(2);
    results.emplace_back(name, std::move(dest));
    results.emplace_back(name + "_particleIDs", std::move(particleIDs));
    return results;
  }

  template<typename VertexMapT>
  std::unique_ptr<edm4hep::VertexCollection>
  convertVertices(const std::string& name, EVENT::LCCollection* LCCollection, VertexMapT& vertexMap)
  {
    auto dest = std::make_unique<edm4hep::VertexCollection>();
    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      auto* rval = static_cast<EVENT::Vertex*>(LCCollection->getElementAt(i));
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

      const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, vertexMap);
      if (!inserted) {
        auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }
    return dest;
  }

  template<typename SimTrHitMapT>
  std::unique_ptr<edm4hep::SimTrackerHitCollection>
  convertSimTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, SimTrHitMapT& SimTrHitMap)
  {
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

  template<typename HitMapT>
  std::unique_ptr<edm4hep::RawTimeSeriesCollection>
  convertTPCHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TPCHitMap)
  {
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

  template<typename HitMapT>
  std::unique_ptr<edm4hep::TrackerHit3DCollection>
  convertTrackerHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitMap)
  {
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

  template<typename HitMapT>
  std::unique_ptr<edm4hep::TrackerHitPlaneCollection>
  convertTrackerHitPlanes(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& TrackerHitPlaneMap)
  {
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

  template<typename TrackMapT>
  std::unique_ptr<edm4hep::TrackCollection>
  convertTracks(const std::string& name, EVENT::LCCollection* LCCollection, TrackMapT& TrackMap)
  {
    auto dest = std::make_unique<edm4hep::TrackCollection>();

    for (unsigned i = 0, N = LCCollection->getNumberOfElements(); i < N; ++i) {
      auto* rval = static_cast<EVENT::Track*>(LCCollection->getElementAt(i));
      auto lval = dest->create();

      lval.setType(rval->getType());
      lval.setChi2(rval->getChi2());
      lval.setNdf(rval->getNdf());
      lval.setDEdx(rval->getdEdx());
      lval.setDEdxError(rval->getdEdxError());
      lval.setRadiusOfInnermostHit(rval->getRadiusOfInnermostHit());

      auto subdetectorHitNum = rval->getSubdetectorHitNumbers();
      for (auto hitNum : subdetectorHitNum) {
        lval.addToSubdetectorHitNumbers(hitNum);
      }
      auto& trackStates = rval->getTrackStates();
      for (auto& trackState : trackStates) {
        lval.addToTrackStates(convertTrackState(trackState));
      }
      auto quantities = edm4hep::Quantity {};
      quantities.value = rval->getdEdx();
      quantities.error = rval->getdEdxError();
      lval.addToDxQuantities(quantities);
      const auto [iterator, inserted] = k4EDM4hep2LcioConv::detail::mapInsert(rval, lval, TrackMap);
      if (!inserted) {
        auto existing = k4EDM4hep2LcioConv::detail::getMapped(iterator);
        const auto existingId = existing.id();
        std::cerr << "EDM4hep element  " << existingId << " did not get inserted. It belongs to the " << name
                  << " collection" << std::endl;
      }
    }

    return dest;
  }

  template<typename HitMapT>
  std::unique_ptr<edm4hep::SimCalorimeterHitCollection>
  convertSimCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& SimCaloHitMap)
  {
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

  template<typename HitMapT>
  std::unique_ptr<edm4hep::RawCalorimeterHitCollection>
  convertRawCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& rawCaloHitMap)
  {
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

  template<typename HitMapT>
  std::unique_ptr<edm4hep::CalorimeterHitCollection>
  convertCalorimeterHits(const std::string& name, EVENT::LCCollection* LCCollection, HitMapT& caloHitMap)
  {
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

  template<typename ClusterMapT, typename PIDMapT>
  std::vector<CollNamePair> convertClusters(
    const std::string& name,
    EVENT::LCCollection* LCCollection,
    ClusterMapT& clusterMap,
    PIDMapT& particleIDMap)
  {
    auto particleIDs = std::make_unique<edm4hep::ParticleIDCollection>();
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

      // Need to convert the particle IDs here, since they are part of the cluster
      // collection in LCIO
      for (const auto lcioPid : rval->getParticleIDs()) {
        auto pid = convertParticleID(lcioPid);
        const auto [pidIt, pidInserted] = k4EDM4hep2LcioConv::detail::mapInsert(lcioPid, pid, particleIDMap);
        if (pidInserted) {
          lval.addToParticleIDs(pid);
          particleIDs->push_back(pid);
        }
        else {
          lval.addToParticleIDs(k4EDM4hep2LcioConv::detail::getMapped(pidIt));
        }
      }
    }
    std::vector<CollNamePair> results;
    results.reserve(2);
    results.emplace_back(name, std::move(dest));
    results.emplace_back(name + "_particleIDs", std::move(particleIDs));
    return results;
  }

  template<typename ObjectMappingT>
  std::vector<CollNamePair>
  convertCollection(const std::string& name, EVENT::LCCollection* LCCollection, ObjectMappingT& typeMapping)
  {
    const auto& type = LCCollection->getTypeName();
    std::vector<CollNamePair> retColls;
    if (type == "MCParticle") {
      retColls.emplace_back(name, convertMCParticles(name, LCCollection, typeMapping.mcParticles));
    }
    else if (type == "ReconstructedParticle") {
      return convertReconstructedParticles(name, LCCollection, typeMapping.recoParticles, typeMapping.particleIDs);
    }
    else if (type == "Vertex") {
      retColls.emplace_back(name, convertVertices(name, LCCollection, typeMapping.vertices));
    }
    else if (type == "Track") {
      retColls.emplace_back(name, convertTracks(name, LCCollection, typeMapping.tracks));
    }
    else if (type == "Cluster") {
      return convertClusters(name, LCCollection, typeMapping.clusters, typeMapping.particleIDs);
    }
    else if (type == "SimCalorimeterHit") {
      retColls.emplace_back(name, convertSimCalorimeterHits(name, LCCollection, typeMapping.simCaloHits));
    }
    else if (type == "RawCalorimeterHit") {
      retColls.emplace_back(name, convertRawCalorimeterHits(name, LCCollection, typeMapping.rawCaloHits));
    }
    else if (type == "CalorimeterHit") {
      retColls.emplace_back(name, convertCalorimeterHits(name, LCCollection, typeMapping.caloHits));
    }
    else if (type == "SimTrackerHit") {
      retColls.emplace_back(name, convertSimTrackerHits(name, LCCollection, typeMapping.simTrackerHits));
    }
    else if (type == "TPCHit") {
      retColls.emplace_back(name, convertTPCHits(name, LCCollection, typeMapping.tpcHits));
    }
    else if (type == "TrackerHit") {
      retColls.emplace_back(name, convertTrackerHits(name, LCCollection, typeMapping.trackerHits));
    }
    else if (type == "TrackerHitPlane") {
      retColls.emplace_back(name, convertTrackerHitPlanes(name, LCCollection, typeMapping.trackerHitPlanes));
    }
    else if (type == "LCIntVec") {
      return convertLCVec<EVENT::LCIntVec>(name, LCCollection);
    }
    else if (type == "LCFloatVec") {
      return convertLCVec<EVENT::LCFloatVec>(name, LCCollection);
    }
    else if (type != "LCRelation") {
      std::cerr << type << " is a collection type for which no known conversion exists." << std::endl;
    }
    return retColls;
  }

  template<typename HitMapT, typename MCParticleMapT>
  std::unique_ptr<edm4hep::CaloHitContributionCollection> createCaloHitContributions(
    HitMapT& SimCaloHitMap,
    const MCParticleMapT& mcparticlesMap)
  {
    auto contrCollection = std::make_unique<edm4hep::CaloHitContributionCollection>();
    for (auto& [lcioHit, edmHit] : SimCaloHitMap) {
      auto NMCParticle = lcioHit->getNMCParticles();
      for (unsigned j = 0; j < NMCParticle; j++) {
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
          }
          else {
            std::cerr << "Cannot find corresponding EDM4hep MCParticle for a LCIO MCParticle, "
                      << "while trying to build CaloHitContributions " << std::endl;
          }
        }
      }
    }
    return contrCollection;
  }

  template<typename MCParticleMapT, typename MCParticleLookupMapT>
  void resolveRelationsMCParticles(MCParticleMapT& mcparticlesMap, const MCParticleLookupMapT& lookupMap)
  {
    int edmnum = 1;
    for (auto& [lcio, edm] : mcparticlesMap) {
      edmnum++;
      auto daughters = lcio->getDaughters();
      auto parents = lcio->getParents();

      for (auto d : daughters) {
        if (d == nullptr) {
          continue;
        }
        if (const auto edmD = k4EDM4hep2LcioConv::detail::mapLookupTo(d, lookupMap)) {
          edm.addToDaughters(edmD.value());
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep MCParticle for an LCIO MCParticle, "
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
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep MCParticle for the LCIO MCParticle, "
                       "while trying to resolve the parents of MCParticles Collections"
                    << std::endl;
        }
      }
    }
  }

  template<typename HitMapT, typename MCParticleMapT>
  void resolveRelationsSimTrackerHits(HitMapT& SimTrHitMap, const MCParticleMapT& mcparticlesMap)
  {
    for (auto& [lcio, edm] : SimTrHitMap) {
      auto mcps = lcio->getMCParticle();
      if (mcps == nullptr) {
        continue;
      }
      if (const auto edmP = k4EDM4hep2LcioConv::detail::mapLookupTo(mcps, mcparticlesMap)) {
        edm.setMCParticle(edmP.value());
      }
      else {
        std::cerr << "Cannot find corresponding EDM4hep MCParticle for the LCIO MCParticle, "
                     "while trying to resolve the SimTrackHit Relations"
                  << std::endl;
      }
    }
  }

  template<
    typename RecoParticleMapT,
    typename RecoParticleLookupMapT,
    typename VertexMapT,
    typename ClusterMapT,
    typename TrackMapT>
  void resolveRelationsRecoParticles(
    RecoParticleMapT& recoparticlesMap,
    const RecoParticleLookupMapT& recoLookupMap,
    const VertexMapT& vertexMap,
    const ClusterMapT& clusterMap,
    const TrackMapT& tracksMap)
  {
    int edmnum = 1;
    for (auto& [lcio, edm] : recoparticlesMap) {
      edmnum++;

      const auto& vertex = lcio->getStartVertex();
      if (vertex != nullptr) {
        if (const auto edmV = k4EDM4hep2LcioConv::detail::mapLookupTo(vertex, vertexMap)) {
          edm.setStartVertex(edmV.value());
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep Vertex for a LCIO Vertex, "
                       "while trying to resolve the ReconstructedParticle Relations "
                    << std::endl;
        }
      }

      auto clusters = lcio->getClusters();
      for (auto c : clusters) {
        if (c == nullptr) {
          continue;
        }
        if (const auto edmC = k4EDM4hep2LcioConv::detail::mapLookupTo(c, clusterMap)) {
          edm.addToClusters(edmC.value());
        }
        else {
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
        }
        else {
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
        }
        else {
          std::cerr << "Cannot find corresponding EDM4hep RecoParticle for a LCIO RecoParticle, "
                       "while trying to resolve the ReconstructedParticles parents Relations"
                    << std::endl;
        }
      }
    }
  }

  template<typename ClusterMapT, typename CaloHitMapT>
  void resolveRelationsClusters(ClusterMapT& clustersMap, const CaloHitMapT& caloHitMap)
  {
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
        }
        else {
          std::cerr << "Couldn't find cluster to add to Relations in edm" << std::endl;
        }
      }
      for (auto cal : calohits) {
        if (cal == nullptr) {
          continue;
        }
        if (const auto edmCaloHit = k4EDM4hep2LcioConv::detail::mapLookupTo(cal, caloHitMap)) {
          edm.addToHits(edmCaloHit.value());
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

  template<typename TrackMapT, typename TrackHitMapT, typename THPlaneHitMapT, typename TPCHitMapT>
  void resolveRelationsTracks(
    TrackMapT& tracksMap,
    const TrackHitMapT& trackerHitMap,
    const THPlaneHitMapT& trackerHitPlaneMap,
    const TPCHitMapT&)
  {
    for (auto& [lcio, edm] : tracksMap) {
      auto tracks = lcio->getTracks();
      for (auto t : tracks) {
        if (t == nullptr) {
          continue;
        }
        if (const auto track = k4EDM4hep2LcioConv::detail::mapLookupTo(t, tracksMap)) {
          edm.addToTracks(track.value());
        }
        else {
          // std::cerr << "Couldn't find tracks to add to Tracks Relations in edm" << std::endl;
        }
      }
      const auto trackerHits = lcio->getTrackerHits();
      for (const auto th : trackerHits) {
        bool found = false;
        if (th == nullptr) {
          continue;
        }
        if (const auto typedTH = dynamic_cast<EVENT::TrackerHitPlane*>(th)) {
          if (const auto trHit = k4EDM4hep2LcioConv::detail::mapLookupTo(typedTH, trackerHitPlaneMap)) {
            edm.addToTrackerHits(trHit.value());
            found = true;
          }
        }
        else if (auto typedTH = dynamic_cast<EVENT::TrackerHit*>(th)) {
          if (const auto trHit = k4EDM4hep2LcioConv::detail::mapLookupTo(typedTH, trackerHitMap)) {
            edm.addToTrackerHits(trHit.value());
            found = true;
          }
        }

        if (!found) {
          std::cerr << "Couldn't find a edm4hep TrackerHit for an LCIO TrackerHit when resolving "
                    << "relations for a Track" << std::endl;
        }
      }
    }
  }

  template<typename VertexMapT, typename RecoParticleMapT>
  void resolveRelationsVertices(VertexMapT& vertexMap, const RecoParticleMapT& recoparticleMap)
  {
    for (auto& [lcio, edm] : vertexMap) {
      auto recoparticle = lcio->getAssociatedParticle();
      if (recoparticle == nullptr) {
        continue;
      }
      if (const auto recoP = k4EDM4hep2LcioConv::detail::mapLookupTo(recoparticle, recoparticleMap)) {
        edm.setAssociatedParticle(recoP.value());
      }
      else {
        std::cerr << "Couldn't find associated Particle to add to Vertex "
                  << "Relations in edm" << std::endl;
      }
    }
  }

  template<typename ObjectMappingT>
  void resolveRelations(ObjectMappingT& typeMapping)
  {
    resolveRelations(typeMapping, typeMapping);
  }

  template<typename ObjectMappingT, typename ObjectMappingU>
  void resolveRelations(ObjectMappingT& updateMaps, const ObjectMappingU& lookupMaps)
  {
    resolveRelationsMCParticles(updateMaps.mcParticles, lookupMaps.mcParticles);
    resolveRelationsRecoParticles(
      updateMaps.recoParticles, lookupMaps.recoParticles, lookupMaps.vertices, lookupMaps.clusters, lookupMaps.tracks);
    resolveRelationsSimTrackerHits(updateMaps.simTrackerHits, lookupMaps.mcParticles);
    resolveRelationsClusters(updateMaps.clusters, lookupMaps.caloHits);
    resolveRelationsTracks(updateMaps.tracks, lookupMaps.trackerHits, lookupMaps.trackerHitPlanes, lookupMaps.tpcHits);
    resolveRelationsVertices(updateMaps.vertices, lookupMaps.recoParticles);
  }

  template<
    typename CollT,
    bool Reverse,
    typename FromMapT,
    typename ToMapT,
    typename FromLCIOT,
    typename ToLCIOT,
    typename FromEDM4hepT,
    typename ToEDM4hepT>
  std::unique_ptr<CollT>
  createAssociationCollection(EVENT::LCCollection* relations, const FromMapT& fromMap, const ToMapT& toMap)
  {
    auto assocColl = std::make_unique<CollT>();
    auto relIter = UTIL::LCIterator<EVENT::LCRelation>(relations);

    while (const auto rel = relIter.next()) {
      auto assoc = assocColl->create();
      assoc.setWeight(rel->getWeight());
      const auto lcioTo = static_cast<ToLCIOT*>(rel->getTo());
      const auto lcioFrom = static_cast<FromLCIOT*>(rel->getFrom());
      const auto edm4hepTo = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioTo, toMap);
      const auto edm4hepFrom = k4EDM4hep2LcioConv::detail::mapLookupTo(lcioFrom, fromMap);
      if (edm4hepTo.has_value() && edm4hepFrom.has_value()) {
        if constexpr (Reverse) {
          if constexpr (std::is_same_v<k4EDM4hep2LcioConv::detail::mutable_t<ToEDM4hepT>, edm4hep::MutableVertex>) {
            assoc.setVertex(*edm4hepTo);
          }
          else {
            assoc.setSim(*edm4hepTo);
          }
          assoc.setRec(*edm4hepFrom);
        }
        else {
          if constexpr (std::is_same_v<k4EDM4hep2LcioConv::detail::mutable_t<FromEDM4hepT>, edm4hep::MutableVertex>) {
            assoc.setVertex(*edm4hepFrom);
          }
          else {
            assoc.setSim(*edm4hepFrom);
          }
          assoc.setRec(*edm4hepTo);
        }
      }
    }

    return assocColl;
  }

  template<typename ObjectMappingT>
  std::vector<CollNamePair> createAssociations(
    const ObjectMappingT& typeMapping,
    const std::vector<std::pair<std::string, EVENT::LCCollection*>>& LCRelation)
  {
    std::vector<CollNamePair> assoCollVec;
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
      else if (fromType == "Vertex" && toType == "ReconstructedParticle") {
        auto mc_a = createAssociationCollection<edm4hep::RecoParticleVertexAssociationCollection, false>(
          relations, typeMapping.vertices, typeMapping.recoParticles);
        assoCollVec.emplace_back(name, std::move(mc_a));
      }
      else if (fromType == "CalorimeterHit" && toType == "MCParticle") {
        auto assoc = createAssociationCollection<edm4hep::MCRecoCaloParticleAssociationCollection, true>(
          relations, typeMapping.caloHits, typeMapping.mcParticles);
        assoCollVec.emplace_back(name, std::move(assoc));
      }
      else if (fromType == "MCParticle" && toType == "CalorimeterHit") {
        auto assoc = createAssociationCollection<edm4hep::MCRecoCaloParticleAssociationCollection, false>(
          relations, typeMapping.mcParticles, typeMapping.caloHits);
        assoCollVec.emplace_back(name, std::move(assoc));
      }
      else {
        std::cout << "Relation from: " << fromType << " to: " << toType << " (" << name
                  << ") is not beeing handled during creation of associations" << std::endl;
      }
    }

    return assoCollVec;
  }

  template<typename ObjectMappingT>
  std::unique_ptr<podio::CollectionBase>
  fillSubset(EVENT::LCCollection* LCCollection, const ObjectMappingT& typeMapping, const std::string& type)
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
      return handleSubsetColl<edm4hep::RawTimeSeriesCollection>(LCCollection, typeMapping.tpcHits);
    }
    else if (type == "TrackerHit") {
      return handleSubsetColl<edm4hep::TrackerHit3DCollection>(LCCollection, typeMapping.trackerHits);
    }
    else if (type == "TrackerHitPlane") {
      return handleSubsetColl<edm4hep::TrackerHitPlaneCollection>(LCCollection, typeMapping.trackerHitPlanes);
    }
    else {
      return nullptr;
    }
  }

} // namespace LCIO2EDM4hepConv
