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
      covMatrix[14],
      0,
      0,
      0,
      0,
      0,
      0};

    return edmtrackState;
  }

  edm4hep::MutableParticleID convertParticleID(const EVENT::ParticleID* pid)
  {
    auto result = edm4hep::MutableParticleID {};
    result.setType(pid->getType());
    result.setPDG(pid->getPDG());
    result.setAlgorithmType(pid->getAlgorithmType());
    result.setLikelihood(pid->getLikelihood());

    for (auto v : pid->getParameters()) {
      result.addToParameters(v);
    }

    return result;
  }

  std::unique_ptr<edm4hep::EventHeaderCollection> createEventHeader(const EVENT::LCEvent* evt)
  {
    auto headerColl = std::make_unique<edm4hep::EventHeaderCollection>();
    auto header = headerColl->create();

    header.setEventNumber(evt->getEventNumber());
    header.setRunNumber(evt->getRunNumber());
    header.setTimeStamp(evt->getTimeStamp());
    header.setWeight(evt->getWeight());
    return headerColl;
  }

  podio::Frame convertEvent(EVENT::LCEvent* evt, const std::vector<std::string>& collsToConvert)
  {
    auto typeMapping = LcioEdmTypeMapping {};
    std::vector<CollNamePair> edmevent;
    std::vector<std::pair<std::string, EVENT::LCCollection*>> LCRelations;

    const auto& lcioNames = [&collsToConvert, &evt]() {
      if (collsToConvert.empty()) {
        return *evt->getCollectionNames();
      }
      return collsToConvert;
    }();

    // In this loop the data gets converted.
    for (const auto& lcioname : lcioNames) {
      const auto& lcioColl = evt->getCollection(lcioname);
      const auto& lciotype = lcioColl->getTypeName();
      if (lciotype == "LCRelation") {
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
    for (const auto& lcioname : lcioNames) {

      auto lcioColl = evt->getCollection(lcioname);
      if (lcioColl->isSubset()) {
        const auto& lciotype = lcioColl->getTypeName();
        auto edmColl = fillSubset(lcioColl, typeMapping, lciotype);
        if (edmColl != nullptr) {
          edmevent.emplace_back(lcioname, std::move(edmColl));
        }
      }
    }
    // Filling all the OneToMany and OneToOne Relations and creating the AssociationCollections.
    resolveRelations(typeMapping);
    auto assoCollVec = createAssociations(typeMapping, LCRelations);
    auto headerColl = createEventHeader(evt);

    // Now everything is done and we simply populate a Frame
    podio::Frame event;
    // convert put the event parameters into the frame
    convertObjectParameters<EVENT::LCEvent>(evt, event);

    // only create CaloHitContributions if necessary (i.e. if we have converted
    // SimCalorimeterHits)
    if (not typeMapping.simCaloHits.empty()) {
      auto calocontr = createCaloHitContributions(typeMapping.simCaloHits, typeMapping.mcParticles);
      event.put(std::move(calocontr), "AllCaloHitContributionsCombined");
    }
    event.put(std::move(headerColl), "EventHeader");
    for (auto& [name, coll] : edmevent) {
      event.put(std::move(coll), name);
    }
    for (auto& [name, coll] : assoCollVec) {
      event.put(std::move(coll), name);
    }
    return event;
  }

  podio::Frame convertRunHeader(EVENT::LCRunHeader* rheader)
  {
    podio::Frame runHeaderFrame;
    runHeaderFrame.putParameter("runNumber", rheader->getRunNumber());
    runHeaderFrame.putParameter("detectorName", rheader->getDetectorName());
    runHeaderFrame.putParameter("description", rheader->getDescription());
    auto subdetectors = rheader->getActiveSubdetectors();
    runHeaderFrame.putParameter("activeSubdetectors", *subdetectors);

    // convert everything set as a parameter
    convertObjectParameters<EVENT::LCRunHeader>(rheader, runHeaderFrame);

    return runHeaderFrame;
  }

} // namespace LCIO2EDM4hepConv
