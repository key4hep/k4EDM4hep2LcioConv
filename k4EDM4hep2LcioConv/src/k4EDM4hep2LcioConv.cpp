#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.h"

#include "edm4hep/Constants.h"
#include "edm4hep/utils/ParticleIDUtils.h"

#include "UTIL/PIDHandler.h"
#include <edm4hep/ParticleIDCollection.h>

#include <limits>
#include <algorithm>

namespace EDM4hep2LCIOConv {

  void convEventHeader(const edm4hep::EventHeaderCollection* const header_coll, lcio::LCEventImpl* const lcio_event)
  {
    convertEventHeader(header_coll, lcio_event);
  }

  void convertEventHeader(const edm4hep::EventHeaderCollection* const header_coll, lcio::LCEventImpl* const lcio_event)
  {
    if (header_coll->size() != 1) {
      return;
    }

    const auto& header = header_coll->at(0);
    lcio_event->setEventNumber(header.getEventNumber());
    lcio_event->setRunNumber(header.getRunNumber());
    lcio_event->setTimeStamp(header.getTimeStamp());
    lcio_event->setWeight(header.getWeight());
  }

  // Check if a collection is already in the event by its name
  bool collectionExist(const std::string& collection_name, const lcio::LCEventImpl* lcio_event)
  {
    const auto* coll = lcio_event->getCollectionNames();
    return std::find(coll->begin(), coll->end(), collection_name) != coll->end();
  }

  std::tuple<std::string, std::string> getPidAlgoName(const std::string& collectionName)
  {
    const auto pidPos = collectionName.find("_PID_");
    return {collectionName.substr(pidPos + 5), collectionName.substr(0, pidPos)};
  }

  void sortParticleIDs(std::vector<ParticleIDConvData>& pidCollections)
  {
    std::sort(pidCollections.begin(), pidCollections.end(), [](const auto& pid1, const auto& pid2) {
      static const auto defaultPidMeta = edm4hep::utils::ParticleIDMeta {"", std::numeric_limits<int>::max(), {}};
      return pid1.metadata.value_or(defaultPidMeta).algoType < pid2.metadata.value_or(defaultPidMeta).algoType;
    });
  }

  std::optional<int32_t> attachParticleIDMetaData(
    IMPL::LCEventImpl* lcEvent,
    const podio::Frame& edmEvent,
    const ParticleIDConvData& pidCollMetaInfo)
  {
    const auto& [name, coll, pidMetaInfo] = pidCollMetaInfo;
    const auto recoName = edmEvent.getName((*coll)[0].getParticle().id().collectionID);
    // If we can't get the reconstructed particle collection name there is not
    // much we can do
    if (!recoName.has_value()) {
      return std::nullopt;
    }
    // If we can't get meta data information there is not much we can do either
    if (!pidMetaInfo.has_value()) {
      return std::nullopt;
    }
    if (pidMetaInfo.has_value() && !recoName.has_value()) {
      return pidMetaInfo->algoType;
    }

    UTIL::PIDHandler pidHandler(lcEvent->getCollection(recoName.value()));
    return pidHandler.addAlgorithm(pidMetaInfo->algoName, pidMetaInfo->paramNames);
  }

  std::unique_ptr<lcio::LCEventImpl> convertEvent(const podio::Frame& edmEvent, const podio::Frame& metadata)
  {
    auto lcioEvent = std::make_unique<lcio::LCEventImpl>();
    auto objectMappings = CollectionsPairVectors {};

    // We have to convert these after all other (specifically
    // ReconstructedParticle) collections have been converted. Otherwise we will
    // not be able to set all the metadata for the PIDHandler (LCIO) to work
    // properly. Here we store the name, the collection as well as potentially
    // available meta information that we obtain when we first come across a
    // collection
    std::vector<ParticleIDConvData> pidCollections {};

    const auto& collections = edmEvent.getAvailableCollections();
    for (const auto& name : collections) {
      const auto edmCollection = edmEvent.get(name);

      const auto& cellIDStr =
        metadata.getParameter<std::string>(podio::collMetadataParamName(name, edm4hep::CellIDEncoding));

      if (auto coll = dynamic_cast<const edm4hep::TrackCollection*>(edmCollection)) {
        auto lcColl = convertTracks(coll, objectMappings.tracks);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::TrackerHit3DCollection*>(edmCollection)) {
        auto lcColl = convertTrackerHits(coll, cellIDStr, objectMappings.trackerHits);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::TrackerHitPlaneCollection*>(edmCollection)) {
        auto lcColl = convertTrackerHitPlanes(coll, cellIDStr, objectMappings.trackerHitPlanes);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::SimTrackerHitCollection*>(edmCollection)) {
        auto lcColl = convertSimTrackerHits(coll, cellIDStr, objectMappings.simTrackerHits);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::CalorimeterHitCollection*>(edmCollection)) {
        auto lcColl = convertCalorimeterHits(coll, cellIDStr, objectMappings.caloHits);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::RawCalorimeterHitCollection*>(edmCollection)) {
        auto lcColl = convertRawCalorimeterHits(coll, objectMappings.rawCaloHits);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::SimCalorimeterHitCollection*>(edmCollection)) {
        auto lcColl = convertSimCalorimeterHits(coll, cellIDStr, objectMappings.simCaloHits);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::RawTimeSeriesCollection*>(edmCollection)) {
        auto lcColl = convertTPCHits(coll, objectMappings.tpcHits);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::ClusterCollection*>(edmCollection)) {
        auto lcColl = convertClusters(coll, objectMappings.clusters);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::VertexCollection*>(edmCollection)) {
        auto lcColl = convertVertices(coll, objectMappings.vertices);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::MCParticleCollection*>(edmCollection)) {
        auto lcColl = convertMCParticles(coll, objectMappings.mcParticles);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::ReconstructedParticleCollection*>(edmCollection)) {
        auto lcColl = convertReconstructedParticles(coll, objectMappings.recoParticles);
        lcioEvent->addCollection(lcColl.release(), name);
      }
      else if (auto coll = dynamic_cast<const edm4hep::EventHeaderCollection*>(edmCollection)) {
        convertEventHeader(coll, lcioEvent.get());
      }
      else if (auto coll = dynamic_cast<const edm4hep::ParticleIDCollection*>(edmCollection)) {
        pidCollections.emplace_back(name, coll, edm4hep::utils::PIDHandler::getAlgoInfo(metadata, name));
      }
      else if (dynamic_cast<const edm4hep::CaloHitContributionCollection*>(edmCollection)) {
        // "converted" during relation resolving later
        continue;
      }
      else {
        std::cerr << "Error trying to convert requested " << edmCollection->getValueTypeName() << " with name " << name
                  << "\n"
                  << "List of supported types: "
                  << "Track, TrackerHit, TrackerHitPlane, SimTrackerHit, "
                  << "Cluster, CalorimeterHit, RawCalorimeterHit, "
                  << "SimCalorimeterHit, Vertex, ReconstructedParticle, "
                  << "MCParticle." << std::endl;
      }
    }

    sortParticleIDs(pidCollections);
    for (const auto& pidCollMeta : pidCollections) {
      // Use -1 as a somewhat easy to identify value of missing reco collections
      // or pid metadata
      const auto algoId = attachParticleIDMetaData(lcioEvent.get(), edmEvent, pidCollMeta).value_or(-1);
      convertParticleIDs(pidCollMeta.coll, objectMappings.particleIDs, algoId);
    }

    resolveRelations(objectMappings);

    return lcioEvent;
  }

} // namespace EDM4hep2LCIOConv
