#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.h"
#include "edm4hep/Constants.h"
#include "UTIL/PIDHandler.h"

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

  std::unique_ptr<lcio::LCEventImpl> convertEvent(const podio::Frame& edmEvent, const podio::Frame& metadata)
  {
    auto lcioEvent = std::make_unique<lcio::LCEventImpl>();
    auto objectMappings = CollectionsPairVectors {};

    // We simply collect all pairs of PID algorithm names and the RecoParticle
    // collection to which they belong (assuming the naming convention
    // introduced in the LCIO to EDM4hep conversion)
    std::vector<std::tuple<std::string, std::string>> pidAlgoNames;

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
        convertParticleIDs(coll, objectMappings.particleIDs, pidAlgoNames.size());
        pidAlgoNames.emplace_back(getPidAlgoName(name));
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

    for (const auto& [algoName, recoCollName] : pidAlgoNames) {
      std::cout << "Adding PID algo " << algoName << " to reco collection " << recoCollName << '\n';
      UTIL::PIDHandler pidHandler(lcioEvent->getCollection(recoCollName));
      std::cout << "assigned algo id " << pidHandler.addAlgorithm(algoName, {}) << '\n';
    }

    resolveRelations(objectMappings);

    return lcioEvent;
  }

} // namespace EDM4hep2LCIOConv
