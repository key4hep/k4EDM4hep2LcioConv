#include "CompareEDM4hepLCIO.h"
#include "EDM4hep2LCIOUtilities.h"
#include "ObjectMapping.h"

#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.h"

#include "Exceptions.h"

int main() {
  const auto& [edmEvent, metadata] = createExampleEvent();
  const auto lcioEvent = EDM4hep2LCIOConv::convertEvent(edmEvent, metadata);

  if (!compareEventHeader(lcioEvent.get(), &edmEvent)) {
    return 1;
  }

  for (const auto& name : edmEvent.getAvailableCollections()) {
    const auto edmColl = edmEvent.get(name);
    const auto typeName = edmColl->getValueTypeName();
    if (typeName == "edm4hep::CaloHitContribution" || typeName == "edm4hep::ParticleID" ||
        typeName == "edm4hep::EventHeader" || typeName == "edm4hep::RecDqdx" ||
        typeName == "edm4hep::VertexRecoParticleLink") {
      continue;
    }
    try {
      const auto* lcColl = lcioEvent->getCollection(name);
      if ((unsigned)lcColl->getNumberOfElements() != edmColl->size()) {
        std::cerr << "Collection " << name << " has different sizes. EDM4hep: " << edmColl->size()
                  << ", LCIO: " << lcColl->getNumberOfElements() << std::endl;
        return 1;
      }
    } catch (EVENT::DataNotAvailableException& ex) {
      std::cerr << "Collection " << name << " is not available in LCEvent" << std::endl;
      return 1;
    }
  }

  const auto objectMapping = ObjectMappings::fromEvent(lcioEvent.get(), edmEvent);

  for (const auto& name : edmEvent.getAvailableCollections()) {
    const auto type = edmEvent.get(name)->getTypeName();
    if (type == "edm4hep::CaloHitContributionCollection" || type == "edm4hep::ParticleIDCollection" ||
        type == "edm4hep::EventHeaderCollection" || type == "edm4hep::RecDqdxCollection") {
      continue;
    }

    if (type == "edm4hep::VertexRecoParticleLinkCollection") {
      // TODO: Dispatch to correct comparison function depending on name (which
      // in turn again depends on a convention in the conversion)

      continue; // avoid exception from LCIO because this is not available there
    }

    const auto* lcioColl = lcioEvent->getCollection(name);

    ASSERT_COMPARE_OR_EXIT(edm4hep::MCParticleCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::ReconstructedParticleCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::TrackCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::TrackerHit3DCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::TrackerHitPlaneCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::SimTrackerHitCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::CalorimeterHitCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::RawCalorimeterHitCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::SimCalorimeterHitCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::RawTimeSeriesCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::ClusterCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::VertexCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::RecoMCParticleLinkCollection)
    ASSERT_COMPARE_OR_EXIT(edm4hep::CaloHitSimCaloHitLinkCollection)

    // TODO: start vertex association (EDM4hep) vs getStartVertex in LCIO
  }

  return 0;
}
