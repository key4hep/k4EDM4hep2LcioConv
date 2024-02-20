#include "CompareEDM4hepLCIO.h"
#include "ObjectMapping.h"
#include "ComparisonUtils.h"

#include "podio/podioVersion.h"
#if PODIO_BUILD_VERSION >= PODIO_VERSION(0, 99, 0)
#include "podio/ROOTReader.h"
#else
#include "podio/ROOTFrameReader.h"
namespace podio {
  using ROOTReader = podio::ROOTFrameReader;
}
#endif

#include "podio/Frame.h"

#include <IOIMPL/LCFactory.h>

#include <iostream>
#include <string_view>

constexpr auto usageMsg = R"(usage: compare-contents lciofile edm4hepfile)";

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << usageMsg << std::endl;
    return 1;
  }

  auto lcreader = IOIMPL::LCFactory::getInstance()->createLCReader();
  lcreader->open(argv[1]);

  auto edmreader = podio::ROOTReader();
  edmreader.openFile(argv[2]);

  // loop going over every name of the lcio collection and checking if a
  // collection with the same name exists in the edm file.
  std::cout << "number of Events" << edmreader.getEntries("events") << std::endl;
  for (size_t n = 0; n < edmreader.getEntries("events"); ++n) {
    if (n % 10 == 0) {
      std::cout << "Event number: " << n << std::endl;
    }

    auto lcEvent = lcreader->readNextEvent();
    auto edmEvent = podio::Frame(edmreader.readNextEntry("events"));

    if (!compareEventHeader(lcEvent, &edmEvent)) {
      return 1;
    }

    for (const auto& name : *(lcEvent->getCollectionNames())) {
      const auto lcioColl = lcEvent->getCollection(name);
      // TODO: The Frame needs to improve here in order to get to the type
      // without retrieving the collection
      if (lcioColl->getTypeName() == "LCRelation") {

        const auto& params = lcioColl->getParameters();
        const auto& fromType = params.getStringVal("FromType");
        if (fromType.length() == 0) {
          // std::cout<<"WARNING: LCRelations "<< name <<" has no 'to' or 'from' set!"<< std::endl;
          continue;
        }
      }
      const auto coll = edmEvent.get(name);
      if (!coll) {
        std::cerr << "Collection " << name << " not present in edm4hep file" << std::endl;
        return 1;
      }

      if (edmEvent.get(name)->size() != (unsigned) lcioColl->getNumberOfElements()) {
        std::cerr << "Collection " << name << " has different sizes. LCIO: " << lcioColl->getNumberOfElements()
                  << ", EDM4hep: " << coll->size() << std::endl;
        return 1;
      }
    }

    const auto objectMapping = ObjectMappings::fromEvent(lcEvent, edmEvent);

    for (const auto& name : *(lcEvent->getCollectionNames())) {
      const auto lcioColl = lcEvent->getCollection(name);
      const auto type = [&edmEvent, &name]() {
        const auto coll = edmEvent.get(name);
        if (coll) {
          return coll->getTypeName();
        }
        static const decltype(coll->getTypeName()) empty = "";
        return empty;
      }();

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
    }
  }

  return 0;
}
