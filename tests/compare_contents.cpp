#include "CompareEDM4hepLCIO.h"

#include "podio/ROOTFrameReader.h"
#include "podio/Frame.h"

#include <IOIMPL/LCFactory.h>

#include <iostream>

#define ASSERT_COMPARE_OR_EXIT(collType)                   \
  if (type == #collType) {                                 \
    auto& edmcoll = edmEvent.get<collType>(name);          \
    if (!compare(lcioColl, edmcoll)) {                     \
      std::cerr << "in collection: " << name << std::endl; \
      return 1;                                            \
    }                                                      \
  }

constexpr auto usageMsg = R"(usage: compare-contents lciofile edm4hepfile)";


int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << usageMsg << std::endl;
    return 1;
  }

  auto lcreader = IOIMPL::LCFactory::getInstance()->createLCReader();
  lcreader->open(argv[1]);

  auto edmreader = podio::ROOTFrameReader();
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

    if (!compareEventHeader(lcEvent,&edmEvent)){
      return 1;
    }

    for (const auto& name : *(lcEvent->getCollectionNames())) {
      const auto lcioColl = lcEvent->getCollection(name);
      // TODO: The Frame needs to improve here in order to get to the type
      // without retrieving the collection
      if (lcioColl->getTypeName() == "LCRelation"){

        const auto& params = lcioColl->getParameters();
        const auto& fromType = params.getStringVal("FromType");
        if (fromType.length() == 0){
          //std::cout<<"WARNING: LCRelations "<< name <<" has no 'to' or 'from' set!"<< std::endl;
          continue;
       }
      }
      const auto& type = [&edmEvent, &name]() {
        const auto coll = edmEvent.get(name);
        if (coll) {
          return coll->getTypeName();
        }
        static std::string empty = "";
        return empty;
      }();
      if (type.empty()) {
        std::cerr << "Collection " << name << " not present in edm4hep file" << std::endl;
        return 1;
      }

      ASSERT_COMPARE_OR_EXIT(edm4hep::MCParticleCollection)
      ASSERT_COMPARE_OR_EXIT(edm4hep::ReconstructedParticleCollection)
      ASSERT_COMPARE_OR_EXIT(edm4hep::TrackCollection)
      ASSERT_COMPARE_OR_EXIT(edm4hep::TrackerHitCollection)
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
