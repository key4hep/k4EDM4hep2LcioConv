#include "EDM4hep2LCIOUtilities.h"
#include "CompareEDM4hepLCIO.h"

#include "k4EDM4hep2LcioConv/k4EDM4hep2LcioConv.h"

#include "podio/Frame.h"

#include <type_traits>

template<typename>
struct TD;

int main()
{
  auto convObjMaps = EDM4hep2LCIOConv::CollectionsPairVectors {};
  auto objectMaps = ObjectMappings {};

  const auto caloHitColl = createCalorimeterHits(10);
  const auto lcioColl = EDM4hep2LCIOConv::convCalorimeterHits(&caloHitColl, "", convObjMaps.calohits);

  const auto edmEvent = createExampleEvent();

  if (!compare(lcioColl, caloHitColl, objectMaps)) {
    return 1;
  }

  return 0;
}
