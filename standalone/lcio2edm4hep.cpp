#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <UTIL/CheckCollections.h>

#include "podio/ROOTFrameWriter.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

std::vector<std::pair<std::string, std::string>> getNamesAndTypes(const std::string& collTypeFile)
{
  std::ifstream input_file(collTypeFile);
  std::vector<std::pair<std::string, std::string>> names_types;

  if (!input_file.is_open()) {
    std::cerr << "Failed to open file countaining the names and types of the LCIO Collections." << std::endl;
  }
  std::string line;
  while (std::getline(input_file, line)) {
    size_t delimiter_pos = line.find(' ');
    if (delimiter_pos == std::string::npos) {
      std::cerr << "Invalid input format" << std::endl;
      input_file.close();
      return {};
    }
    std::string name = line.substr(0, delimiter_pos);
    std::string type = line.substr(delimiter_pos + 1);
    names_types.emplace_back(name, type);
  }

  input_file.close();

  return names_types;
}

constexpr auto usageMsg = R"(usage: lcio2edm4hep [-h] inputfile outputfile [colltypefile])";
constexpr auto helpMsg = R"(
Convert an LCIO file to EDM4hep

positional arguments:
  inputfile         the input LCIO file
  outputfile        the output EDM4hep file that will be created
  colltypefile      An optional input file that specifies the names and types of all
                    collections that should be present in the output.

optional arguments:
  -h, --help        show this help message and exit
)";

int main(int argc, char* argv[])
{
  if (argc == 2 && (argv[1] == std::string("-h") || argv[1] == std::string("--help"))) {
    std::cerr << usageMsg << '\n' << helpMsg << std::endl;
    return 0;
  }
  if (argc < 3) {
    std::cerr << usageMsg << std::endl;
    return 1;
  }
  const auto outputFile = std::string(argv[2]);

  bool patching = false;
  UTIL::CheckCollections colPatcher {};
  if (argc == 4) {
    const auto namesTypes = getNamesAndTypes(argv[3]);
    if (namesTypes.empty()) {
      std::cerr << "The provided list of collection names and types does not satisfy the required format: Pair of Name "
                   "and Type per line separated by space"
                << std::endl;
      return 1;
    }
    colPatcher.addPatchCollections(namesTypes);
    patching = true;
  }

  auto lcreader = IOIMPL::LCFactory::getInstance()->createLCReader();
  lcreader->open(argv[1]);
  std::cout << "Number of events: " << lcreader->getNumberOfEvents() << '\n';

  podio::ROOTFrameWriter writer(outputFile);

  std::cout << "starting Conversion" << std::endl;
  for (auto i = 0u; i < lcreader->getNumberOfEvents(); ++i) {
  //for (auto i = 0u; i < 4; ++i) {
    if (i % 10 == 0) {
      std::cout << "processing Event: " << i << std::endl;
    }
    auto evt = lcreader->readNextEvent();
    std::cout<<"read Event"<<std::endl;
    // Patching the Event to make sure all events contain the same Collections.
    if (patching == true) {
      std::cout<<"starting to patch"<<std::endl;
      colPatcher.patchCollections(evt);
      std::cout<<"patched"<<std::endl;
    }
    const auto edmEvent = LCIO2EDM4hepConv::convertEvent(evt);
    writer.writeFrame(edmEvent, "events");
  }

  writer.finish();

  return 0;
}
