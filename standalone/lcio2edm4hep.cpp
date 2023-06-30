#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

#include <EVENT/LCIO.h>
#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <UTIL/CheckCollections.h>

#include "podio/ROOTFrameWriter.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
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
    std::stringstream sline(std::move(line));
    std::string name, type;
    // This only looks for the first two words in the line and ignores everything that comes after that.
    if (!(sline >> name >> type)) {
      std::cerr << "need a name and a type per line" << std::endl;
      return {};
    }
    names_types.emplace_back(std::move(name), std::move(type));
  }

  input_file.close();

  return names_types;
}

constexpr auto usageMsg = R"(usage: lcio2edm4hep [-h] inputfile outputfile [colltypefile] [-s] [n])";
constexpr auto helpMsg = R"(
Convert an LCIO file to EDM4hep

positional arguments:
  inputfile         the input LCIO file
  outputfile        the output EDM4hep file that will be created
  colltypefile      an optional input file that specifies the names and types of collections
  -s                if -s is appended, only the collections in 'colltypefile' are written
                    otherwise all the specified collections are written
  n                 optionally limit the number of events to n (only w/ -s option)

optional arguments:
  -h, --help        show this help message and exit

Examples:
- print this message:
lcio2edm4hep -h
- convert complete file (needs all lcio collections to be present in all files):
lcio2edm4hep infile.slcio outfile_edm4hep.root
- the same but providing complete set of collections:
lcio2edm4hep infile.slcio outfile_edm4hep.root coltype.txt
- write only 42 events with the subset of collections specified in coltype.txt:
lcio2edm4hep infile.slcio outfile_edm4hep.root coltype.txt -s 42

)";

int main(int argc, char* argv[])
{
  if (argc == 2 && (argv[1] == std::string("-h") || argv[1] == std::string("--help"))) {
    std::cerr << usageMsg << '\n' << helpMsg << std::endl;
    return 0;
  }

  bool createSubset = false;

  if (argc < 3 || (argc > 4 && (argv[4] != std::string("-s")))) {
    std::cerr << usageMsg << std::endl;
    return 1;
  }
  const auto outputFile = std::string(argv[2]);

  bool patching = false;
  UTIL::CheckCollections colPatcher {};
  std::vector<std::pair<std::string, std::string>> namesTypes;
  std::set<std::string> subsetColls;
  int maxEvt = -1;

  if (argc >= 4) {
    namesTypes = getNamesAndTypes(argv[3]);
    if (namesTypes.empty()) {
      std::cerr << "The provided list of collection names and types does not satisfy the required format: Pair of Name "
                   "and Type per line separated by space"
                << std::endl;
      return 1;
    }
    colPatcher.addPatchCollections(namesTypes);
    patching = true;

    if (argc >= 5 && (argv[4] == std::string("-s"))) {
      // new mode: write only a subset of collection as specified in patch file
      patching = false;
      createSubset = true;
      for (auto& nt : namesTypes) {
        subsetColls.insert(nt.first);
      }
    }
    if (argc == 6) maxEvt = std::atoi(argv[5]);
  }

  auto lcreader = IOIMPL::LCFactory::getInstance()->createLCReader();
  lcreader->open(argv[1]);
  std::cout << "Number of events: " << lcreader->getNumberOfEvents() << '\n';
  std::cout << "Number of runs: " << lcreader->getNumberOfRuns() << '\n';

  podio::ROOTFrameWriter writer(outputFile);

  for (auto j = 0u; j < lcreader->getNumberOfRuns(); ++j) {
    if (j % 1 == 0) {
      std::cout << "processing RunHeader: " << j << std::endl;
    }
    auto rhead = lcreader->readNextRunHeader();

    const auto edmRunHeader = LCIO2EDM4hepConv::convertRunHeader(rhead);
    writer.writeFrame(edmRunHeader, "runs");
  }

  int nEvt = maxEvt > 0 ? maxEvt : lcreader->getNumberOfEvents();

  for (auto i = 0u; i < nEvt; ++i) {
    if (i % 10 == 0) {
      std::cout << "processing Event: " << i << std::endl;
    }
    auto evt = lcreader->readNextEvent(EVENT::LCIO::UPDATE);
    // Patching the Event to make sure all events contain the same Collections.
    if (patching == true) {
      colPatcher.patchCollections(evt);
    }
    if (createSubset) {
      auto* colNames = evt->getCollectionNames();
      for (auto& n : *colNames) {
        if (subsetColls.find(n) == subsetColls.end()) evt->removeCollection(n);
      }
    }

    const auto edmEvent = LCIO2EDM4hepConv::convertEvent(evt);
    writer.writeFrame(edmEvent, "events");
  }

  writer.finish();

  return 0;
}
