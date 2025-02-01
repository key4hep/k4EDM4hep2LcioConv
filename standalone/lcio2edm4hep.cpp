#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <UTIL/CheckCollections.h>

#include <edm4hep/utils/ParticleIDUtils.h>

#include "podio/ROOTWriter.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <optional>
#include <string>
#include <utility>
#include <vector>

/// Simple helper struct to group information about a collection to be converted
struct NamesType {
  std::string lcioName{};    ///< The name in the lcio file
  std::string edm4hepName{}; ///< The name in the edm4hep file
  std::string typeName{};    ///< The LCIO type
};

/// Convert a config file line into a NamesType struct
std::optional<NamesType> fromConfigLine(std::string line) {
  NamesType info;
  std::stringstream sline(std::move(line));
  std::string name;
  // This only looks for the first two words in the line and ignores
  // everything that comes after that.
  if (!(sline >> name >> info.typeName)) {
    std::cerr << "need a name (mapping) and a type per line" << std::endl;
    return std::nullopt;
  }

  if (const auto colon = name.find(':'); colon != std::string::npos) {
    info.lcioName = name.substr(0, colon);
    info.edm4hepName = name.substr(colon + 1);
  } else {
    info.lcioName = name;
    info.edm4hepName = name;
  }

  return info;
}

std::vector<NamesType> getNamesAndTypes(const std::string& collTypeFile) {
  std::ifstream input_file(collTypeFile);
  std::vector<NamesType> names_types;

  if (!input_file.is_open()) {
    std::cerr << "Failed to open file containing the names and types of the LCIO Collections." << std::endl;
  }
  std::string line;
  while (std::getline(input_file, line)) {
    auto lineInfo = fromConfigLine(std::move(line));
    if (!lineInfo) {
      return {};
    }
    names_types.emplace_back(std::move(lineInfo.value()));
  }
  input_file.close();

  return names_types;
}

constexpr auto usageMsg = R"(usage: lcio2edm4hep [-h] inputfile outputfile [colltypefile] [-n N])";

constexpr auto helpMsg = R"(
Convert an LCIO file to EDM4hep

positional arguments:
  inputfile         the input LCIO file
  outputfile        the output EDM4hep file that will be created
  colltypefile      An optional input file that specifies the names and types of
                    collections that should be present in the output.

optional arguments:
  -h, --help        show this help message and exit
  -n N              Limit the number of events to convert to N (default = -1, all events)

Examples:
- print this message:
lcio2edm4hep -h
- convert complete file (needs all lcio collections to be present in all files):
lcio2edm4hep infile.slcio outfile_edm4hep.root
- the same but providing complete set of collections (either to patch collections in,
  or to only convert a subset):
lcio2edm4hep infile.slcio outfile_edm4hep.root coltype.txt
)";

struct ParsedArgs {
  std::string inputFile{};
  std::string outputFile{};
  std::string patchFile{};
  int nEvents{-1};
};

void printUsageAndExit() {
  std::cerr << usageMsg << std::endl;
  std::exit(1);
}

ParsedArgs parseArgs(std::vector<std::string> argv) {
  // find help
  if (std::find_if(argv.begin(), argv.end(), [](const auto& elem) { return elem == "-h" || elem == "--help"; }) !=
      argv.end()) {
    std::cerr << usageMsg << '\n' << helpMsg << std::endl;
    std::exit(0);
  }

  int argc = argv.size();
  if (argc < 3 || argc > 7) {
    printUsageAndExit();
  }

  ParsedArgs args;
  auto nEventIt = std::find(argv.begin(), argv.end(), "-n");
  if (nEventIt != argv.end()) {
    const auto index = std::distance(argv.begin(), nEventIt);
    if (index > argc - 2) {
      // No argument left to parse
      printUsageAndExit();
    }
    const auto& value = argv[index + 1]; // get the actual value
    try {
      args.nEvents = std::stoi(value);
    } catch (std::invalid_argument& err) {
      std::cerr << "Cannot parse " << value << " as an integer" << std::endl;
      printUsageAndExit();
    }
    argv.erase(nEventIt, nEventIt + 2);
  }

  argc = argv.size();
  if (argc < 3) {
    printUsageAndExit();
  }
  args.inputFile = argv[1];
  args.outputFile = argv[2];
  if (argc == 4) {
    args.patchFile = argv[3];
  }
  return args;
}

int main(int argc, char* argv[]) {
  const auto args = parseArgs({argv, argv + argc});

  UTIL::CheckCollections colPatcher{};
  std::vector<NamesType> namesTypes{};
  const bool patching = !args.patchFile.empty();
  if (patching) {
    namesTypes = getNamesAndTypes(args.patchFile);
    if (namesTypes.empty()) {
      std::cerr << "The provided list of collection names and types does not satisfy the required format: Pair of Name "
                   "(mapping) and Type per line separated by space"
                << std::endl;
      return 1;
    }
    std::vector<std::pair<std::string, std::string>> patchNamesTypes{};
    patchNamesTypes.reserve(namesTypes.size());
    for (const auto& [name, _, type] : namesTypes) {
      patchNamesTypes.emplace_back(name, type);
    }
    colPatcher.addPatchCollections(patchNamesTypes);
  }
  // Construct a vector of collections to convert. If namesTypes is empty, this
  // will be empty, and convertEvent will fall back to use the collections in
  // the event
  const auto collsToConvert = [&namesTypes]() {
    std::vector<std::pair<std::string, std::string>> names{};
    names.reserve(namesTypes.size());
    for (const auto& [lcioName, edm4hepName, type] : namesTypes) {
      // filter out the ParticleID patching from collection names to convert
      if (type.find('|') == std::string::npos) {
        names.emplace_back(lcioName, edm4hepName);
      }
    }
    return names;
  }();

  auto lcreader = IOIMPL::LCFactory::getInstance()->createLCReader();
  lcreader->open(args.inputFile);
  std::cout << "Number of events in file: " << lcreader->getNumberOfEvents() << '\n';
  std::cout << "Number of runs in file: " << lcreader->getNumberOfRuns() << '\n';

  podio::ROOTWriter writer(args.outputFile);

  podio::Frame metadata{};

  for (int j = 0; j < lcreader->getNumberOfRuns(); ++j) {
    if (j % 1 == 0) {
      std::cout << "processing RunHeader: " << j << std::endl;
    }
    auto rhead = lcreader->readNextRunHeader();

    const auto edmRunHeader = LCIO2EDM4hepConv::convertRunHeader(rhead);
    writer.writeFrame(edmRunHeader, "runs");
  }

  const int nEvt = args.nEvents > 0 ? args.nEvents : lcreader->getNumberOfEvents();
  bool haveSimCaloHits{false};
  for (int i = 0; i < nEvt; ++i) {
    int tenPercent = nEvt / 10;
    if ((i + 1) % tenPercent == 0) {
      std::cout << "processed amount of events: " << percEvt << "% (event: " << i << ")" << std::endl;
    }
    auto evt = lcreader->readNextEvent();
    // Patching the Event to make sure all events contain the same Collections.
    if (patching == true) {
      colPatcher.patchCollections(evt);
    }
    if (i == 0) {
      for (const auto& name : *evt->getCollectionNames()) {
        if (evt->getCollection(name)->getTypeName() == "SimCalorimeterHit") {
          haveSimCaloHits = true;
          break;
        }
      }
    }

    auto edmEvent = LCIO2EDM4hepConv::convertEvent(evt, collsToConvert);
    if (haveSimCaloHits && edmEvent.get("AllCaloHitContributionsCombined") == nullptr) {
      edmEvent.put(edm4hep::CaloHitContributionCollection(), "AllCaloHitContributionsCombined");
    }

    // For the first event we also convert some meta information for the
    // ParticleID handling
    if (i == 0) {
      for (const auto& name : *evt->getCollectionNames()) {
        auto coll = evt->getCollection(name);
        if (coll->getTypeName() == "ReconstructedParticle") {
          for (const auto& pidInfo : LCIO2EDM4hepConv::getPIDMetaInfo(coll)) {
            edm4hep::utils::PIDHandler::setAlgoInfo(metadata, LCIO2EDM4hepConv::getPIDCollName(name, pidInfo.algoName),
                                                    pidInfo);
          }
        }
      }
    }

    writer.writeFrame(edmEvent, "events");
  }

  writer.writeFrame(metadata, podio::Category::Metadata);

  writer.finish();

  return 0;
}
