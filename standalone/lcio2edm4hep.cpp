#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

#include <IO/LCReader.h>
#include <IOIMPL/LCFactory.h>
#include <UTIL/CheckCollections.h>

#include "podio/ROOTFrameWriter.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>

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
  std::string inputFile {};
  std::string outputFile {};
  std::string patchFile {};
  int nEvents {-1};
};

void printUsageAndExit()
{
  std::cerr << usageMsg << std::endl;
  std::exit(1);
}

ParsedArgs parseArgs(std::vector<std::string> argv)
{
  // find help
  if (std::find_if(argv.begin(), argv.end(), [](const auto& elem) {
        return elem == "-h" || elem == "--help";
      }) != argv.end()) {
    std::cerr << usageMsg << '\n' << helpMsg << std::endl;
    std::exit(0);
  }

  auto argc = argv.size();
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

int main(int argc, char* argv[])
{
  const auto args = parseArgs({argv, argv + argc});

  UTIL::CheckCollections colPatcher {};
  std::vector<std::pair<std::string, std::string>> namesTypes {};
  const bool patching = !args.patchFile.empty();
  if (patching) {
    namesTypes = getNamesAndTypes(args.patchFile);
    if (namesTypes.empty()) {
      std::cerr << "The provided list of collection names and types does not satisfy the required format: Pair of Name "
                   "and Type per line separated by space"
                << std::endl;
      return 1;
    }
    colPatcher.addPatchCollections(namesTypes);
  }
  // Construct a vector of collections to convert. If namesTypes is empty, this
  // will be empty, and convertEvent will fall back to use the collections in
  // the event
  const auto collsToConvert = [&namesTypes]() {
    std::vector<std::string> names;
    names.reserve(namesTypes.size());
    for (const auto& [name, type] : namesTypes) {
      names.emplace_back(name);
    }
    return names;
  }();

  auto lcreader = IOIMPL::LCFactory::getInstance()->createLCReader();
  lcreader->open(args.inputFile);
  std::cout << "Number of events in file: " << lcreader->getNumberOfEvents() << '\n';
  std::cout << "Number of runs in file: " << lcreader->getNumberOfRuns() << '\n';

  podio::ROOTFrameWriter writer(args.outputFile);

  for (auto j = 0u; j < lcreader->getNumberOfRuns(); ++j) {
    if (j % 1 == 0) {
      std::cout << "processing RunHeader: " << j << std::endl;
    }
    auto rhead = lcreader->readNextRunHeader();

    const auto edmRunHeader = LCIO2EDM4hepConv::convertRunHeader(rhead);
    writer.writeFrame(edmRunHeader, "runs");
  }

  const int nEvt = args.nEvents > 0 ? args.nEvents : lcreader->getNumberOfEvents();
  for (auto i = 0u; i < nEvt; ++i) {
    if (i % 10 == 0) {
      std::cout << "processing Event: " << i << std::endl;
    }
    auto evt = lcreader->readNextEvent();
    // Patching the Event to make sure all events contain the same Collections.
    if (patching == true) {
      colPatcher.patchCollections(evt);
    }
    const auto edmEvent = LCIO2EDM4hepConv::convertEvent(evt, collsToConvert);
    writer.writeFrame(edmEvent, "events");
  }

  writer.finish();

  return 0;
}
