# v00-08-01

* 2024-02-07 jmcarcell ([PR#48](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/48))
  - Delete build workflow since we have another one for key4hep that covers builds for nightlies, releases and all the operating systems we support

* 2024-02-07 jmcarcell ([PR#47](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/47))
  - Change ROOTFrame{Writer,Reader} to ROOT{Writer,Reader} following Remove the Frame from the default readers and writers following https://github.com/AIDASoft/podio/pull/549

* 2024-02-07 jmcarcell ([PR#46](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/46))
  - Fix compiler warnings related to double - float that appear after https://github.com/key4hep/EDM4hep/pull/237
  - Switch to non-deprecated name for `ROOTWriter` (https://github.com/AIDASoft/podio/pull/549)

# v00-08

* 2024-01-16 jmcarcell ([PR#44](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/44))
  - Change ${LCIO_LIBRARIES} to LCIO::lcio
  - Remove unnecessary LCIO_INCLUDE_DIRS

* 2023-12-03 Mateusz Jakub Fila ([PR#43](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/43))
  - Fix typos in documentation

* 2023-11-30 tmadlener ([PR#42](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/42))
  - Add conversion of TrackerHitPlane from EDM4hep to LCIO.
    - NOTE: The covariance matrix is not set, because there is no public setter available to do so in LCIO.

* 2023-11-28 jmcarcell ([PR#40](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/40))
  - Fix downloading input data for tests (it wasn't happening, at the very least not always).
  - Also compare for NaN values, which can be present; there is at least one value in one of the current test files
  - Make sure that tests are run in CI

# v00-07

* 2023-11-14 tmadlener ([PR#38](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/38))
  - Remove compatibility with all legacy versions of EDM4hep since #34 made the minimum version 0.10.1 in any case.

* 2023-11-08 jmcarcell ([PR#34](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/34))
  - Use the new `edm4hep::CellIDEncoding` for consistency. Needs https://github.com/key4hep/EDM4hep/pull/234

# v00-06-01

* 2023-11-06 tmadlener ([PR#36](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/36))
  - Make sure to convert the full content of `edm4hep::Clusters` to LCIO clusters, including `subdetectorEnergies` and the related calorimeter hits. Set the `contribution` 1.0 because that seems to be the only value in use.
  - Add tests that cover this part of the converter.

* 2023-11-02 jmcarcell ([PR#33](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/33))
  - Use ExternalData to download the test data at build time. See https://github.com/AIDASoft/podio/pull/508

# v00-06

* 2023-10-19 jmcarcell ([PR#31](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/31))
  - Do not forward declare `podio::ObjectID` since this doesn't always build
  - Add a library alias for the `k4EDM4hep2LcioConv` target
  - Add `podio` to the list of required packages (it's already there in the spack recipe)

* 2023-10-05 tmadlener ([PR#29](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/29))
  - Add the existing EDM4hep to LCIO conversion tests from MarlinWrapper
  - Fix minor issues in conversion that were uncovered during this

* 2023-10-05 tmadlener ([PR#21](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/21))
  - Generalize the conversion functionality to make it possible to use "generic maps" (i.e. `vector<tuple<K, V>>` or a proper `map<K, V>`).
    - This makes most of the conversion a template (i.e. header) library.
    - Keep the current behavior by specifying suitable defaults for all the templates.
  - This is necessary to support the introduction of a shared global map in [key4hep/k4MarlinWrapper#147](https://github.com/key4hep/k4MarlinWrapper/pull/147)

* 2023-10-04 tmadlener ([PR#28](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/28))
  - Add test setup to run the standalone converter on ILD REC and DST files with a comparison between the original and the converted file afterwards
  - Make sure to check relations in converted objects
  - Fix a bug in the relation resolution of the ReconstructedParticle that was uncovered.
  - Make `RelWithDebInfo` the default build type.

* 2023-09-12 tmadlener ([PR#27](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/27))
  - Introduce the `EDM4hep2LCIOConv` namespace for the EDM4hep to LCIO conversion functionality to avoid polluting the global namespace with too many (rather generically named) symbols.
  - Define a `EDM4HEP2LCIOCONV_NAMESPACE` preprocessor "symbol" that allows downstream users to make for a slightly smoother transition.

* 2023-07-31 jmcarcell ([PR#26](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/26))
  - `CMAKE_PROJECT_NAME` to `PROJECT_NAME` since `CMAKE_PROJECT_NAME` is the name of the top-level project and `PROJECT_NAME` is the name of the current project.

* 2023-07-31 Thomas Madlener ([PR#25](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/25))
  - Make an output message less confusing

# v00-05

* 2023-07-10 tmadlener ([PR#23](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/23))
  - Remove the explicit `clang-format-check` workflow as it is also covered by the `pre-commit` workflow.

* 2023-07-10 Frank Gaede ([PR#20](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/20))
  - Add the possibility to convert only a subset of collections and events.
  - Fix minor bug in TrackState conversion (covMatrix[15])
  - Write `AllCaloHitContributionsCombined` only if needed, i.e. `SimCalorimeterHits` are present in the events in lcio2edm4hep

* 2023-06-13 Finn Johannsen ([PR#11](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/11))
  - Add LCIO to EDM4hep conversion functionality with a similar interface as the one that is already present for the other direction.

* 2023-06-07 Thomas Madlener ([PR#16](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/16))
  - Match the renaming of `subdetectorHitNumbers` in EDM4hep (cf. [key4hep/EDM4hep#188](https://github.com/key4hep/EDM4hep/pull/188))
  - Introduce some guards for the `TPCHit` -> `RawTimeSeries` change in EDM4hep (cf. [key4hep/EDM4hep#179](https://github.com/key4hep/EDM4hep/pull/179))

* 2023-04-19 Leonhard Reichenbach ([PR#13](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/13))
  - Add support for EventHeader conversion

# v00-04

* 2023-03-02 jmcarcell ([PR#8](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/8))
  - Rename TPCHit -> RawTimeSeries

* 2023-02-22 Thomas Madlener ([PR#9](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/9))
  - Fix the pre-commit workflow and update it to run using a newer LCG release

* 2023-02-22 Finn Johannsen ([PR#4](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/4))
  - Fix the setting of the Track Type in the Tracks Collection in LCIO. It is now set bit wise as required by LCIO

# v00-03

* 2022-11-30 Andre Sailer ([PR#7](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/7))
  - SimCaloHitContribution Conversion: fix an issue when the MCParticles were not converted before the SimCalorimeterHits, which is happening now that the map of collections is sorted in some way according to c++ implementation. Belongs to key4hep/k4MarlinWrapper#99 Fixes key4hep/k4MarlinWrapper#98

* 2022-10-06 Andre Sailer ([PR#6](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/6))
  - CI: add pre-commit workflow with clang-format and whitespace checker

* 2022-10-06 Andre Sailer ([PR#5](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/5))
  - CI: add test for compiling this project

# v00-02-01

* 2022-06-16 Thomas Madlener ([PR#3](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/3))
  - Update cov matrix type after key4hep/EDM4hep#138

# v00-02

* 2022-05-30 Valentin Volkl ([PR#2](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/2))
  - HitContributions were removed from Clusters in https://github.com/key4hep/EDM4hep/pull/140

# v00-01

# v00-00

* Add versioning
