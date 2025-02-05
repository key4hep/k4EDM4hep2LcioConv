# v00-10

* 2025-02-04 Thomas Madlener ([PR#109](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/109))
  - Update the documentation to include the new capabilities of patching subset collections, introduced in [iLCSoft/LCIO#201](https://github.com/iLCSoft/LCIO/pull/201)

* 2025-02-03 jmcarcell ([PR#108](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/108))
  - Fix compilation warning about comparing integers with different signs, related to https://github.com/key4hep/EDM4hep/pull/398

* 2024-12-19 Thomas Madlener ([PR#104](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/104))
  - Remove the conversion of the `colorFlow` from the MCParticle. For LCIO the color flow will be set to `{0, 0}` during the conversion. (`colorFlow` is removed from the `edm4hep::MCParticle` in [EDM4hep#389](https://github.com/key4hep/EDM4hep/pull/389))

* 2024-12-18 Thomas Madlener ([PR#103](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/103))
  - Remove the inclusion of a no longer existing header in the tests.

* 2024-12-17 Mateusz Jakub Fila ([PR#101](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/101))
  - Replaced link to local file with link to github, so the page build correctly in keyhep documentation

* 2024-12-10 Thomas Madlener ([PR#100](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/100))
  - Make the pre-commit CI workflow run on EL9

* 2024-12-10 jmcarcell ([PR#99](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/99))
  - Remove the check for TrackerHit3D from edm4hep
  - Remove the check for `CovMatrix3f.h`

* 2024-09-25 tmadlener ([PR#98](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/98))
  - Update README and remove outdated content
    - Add zenodo and CI status badges
    - Add description of how to build

* 2024-09-12 tmadlener ([PR#97](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/97))
  - Exclude the release notes from pre-commit check for trailing whitespaces

# v00-09

* 2024-09-10 jmcarcell ([PR#96](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/96))
  - Use the Key4hepConfig flag to set the standard, compiler flags and rpath magic.
  - Fix warnings about shadowing variables that were not there before

* 2024-09-09 tmadlener ([PR#91](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/91))
  - Remove deprecated conversion functions now that all downstream consumers have switched.
  - Remove *association* (and derived terms) from documentation, variable names and comments.

* 2024-09-05 tmadlener ([PR#95](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/95))
  - Convert `Nholes` and `subdetectorHoleNumbers` for tracks
  - Bump the minimum LCIO version to `2.22` since this information is not available before that

* 2024-09-05 tmadlener ([PR#94](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/94))
  - Adjust the standalone conversion for the new capabilities of LCIO to also patch ParticleID (meta) information on the fly (see [LCIO#193](https://github.com/iLCSoft/LCIO/pull/193)).

* 2024-08-10 jmcarcell ([PR#89](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/89))
  - Don't throw a pedantic warning because of #warning

* 2024-08-10 jmcarcell ([PR#80](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/80))
  - Delete the version checks for Podio before 1.0

* 2024-08-09 Andre Sailer ([PR#88](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/88))
  - LCIO2EDM4hep collection: fix linking of simtrackerhits and trackerhitplanes. The wrong map was used for associating the converted (edm4hep) trackerhitplanes, since the other map is for TrackerHit3D, fixes https://github.com/key4hep/CLDConfig/issues/48

* 2024-08-02 tmadlener ([PR#87](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/87))
  - Depend on EDM4hep v00-99 in CMake configuration

* 2024-07-31 tmadlener ([PR#86](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/86))
  - Switch from `Association` to the newer `Link` types (key4hep/EDM4hep#341) and use the non-deprecated methods on them
  - Update some conversion functions and docstrings to also reflect this renaming

* 2024-07-25 tmadlener ([PR#84](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/84))
  - Make sure that all the necessary definitions are available
  - Force the evaluation of some string constants to compile time to detect undefined functions earlier

* 2024-07-24 Leonhard Reichenbach ([PR#83](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/83))
  - Fix a typo in the type checks for creating Associations between `TrackerHit`s and `SimTrackerHit`s.

* 2024-07-19 tmadlener ([PR#82](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/82))
  - Adapt the conversion of `Vertex` and `ReconstructedParticle` to follow the new relation schema introduce in [key4hep/EDM4hep#332](https://github.com/key4hep/EDM4hep/pull/332). **NOTE: This is not yet 100 % debugged and tested, due to the inherent complexities of getting both conventions to work simultaneously and together**. Please report any issues that you observe.

* 2024-07-16 tmadlener ([PR#69](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/69))
  - Make the conversion handle the fact that EDM4hep tracks do no longer have any dQ/dx information but rather store this in a separate `RecDqdx` collection.
    - From LCIO to EDM4hep: Create a `RecDqdx` collection for every converted track collection with the suffix `_dQdx`
    - From EDM4hep to LCIO: Offer functionality to attach the information stored in `RecDqdx` collections to converted tracks via the `attachDqdxInfo` function(s).

* 2024-07-12 tmadlener ([PR#81](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/81))
  - Add conversion of EDM4hep `Association` collections to `LCRelation` collections in LCIO.
    - The order of the relation in LCIO follows the following pattern: `From` will be the reconstruction part and `To` will be the simulation / mc part.

* 2024-07-08 tmadlener ([PR#74](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/74))
  - Switch to the new `Vertex::isPrimary` functionality after its introduction in [EDM4hep#329](https://github.com/key4hep/EDM4hep/pull/329)

* 2024-07-04 tmadlener ([PR#79](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/79))
  - Always put a `CaloHitContribution` collection into the event, even if no `SimCalorimeterHit`s have been converted in `lcio2edm4hep`. Fixes [#78](https://github.com/key4hep/k4EDM4hep2LcioConv/issues/78)

* 2024-07-03 jmcarcell ([PR#72](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/72))
  - Don't use radiusOfInnermostHit for EDM4hep tracks, compute it from track state at first hit when going from LCIO to EDM4hep.

* 2024-07-02 tmadlener ([PR#76](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/76))
  - Make sure to still convert assocations between `TrackerHitPlane` and `SimTrackerHit`

* 2024-07-01 jmcarcell ([PR#77](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/77))
  - Find `MathCore` from ROOT since it is being linked to later

* 2024-06-28 jmcarcell ([PR#71](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/71))
  - Use ndf instead of probability for vertexes in EDM4hep
  - Add a utility function to find ndf from chi^2 and the probability to go from LCIO to EDM4hep
  - Add a utility macro to compare float values

* 2024-06-27 tmadlener ([PR#73](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/73))
  - Remove the `MCRecoTrackerHitPlaneAssociation` since it has been / will be removed from EDM4hep in [EDM4hep#331](https://github.com/key4hep/EDM4hep/pull/331)

* 2024-06-20 jmcarcell ([PR#70](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/70))
  - Use edm4hep::labels for `cellIDEncoding`

* 2024-06-11 tmadlener ([PR#68](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/68))
  - Switch to member function access for `algoType` in ParticleID conversion (See https://github.com/key4hep/EDM4hep/pull/307 for more details)

* 2024-06-03 Leonhard Reichenbach ([PR#66](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/66))
  - Standalone printout only for every 10% of processed events instead of every 10 events

* 2024-05-16 tmadlener ([PR#64](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/64))
  - Introduce pre-processor checks to transparently switch to the new `std::optional` return values of `podio::Frame::getParameter` (introduced with [AIDASoft/podio#580](https://github.com/AIDASoft/podio/pull/580))

* 2024-05-16 tmadlener ([PR#63](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/63))
  - Format the `.ipp` files according to the new `.clang-format` configuration and make sure that `pre-commit` enforces it.

* 2024-05-16 tmadlener ([PR#61](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/61))
  - Introduce possibility to remap collection names during the standalone conversion from LCIO to EDM4hep (fixes [#58](https://github.com/key4hep/k4EDM4hep2LcioConv/issues/58))
    - Make the patch file grammar accept an optional `[:output-name]` as part of the collection name that will be used for the output collection

* 2024-05-03 tmadlener ([PR#60](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/60))
  - Synchronize `.clang-format` configuration with the one from `key4hep-dev-utils/defaults` to have the same configuration as in other Key4hep repositories.
  - Make necessary format changes.
  - Switch the pre-commit CI workflow to use the key4hep nightlies as environment.

* 2024-05-03 jmcarcell ([PR#59](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/59))
  - Fix narrowing in Cov3f from int to float. This triggers a warning in GCC 13 and a compiler error in Clang 17. 
  - Fix a compiler error and a warning with Clang. The first one about constness of the LCIO pid handler and the other one about a lambda parameter not being used

* 2024-05-01 tmadlener ([PR#56](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/56))
  - Make the introduction of covariance matrix [key4hep/EDM4hep#287](https://github.com/key4hep/EDM4hep/pull/287) components transparent

* 2024-05-01 tmadlener ([PR#51](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/51))
  - Adapt the conversion of `ReconstructedParticle` and `ParticleID` after the reversal of the relations in https://github.com/key4hep/EDM4hep/pull/268

* 2024-03-27 Leonhard Reichenbach ([PR#57](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/57))
  - Fix typo in run header parameter `detectoName` to `detectorName` in k4Lcio2EDM4hepConv.cpp

* 2024-03-12 tmadlener ([PR#54](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/54))
  - Rename all `convXYZ` methods to `convertXYZ` and deprecate the former versions for the EDM4hep to LCIO direction to make it more consistent with the other direction.
  - Cleanly split the EDM4hep to LCIO conversion into two steps, doing data conversion only in the first step and relation resolving in a second step.
    - Introduce `resolveXYZRelations` functions for all types where it is necessary and remove preliminary relation resolving from data conversion functions
    - All `convertXYZ` functions now only take one map that is populated as they no longer have to do any relation resolving
  - Add roundtrip tests for ReconstructedParticle conversion
    - Fix issue that was lurking here in passing by cleanly splitting conversion process into two steps

* 2024-03-12 tmadlener ([PR#53](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/53))
  - Add compiler warnings to the build process
  - Fix compiler warnings

* 2024-03-11 tmadlener ([PR#55](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/55))
  - Switch to non-deprecated access methods after original methods have been deprecated in [EDM4hep#267](https://github.com/key4hep/EDM4hep/pull/267), [EDM4hep#256](https://github.com/key4hep/EDM4hep/pull/256) and [EDM4hep#273](https://github.com/key4hep/EDM4hep/pull/273)

* 2024-02-23 tmadlener ([PR#49](https://github.com/key4hep/k4EDM4hep2LcioConv/pull/49))
  - Convert all TrackerHit types now that EDM4hep has a TrackerHit interface

# v00-08-02

* 2024-02-07 tmadlener ([PR#50](https://github.com/key4hep/k4edm4hep2lcioconv/pull/50))
  - Make it possible to keep working with MCParticle momenta based on floats for now. See https://github.com/key4hep/EDM4hep/pull/266

* 2024-02-07 tmadlener ([PR#45](https://github.com/key4hep/k4edm4hep2lcioconv/pull/45))
  - Make the upcoming renaming of `TrackerHit` -> `TrackerHit3D` transparent

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
