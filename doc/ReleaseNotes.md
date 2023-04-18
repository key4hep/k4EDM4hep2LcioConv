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
