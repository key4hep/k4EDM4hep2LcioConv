# k4EDM4hep2LcioConv

[![Key4hep build](https://github.com/key4hep/k4EDM4hep2LcioConv/actions/workflows/key4hep-build.yaml/badge.svg)](https://github.com/key4hep/k4EDM4hep2LcioConv/actions/workflows/key4hep-build.yaml)
[![DOI](https://zenodo.org/badge/478554694.svg)](https://zenodo.org/doi/10.5281/zenodo.13837370)

<p align="center">
  <img src="doc/k4EDM4hep2LcioConv_logo.svg"/>
</p>


Converter library to convert between the EDM4hep and LCIO event data models.
Supports in-memory conversion in both directions and provides a standalone
conversion tool for LCIO to EDM4hep.

## Dependencies
- LCIO >= v02-22
- EMD4hep >= v00-99
- podio >= v01-00
- ROOT

## Build and install

If you have an environment that fulfils all dependencies (e.g. a Key4hep stack), simply do

- Get the sources to build from
```bash
git clone https://github.com/key4hep/k4EDM4hep2LcioConv
cd k4EDM4hep2LcioConv
```
- Run CMake and configure it to use `install` in the current directory as install prefix
```bash
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=$(pwd)/install
```
- Build the library and tests and run the tests
```bash
cmake --build build
ctest --test-dir build
```
- Install the library and the standalone conversion tools
```bash
cmake --build build --target install
```
