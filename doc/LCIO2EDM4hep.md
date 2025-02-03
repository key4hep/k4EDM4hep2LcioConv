# Standalone conversion from LCIO to EDM4hep
The `lcio2edm4hep` executable reads LCIO (`.slcio`) files and converts its
contents into EDM4hep. Each `LCEvent` of the input file will be put into a
`podio::Frame` in the output file (under the `events` category). The most basic
usage is simply

```bash
lcio2edm4hep <input.slcio> <output.edm4hep.root>
```

## Patching missing collections on the fly
A major difference between LCIO and EDM4hep is that in LCIO an `LCEvent` can
effectively have arbitrary contents, whereas in EDM4hep the assumption is that
each event consists of the same collections (even if some of them are empty).
Hence, it is necessary to either ensure that all events in the LCIO file have
the same contents or or to give `lcio2edm4hep` some additional information such
that it can patch in potentially missing collections on the fly. This additional
information comes in the form of a third argument to `lcio2edm4hep` and is
effectively a list of collection names and their types that comprise the
superset of all collections appearing in at least one event in the input LCIO
file. The format looks like this

```
name[:output-name]  type-name[*]
```

Each collection is a single line containing the name first (including an
optional output name [see below](#renaming-collections-on-the-fly)) and than its
type. An additional `*` signifies a *subset collection*. The simplest form looks like this:

```
SETSpacePoints             TrackerHit
RecoMCTruthLink            LCRelation[ReconstructedParticle,MCParticle]
```

The easiest way to obtain such a file is to use the `check_missing_cols`
executable that comes with LCIO using the `--minimal` flag. The output of this
can be directly consumed by `lcio2edm4hep`

### Patching missing ParticleID information on the fly

EDM4hep also assumes that the `ParticleID` objects that are attached to elements
of a `ReconstructedParticle` are consistent, i.e. the same PID algo names have
been used throughout the processing. In order to guarantee this for the
conversion it is possible to attach missing information on the fly, the grammar
for this is

```
pid-algo-name  reco-coll-name|[parameter-names[,param-names]]
```

This will use the (LCIO) `PIDHandler` to add a PID algorithm with name
`pid-algo-name` to the collection `reco-coll-name`. Optionally if any parameter
names are present it will also set the parameter names for this PID algorithm.

```{note}
This is only available from LCIO versions **larger** than `v02-22-01`!
```

### Patching `LCRelation` collections
For collections of `LCRelation` type it is necessary to define the `FromType` and
`ToType` as well, as otherwise the converter will not be able to create the
correct edm4hep file. The `check_missing_cols` executable will try to determine
these types from the collection parameters and will warn if it cannot do it for
certain collections. In this case it is the **users responsibility to provide
the missing types** as otherwise the conversion will simply skip these
collections, or potentially even crash.

### Patching collections as subset collections
A single `*` (star) after the type name tells the patching to flip the *subset*
collection flag for a newly created collection. This can be necessary for
producing consistent outputs in EDM4hep. As an example

```
FTD_STRIPCollection      SimTrackerHit*
```
will create a `SimTrackerHit` subset collection with the name `FTD_STRIPCollection`.


```{note}
This is only available from LCIO versions **larger** than `v02-22-03`!
```

#### Example:
1. Get the patch file
```bash
check_missing_cols --minimal \
  /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/higgs/ILD_l5_o2_v02/v02-02-01/00015671/000/rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I402005.Pe3e3h.eL.pR.n000_002.d_rec_00015671_493.slcio \
  > patch.txt
```
2. Pass it to `lcio2edm4hep`
```bash
lcio2edm4hep \
  /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/higgs/ILD_l5_o2_v02/v02-02-01/00015671/000/rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I402005.Pe3e3h.eL.pR.n000_002.d_rec_00015671_493.slcio \
  Output.root \
  patch.txt
```


## Converting only a subset of collections
Using the same mechanism as for patching collections it is also possible to only
convert a subset of all available collections. `lcio2edm4hep` uses the contents
of the `colltypefile` to determine the contents of the output. If that contains
only a subset of all collections, only that subset will be converted. Missing
collections will still be patched in, in this case.

## Renaming collections on the fly
The optional `[:output-name]` part of each collection can be used to remap the
names of the collections in the input LCIO file to a different name in the
output EDM4hep file, e.g.

```
MCParticle:MCParticles      MCParticle
```

will read the `MCParticle` collection from the input file but store it as
`MCParticles` in the output file.

# Library usage of the conversion functions
The conversion functions are designed to also be usable as a library. The overall design is to make the conversion a two step process. Step one is converting the data and step two being the resolving of the relations and filling of subset collection.

## Converting collection (data)
The main entry point is `convertCollection` which will automatically dispatch to
the correct conversion function depending on the type information that is stored
in the input `LCCollection`. It is also possible to access the individual
conversion functions for each type. All of the conversion functions take a map
of LCIO to EDM4hep objects of their specific type that will be filled during the
conversion. for convenience all necessary maps are bundled in the
`LcioEdmTypeMapping` struct.

## Handling relations
**Once all necessary collections have been converted, it is necessary to resolve
the relations between the objects.** This is done using the `resolveRelations`
function. This will again dispatch to the correct relation resolving function
for the corresponding types, which can obviously also be invoked directly.

## Handling of subset collections
Subset collections are handled similar to relations using the function
`fillSubset`. Internally this simply forwards to `handleSubsetColl` which
handles all the type details and can obviously also be used directly.

## Handling of `LCRelation`s
`LCRelation` only exist in LCIO and their conversion is limited to what is
available in EDM4hep. They use the `"FromType"` and `"ToType"` collection
parameters to get the necessary type information.

The LinkCollections in EDM4hep are then created using `createLinks`.

## Converting entire events
Converting an entire event can be done calling the `convertEvent`. This can also
be used as an example to guide the implementation of custom conversions using
the available functionality.

## Converting Event parameters
This can be done by calling `convertObjectParameters` that will put all the event parameters into the passed `podio::Frame`.

## Subtle differences between LCIO and EDM4hep
There are a few small differences between LCIO and EDM4hep that shine through in the conversion, these are:

- `CaloHitContributions` are part of the SimCalorimeterHits in LCIO while being their own data type in EDM4hep. They are created by [`createCaloHitContributions`](https://github.com/key4hep/k4EDM4hep2LcioConv/blob/main/k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h).
- The event information like an event number is part of the `LCEvent` in LCIO. In EDM4hep there is a separate EventHeader Collection. It can be created using [`EventHeaderCollection`](https://github.com/key4hep/k4EDM4hep2LcioConv/blob/main/k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h) which is stored under the name `"EventHeader"`.
- Particle IDs are converted during the conversion of the the reconstructed Particle collection.

## Example for a ReconstructedParticle Collection
```cpp
#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"

// the struct defined in the header file is used for the maps linking Lcio particles
// to their EDM counterparts.

auto typeMapping = LcioEdmTypeMapping{};

// We assume that this is a collection of ReconstructedParticles!
LCEVENT::LCCollection* lcCollection;

// Convert the data
auto edmCollections = convertReconstructedParticle("name",
                                                   lcCollection,
                                                   typeMapping.recoParticles,
                                                   typeMapping.particleIDs);

// Resolve relations (only converted objects will be available)
// This has to be called at the very end, after all collection data has been
// converted
resolveRelations(typeMapping);
```
