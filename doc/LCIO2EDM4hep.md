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
superset of all collectoins appearing in at least one event in the input LCIO
file. The format looks like this, where each collection is a single line
containing the name first and than its type, e.g.

```
SETSpacePoints             TrackerHit
RecoMCTruthLink            LCRelation[ReconstructedParticle,MCParticle]
```

The easiest way to obtain such a file is to use the `check_missing_cols`
executable that comes with LCIO using the `--minimal` flag. The output of this
can be directly consumed by `lcio2edm4hep`

Example:
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

### Converting `LCRelation` collections
For collections of `LCRelation` type it is necessary to define the `From` and
`To` type as well, as otherwise the converter will not be able to create the
correct edm4hep file. The `check_missing_cols` executable will try to determine
these types from the collection parameters and will warn if it cannot do it for
certain collections. In this case it is the **users responsibility to provide
the missing types** as otherwise the conversion will simply skip these
collections, or potentially even crash.


# Integrated use of conversion
The functions can be used integrated, making the conversion a two step process. Step one is converting the data and step two being the resolving of therelations and filling of subset collection.
There exists a convert function for every collection not of `LCRelation` type (e.g. [`convertReconstructedParticle`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h)). These need to be called before the relations can be handled, since they fill the maps linking the particle in LCIO to their edm4HEP equivalents. Every type has a seperate map. The maps are grouped in a struct for ease of use. They can be defined and stored seperatly.
The order in which the data is converted does not matter because converting data and resolving relations are two differet steps that are carried out in sequence. 
Subset collections are also handled at the same step as relations using the fuction [`fillSubset`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h). Alternatively [`handleSubsetColl`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h) can also be called to convert a Subset Collection. This way the unique pointers can be obtained directly.
The OneToMany and OnToOne Relations can be resolved using [`resolveRelations`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h).  There is a resolveRelations function for each type.
The AssociationCollections in EDM4hep are then created using [`createAssociations`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h), the exception here are the CaloHitContributions. Since they are part of the SimCalorimeterHits in LCIO while being a seperate Association in EDM4hep. They are created by [`createCaloHitContributions`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h).
The EventHeader Colletion can be created using [`EventHeaderCollection`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h).

Particle IDs are converted during the conversion of the the reconstructed Particle collection.


Converting an entire event can be done calling the [`convertEvent`](../k4EDM4hep2LcioConv/include/k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h).

Example for a ReconstructedParticle Collection
```cpp
//the structs defined in the header file are used for the maps linking Lcio particles to there EDM counterparts.

#include "k4EDM4hep2LcioConv/k4Lcio2EDM4hepConv.h"  

auto convertedReconstructedParticleCollection = convertReconstructedParticle(name, LCCollection, typeMapping.recoParticles, typeMapping.particleIDs)

//If the relations to other data types are supposed to be converted it is necessary that these are converted aswell by calling their convert function. This needs to be done in order to fill the maps used in setting the relations.
//For a collction of the type reconstructedparticle those are the vertex, cluster and track collections containing the data related to the reconstructed particles. 

//next step is resolving the relations.
resolveRelationsRecoParticle(
      typeMapping.recoParticles, typeMapping.vertices, typeMapping.clusters, typeMapping.tracks);

//after this the reconstructed particles in `convertedReconstructedParticleCollection` that was created earlier got their vertecies, clusters and tracks attached.
```