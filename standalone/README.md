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
