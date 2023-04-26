To convert aa LCIO file to EDM4hep either a LCIO file where all collections are present in every event
or a file containing a list of all collections in the file are needed. 
This list can be obtained by executing 'check_missing_cols' from the LCIO package with the '--minimal' flag.

Example:


```bash
check_missing_cols --minimal \
  /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/higgs/ILD_l5_o2_v02/v02-02-01/00015671/000/rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I402005.Pe3e3h.eL.pR.n000_002.d_rec_00015671_493.slcio \
  > patch.txt
```

This file will be the third argument when the converter ist called. 

```bash
lcio2edm4hep \
  /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/higgs/ILD_l5_o2_v02/v02-02-01/00015671/000/rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I402005.Pe3e3h.eL.pR.n000_002.d_rec_00015671_493.slcio \
  Output.root \
  patch.txt
```

The resulting root file will be "Output.root" 


A few things to note:
- LCRelations that do not have a `"ToType"` and `"FromType"` as collection parameters set are currently not being converted. If `check_missing_cols` is run without the `--minimal` flag a warning will be given for each Relation this applies to.
