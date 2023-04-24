To convert an lcio file to edm4HEP either a lcio file where all collections are present in every event
or a file containing a list of all collections in the file are needed. 
This list can be obtained by executing check_missing_cols from the lcio package with the --minimal flag.

Example:
./lcio/install/bin/check_missing_cols /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/higgs/ILD_l5_o2_v02/v02-02-01/00015671/000/rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I402005.Pe3e3h.eL.pR.n000_002.d_rec_00015671_493.slcio > patch.txt

This file will be the third argument when the converter ist called. 

./standalone/lcio2edm4hep /pnfs/desy.de/ilc/prod/ilc/mc-2020/ild/rec/250-SetA/higgs/ILD_l5_o2_v02/v02-02-01/00015671/000/rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I402005.Pe3e3h.eL.pR.n000_002.d_rec_00015671_493.slcio Output.root patch.txt

The resulting root file will be "Output.root" 


A few things to note:
ralations that do not have a "toType" and "fromType" set are currently not beeing converted. If check_missing_cols is run without the minimal flag a warning will be given for each Ralations this applies to.
