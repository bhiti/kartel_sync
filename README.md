# kartel_sync
Test beam synchronization script

Ljubljana beam telescope data conversion and synchronization
Based on aligning the time stamps in two files. Missing events in the telescope file are filled with dummy events (0 tracks).

------------------------------------------------------------
REQUIREMENTS:
Judith is required for converting Kartel raw files to root:
https://svnweb.cern.ch/cern/wsvn/atlasdbm/kartel/Judith/branches/JudithKartelWf/

A Judith Twiki:
https://twiki.cern.ch/twiki/bin/view/Main/JudithSw

PyTables2RootTTree converter (branch judith) is required to convert h5 files (FEI4) to root:
https://github.com/SiLab-Bonn/PyTables2RootTTree/tree/judith

------------------------------------------------------------
WORKFLOW
Telescope data has to be converted to root:
./Judith -c convert -i $RAW.bin -o $CONV.root -n 3000000
input name:   acq<RUN>.bin
output name:  ref<RUN>.root

FEI4 data has to be converted to root:
a) using PyTables2RootTTree
b) In folder h5conv use:
python convert_simple.py <FEI4_RUN_NUMBER> <TELESCOPE_RUN_NUMBER>
input name:   <FEI4_RUN_NUMBER>_proto9_ext_trigger_scan.h5
output name:  anchor<RUN>.root

Telescope and FEI4 have to be synchronized:
In kartel_sync use:
root -l
.L synchronization.C+
synchronize(fname_tel, fname_anchor, fname_synced, Nevents)
