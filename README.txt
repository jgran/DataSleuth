setup:

export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_4_patch1
cd CMSSW_7_4_4_patch1/src
cmsenv
mkdir DataSleuth
cd DataSleuth
git clone git@github.com:jgran/DataSleuth.git
scram b -j5

to run:
cmsRun test/run_cfg.py
