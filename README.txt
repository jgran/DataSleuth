setup:

export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_4_patch1
cd CMSSW_7_4_4_patch1/src
cmsenv
git clone git@github.com:jgran/DataSleuth.git DataSleuth/DataSleuth
pushd DataSleuth/DataSleuth
source setup/setup.sh
cd $CMSSW_BASE/src
scram b -j5

to run:
cmsRun test/run_cfg.py
