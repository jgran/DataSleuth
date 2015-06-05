import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *

process = cms.Process("DATASLEUTH")

process.load("DataSleuth.DataSleuth.eventMaker_cfi")
process.load("DataSleuth.DataSleuth.metMaker_cfi")
process.load("DataSleuth.DataSleuth.caloJetMaker_cfi")
process.load("DataSleuth.DataSleuth.caloTowerMaker_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = "GR_R_74_V12"

#process.source = cms.Source("PoolSource",
#       fileNames = cms.untracked.vstring('file:/nfs-7/userdata/jgran/Run2015A/ExpressPhysics/FEVT/Express-v1/0EA17D6D-B609-E511-9404-02163E014682.root')
#)

process.out = cms.OutputModule("PoolOutputModule",
  fileName     = cms.untracked.string('ntuple.root'),
  dropMetaData = cms.untracked.string("NONE")
)
process.outpath = cms.EndPath(process.out)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_DATASLEUTH*'))

process.p = cms.Path( 
  process.eventMaker *
  process.metMaker *
  process.caloJetMaker *
  process.caloTowerMaker
)

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.eventMaker.isData                        = cms.bool(True)
