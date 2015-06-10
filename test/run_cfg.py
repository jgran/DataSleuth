import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *

process = cms.Process("DATASLEUTH")

process.load("DataSleuth.DataSleuth.eventMaker_cfi")
process.load("DataSleuth.DataSleuth.metMaker_cfi")
process.load("DataSleuth.DataSleuth.caloJetMaker_cfi")
process.load("DataSleuth.DataSleuth.caloTowerMaker_cfi")
process.load("DataSleuth.DataSleuth.hcalNoiseSummaryMaker_cfi")
process.load("DataSleuth.DataSleuth.hltMaker_cfi")
process.load("DataSleuth.DataSleuth.pfClusterMaker_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")

process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.CSCTightHaloFilter.taggingMode = cms.bool(True)

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('RecoMET.METProducers.PFClusterMET_cfi')
process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")
process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.pfClusterRefsForJetsHCAL = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterHCAL'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJetsECAL = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterECAL'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJetsHF = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterHF'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJetsHO = cms.EDProducer("PFClusterRefCandidateProducer",
  src          = cms.InputTag('particleFlowClusterHO'),
  particleType = cms.string('pi+')
)
process.pfClusterRefsForJets = cms.EDProducer("PFClusterRefCandidateMerger",
  src = cms.VInputTag("pfClusterRefsForJetsHCAL", "pfClusterRefsForJetsECAL", "pfClusterRefsForJetsHF", "pfClusterRefsForJetsHO"))

process.pfClusterRefsForJets_step = cms.Sequence(
    process.particleFlowRecHitECAL*
    process.particleFlowRecHitHBHE*
    process.particleFlowRecHitHF*
    process.particleFlowRecHitHO*
    process.particleFlowClusterECALUncorrected*
    process.particleFlowClusterECAL*
    process.particleFlowClusterHBHE*
    process.particleFlowClusterHCAL*
    process.particleFlowClusterHF*
    process.particleFlowClusterHO*
    process.pfClusterRefsForJetsHCAL*
    process.pfClusterRefsForJetsECAL*
    process.pfClusterRefsForJetsHF*
    process.pfClusterRefsForJetsHO*
    process.pfClusterRefsForJets
)

process.GlobalTag.globaltag = "GR_R_74_V12"

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/nfs-7/userdata/jgran/Run2015A/ExpressPhysics/FEVT/Express-v1/0EA17D6D-B609-E511-9404-02163E014682.root')
                            # fileNames = cms.untracked.vstring('file:/home/users/fgolf/run2/met/stripped_events/HighMET_newHcalNoiseFilt_246074-246214.root')
)

process.out = cms.OutputModule("PoolOutputModule",
                               # fileName     = cms.untracked.string('HighMET_newHcalNoiseFilt_246074-246214_ntuple.root'),
                               fileName     = cms.untracked.string('ntuple.root'),                               
                               dropMetaData = cms.untracked.string("NONE")
)
process.outpath = cms.EndPath(process.out)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_DATASLEUTH*'))
# process.out.outputCommands.extend(cms.untracked.vstring('keep *_*_*_DATASLEUTH*'))                                              
# process.out.outputCommands = cms.untracked.vstring( 'keep *' )

process.p = cms.Path( 
    process.eventMaker *
    process.CSCTightHaloFilter *
    process.hcalNoiseSummaryMaker *
    process.HBHENoiseFilterResultProducer *
    process.metMaker *
    process.caloJetMaker *
    process.caloTowerMaker *
    process.hltMaker *
    process.pfClusterRefsForJets_step *
    process.pfClusterMet *
    process.pfClusterMaker
)

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.eventMaker.isData                        = cms.bool(True)
