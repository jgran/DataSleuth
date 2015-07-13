import FWCore.ParameterSet.Config as cms
from Configuration.EventContent.EventContent_cff import *
import os

process = cms.Process("DATASLEUTH")

process.load("DataSleuth.DataSleuth.eventMaker_cfi")
process.load("DataSleuth.DataSleuth.metMaker_cfi")
process.load("DataSleuth.DataSleuth.caloJetMaker_cfi")
process.load("DataSleuth.DataSleuth.pfJetMaker_cfi")
process.load("DataSleuth.DataSleuth.caloTowerMaker_cfi")
process.load("DataSleuth.DataSleuth.hcalNoiseSummaryMaker_cfi")
process.load("DataSleuth.DataSleuth.hltMaker_cfi")
process.load("DataSleuth.DataSleuth.pfClusterMaker_cfi")
process.load("DataSleuth.DataSleuth.pfCandidateMaker_cfi")
process.load("DataSleuth.DataSleuth.pfmetMaker_cfi")
# process.load("DataSleuth.DataSleuth.metFilterMaker_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")

# various met filters
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.CSCTightHaloFilter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
process.load('RecoMET.METFilters.eeNoiseFilter_cfi')
process.eeNoiseFilter.taggingMode = cms.bool(True)
# process.load('RecoMET.METFilters.jetIDFailureFilter_cfi')
# process.jetIDFailure.taggingMode = cms.bool(True)

# https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/RecoMET/METFilters/python/metFilters_cff.py
## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)

process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.taggingMode = cms.bool(True)


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

process.hltParticleFlowBlock = cms.EDProducer("PFBlockProducer",
  debug = cms.untracked.bool(False),
  verbose = cms.untracked.bool(False),
  elementImporters = cms.VPSet(
      cms.PSet(
          source = cms.InputTag("particleFlowClusterECAL"),
          #source = cms.InputTag("particleFlowClusterECALUncorrected"), #we use uncorrected
          importerName = cms.string('GenericClusterImporter')
      ),
      cms.PSet(
          source = cms.InputTag("particleFlowClusterHCAL"),
          importerName = cms.string('GenericClusterImporter')
      ),
      cms.PSet(
          source = cms.InputTag("particleFlowClusterHO"),
          importerName = cms.string('GenericClusterImporter')
      ),
      cms.PSet(
          source = cms.InputTag("particleFlowClusterHF"),
          importerName = cms.string('GenericClusterImporter')
      )
  ),
  linkDefinitions = cms.VPSet(
      cms.PSet(
          linkType = cms.string('ECAL:HCAL'),
          useKDTree = cms.bool(False),
          #linkerName = cms.string('ECALAndHCALLinker')
          linkerName = cms.string('ECALAndHCALCaloJetLinker') #new ECal and HCal Linker for PFCaloJets
      ),
      cms.PSet(
          linkType = cms.string('HCAL:HO'),
          useKDTree = cms.bool(False),
          linkerName = cms.string('HCALAndHOLinker')
      ),
      cms.PSet(
          linkType = cms.string('HFEM:HFHAD'),
          useKDTree = cms.bool(False),
          linkerName = cms.string('HFEMAndHFHADLinker')
      ),
      cms.PSet(
          linkType = cms.string('ECAL:ECAL'),
          useKDTree = cms.bool(False),
          linkerName = cms.string('ECALAndECALLinker')
      )
   )
)

from RecoParticleFlow.PFProducer.particleFlow_cfi import particleFlowTmp
process.hltParticleFlow = particleFlowTmp.clone(
    GedPhotonValueMap = cms.InputTag(""),
    useEGammaFilters = cms.bool(False),
    useEGammaElectrons = cms.bool(False), 
    useEGammaSupercluster = cms.bool(False),
    rejectTracks_Step45 = cms.bool(False),
    usePFNuclearInteractions = cms.bool(False),  
    blocks = cms.InputTag("hltParticleFlowBlock"), 
    egammaElectrons = cms.InputTag(""),
    useVerticesForNeutral = cms.bool(False),
    PFEGammaCandidates = cms.InputTag(""),
    useProtectionsForJetMET = cms.bool(False),
    usePFConversions = cms.bool(False),
    rejectTracks_Bad = cms.bool(False),
    muons = cms.InputTag(""),
    postMuonCleaning = cms.bool(False),
    usePFSCEleCalib = cms.bool(False)
    )

process.load("RecoMET.METProducers.PFMET_cfi")
process.pfCaloMet = process.pfMet.clone(
  src = cms.InputTag("hltParticleFlow"),
  alias = cms.string('pfCaloMet')
)

process.pfCaloMetSequence = cms.Sequence(
    process.hltParticleFlowBlock *
    process.hltParticleFlow *
    process.pfCaloMet
)

process.pfCaloMetMaker = process.pfmetMaker.clone(
    aliasPrefix = cms.untracked.string("pfCaloMet"),
    pfMetInputTag_ = cms.InputTag("pfCaloMet")
)

process.pfChMetMaker = process.pfmetMaker.clone(
    aliasPrefix = cms.untracked.string("pfChMet"),
    pfMetInputTag_ = cms.InputTag("pfChMet")
)

process.GlobalTag.globaltag = "GR_R_74_V12"

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/home/users/namin/public_html/dump/DABF7A8D-5427-E511-AA1F-02163E01280D.root')
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
    process.eeBadScFilter *
    process.EcalDeadCellTriggerPrimitiveFilter *
    process.eeNoiseFilter *
    # process.jetIDFailure *
    process.goodVertices * process.trackingFailureFilter *
    process.hcalNoiseSummaryMaker *
    process.HBHENoiseFilterResultProducer *
    process.metMaker *
    process.pfCandidateMaker *
    process.caloJetMaker *
    process.pfJetMaker *
    process.caloTowerMaker *
    process.hltMaker *
    process.pfClusterRefsForJets_step *
    process.pfClusterMet *
    process.pfCaloMetSequence *
    process.pfCaloMetMaker *
    process.pfChMetMaker *
    process.pfmetMaker *
    process.pfClusterMaker
)

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.eventMaker.isData                        = cms.bool(True)
