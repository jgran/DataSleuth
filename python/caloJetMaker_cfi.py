import FWCore.ParameterSet.Config as cms

caloJetMaker = cms.EDProducer(
  "CaloJetMaker",
  caloJetsInputTag_ = cms.InputTag("ak4CaloJets"),
  # caloJetsCorrectorL1FastL2L3 = cms.string('ak4CaloL1FastL2L3')
)
