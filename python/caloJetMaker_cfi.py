import FWCore.ParameterSet.Config as cms

caloJetMaker = cms.EDProducer(
  "CaloJetMaker",
  caloJetsInputTag_ = cms.InputTag("ak4CaloJets")
)
