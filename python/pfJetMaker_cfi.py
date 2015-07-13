import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets"),
  pfJetsInputTag                   = cms.InputTag("ak4PFJets","","RECO"),
  pfCandidatesTag                  = cms.InputTag("particleFlow"),
  pfJetPtCut                       = cms.double(5.),
)
