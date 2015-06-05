import FWCore.ParameterSet.Config as cms

eventMaker = cms.EDProducer(
  "EventMaker",
  aliasPrefix = cms.untracked.string("evt"),
  datasetName = cms.string("undefined"),
  dcsTag      = cms.InputTag("scalersRawToDigi"),
  isData      = cms.bool(True)
)
