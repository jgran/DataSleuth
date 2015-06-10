import FWCore.ParameterSet.Config as cms

pfClusterMaker = cms.EDProducer("PFClusterMaker",
                                aliasPrefix = cms.untracked.string("pfcluster"),
                                PFClustersECAL = cms.InputTag("particleFlowClusterECAL"),
                                PFClustersHCAL = cms.InputTag("particleFlowClusterHCAL"),
                                PFClustersHF = cms.InputTag("particleFlowClusterHF"),
                                PFClustersHO = cms.InputTag("particleFlowClusterHO"),
                                PFClusterMET = cms.InputTag("pfClusterMet")
)                                
