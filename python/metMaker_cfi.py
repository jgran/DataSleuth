import FWCore.ParameterSet.Config as cms

metMaker = cms.EDProducer("METMaker",
	aliasPrefix = cms.untracked.string("evt"),
                        met_tag_               = cms.InputTag("caloMet"              ),               
                        metHO_tag_             = cms.InputTag("caloMetBEFO"          ),             
                        metNoHF_tag_           = cms.InputTag("caloMetBE"            ),           
                        metNoHFHO_tag_         = cms.InputTag("metNoHFHO"            ),                        
                        metOpt_tag_            = cms.InputTag(""                     ),            
                        metOptHO_tag_          = cms.InputTag(""                     ),          
                        metOptNoHF_tag_        = cms.InputTag(""                     ),        
                        metOptNoHFHO_tag_      = cms.InputTag(""                     ),     
                        corMetGlobalMuons_tag_ = cms.InputTag("caloMetM"             ),
                        MuonJEScorMET_tag_     = cms.InputTag("caloType1CorrectedMet"),
                        muon_tag_              = cms.InputTag("muons"                ),
                        muon_vm_tag_           = cms.InputTag("muonMETValueMapProducer", "muCorrData"),
                        hbheNoiseFilterInputTag_ = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"),
                        towerEtThreshold_      = cms.double(0.3),
                        make_eta_rings_        = cms.bool(True)
)                                                              
