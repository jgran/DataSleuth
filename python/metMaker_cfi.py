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
                        # hbheIsoNoiseFilterInputTag_ = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),
                        hbheNoiseFilterRun1InputTag_ = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun1"),
                        hbheNoiseFilterRun2LooseInputTag_ = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
                        hbheNoiseFilterRun2TightInputTag_ = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
                        cscTightHaloFilterInputTag_ = cms.InputTag("CSCTightHaloFilter", ""),
                        EcalDeadCellTriggerPrimitiveFilterInputTag_ = cms.InputTag("EcalDeadCellTriggerPrimitiveFilter", ""),
                        EcalDeadCellBoundaryEnergyFilterInputTag_ = cms.InputTag("EcalDeadCellBoundaryEnergyFilter", ""),
                        EcalLaserCorrFilterInputTag_ = cms.InputTag("EcalLaserCorrFilter", ""),
                        eeBadScFilterInputTag_ = cms.InputTag("eeBadScFilter", ""),
                        goodVerticesInputTag_ = cms.InputTag("goodVertices", ""),
                        trackingFailureFilterInputTag_ = cms.InputTag("trackingFailureFilter", ""),
                        eeNoiseFilterInputTag_ = cms.InputTag("eeNoiseFilter", ""),

                        towerEtThreshold_      = cms.double(0.3),
                        make_eta_rings_        = cms.bool(True)
)                                                              
