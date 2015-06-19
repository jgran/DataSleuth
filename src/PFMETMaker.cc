//-*- C++ -*-
//
// Package:    PFMETMaker
// Class:      PFMETMaker
// 
/**\class PFMETMaker PFMETMaker.cc CMS2/PFMETMaker/src/PFMETMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: PFMETMaker.cc,v 1.11 2012/05/09 23:41:32 fgolf Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataSleuth/DataSleuth/interface/PFMETMaker.h"

#include "DataFormats/METReco/interface/PFMET.h"

typedef math::XYZTLorentzVectorF LorentzVector;

//
// constructors and destructor
//

PFMETMaker::PFMETMaker(const edm::ParameterSet& iConfig) {

    aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
    std::string branchprefix = aliasprefix_;
     if(branchprefix.find("_") != std::string::npos)
	  branchprefix.replace(branchprefix.find("_"),1,"");

    produces<float> (branchprefix+"met"                 ).setBranchAlias(aliasprefix_+"_met"          );
    produces<float> (branchprefix+"metPhi"              ).setBranchAlias(aliasprefix_+"_metPhi"       );
    produces<float> (branchprefix+"metSig"              ).setBranchAlias(aliasprefix_+"_metSig"       ); //this is just MET/sqrt(sumET). Use evt_metSignificance unless you really want this branch
    produces<float> (branchprefix+"sumet"               ).setBranchAlias(aliasprefix_+"_sumet"        );
    produces<float> (branchprefix+"metSignificance"     ).setBranchAlias(aliasprefix_+"_metSignificance");
    produces<float> (branchprefix+"mettype1cor"         ).setBranchAlias(aliasprefix_+"_met_type1cor");
    produces<float> (branchprefix+"metPhitype1cor"      ).setBranchAlias(aliasprefix_+"_metPhi_type1cor");

    pfMetInputTag = iConfig.getParameter<edm::InputTag>("pfMetInputTag_");
    pfMetCorInputTag = iConfig.getParameter<edm::InputTag>("pfMetCorInputTag_");
}


PFMETMaker::~PFMETMaker() {}

void  PFMETMaker::beginJob() {
}

void PFMETMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void PFMETMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::auto_ptr<float>   evt_met         (new float   );
    std::auto_ptr<float>   evt_metPhi      (new float   );
    std::auto_ptr<float>   evt_metSig      (new float   ); //this is just MET/sqrt(sumET). Use evt_metSignificance unless you really want this branch
    std::auto_ptr<float>   evt_sumet       (new float   );
    std::auto_ptr<float>   evt_metSignificance(new float   );
    std::auto_ptr<float>   evt_met_type1cor         (new float   );
    std::auto_ptr<float>   evt_metPhi_type1cor      (new float   );

    edm::Handle<edm::View<reco::PFMET> > met_h;
    iEvent.getByLabel(pfMetInputTag, met_h);

    edm::Handle<edm::View<reco::PFMET> > metcor_h;
    iEvent.getByLabel(pfMetCorInputTag, metcor_h);

    if( !met_h.isValid() ) {
        edm::LogInfo("OutputInfo") << " failed to retrieve particle-flow MET collection";
        edm::LogInfo("OutputInfo") << " PFMETMaker cannot continue...!";
        return;
    }

    *evt_met    = ( met_h->front() ).et();
    *evt_metPhi = ( met_h->front() ).phi();
    *evt_metSig = ( met_h->front() ).mEtSig();
    *evt_sumet  = ( met_h->front() ).sumEt();       

    try { 
        *evt_metSignificance = ( met_h->front() ).significance();
    }
    catch ( cms::Exception& ex ) {
        *evt_metSignificance = -9999;
    }

    std::string branchprefix = aliasprefix_;
    if(branchprefix.find("_") != std::string::npos)
        branchprefix.replace(branchprefix.find("_"),1,"");

    iEvent.put(evt_met             , branchprefix+"met"      );
    iEvent.put(evt_metPhi          , branchprefix+"metPhi"   );
    iEvent.put(evt_metSig          , branchprefix+"metSig"   );
    iEvent.put(evt_sumet           , branchprefix+"sumet"    );  
    iEvent.put(evt_metSignificance , branchprefix+"metSignificance" );  

    if( !metcor_h.isValid() ) {
        edm::LogInfo("OutputInfo") << " failed to corrected retrieve particle-flow MET collection";
        edm::LogInfo("OutputInfo") << " PFMETMaker cannot continue...!";
        return;
    }

    *evt_met_type1cor = ( metcor_h->front() ).et();
    *evt_metPhi_type1cor = ( metcor_h->front() ).phi();

    iEvent.put(evt_met_type1cor    , branchprefix+"mettype1cor"      );
    iEvent.put(evt_metPhi_type1cor , branchprefix+"metPhitype1cor"   );
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFMETMaker);
