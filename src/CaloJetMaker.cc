#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataSleuth/DataSleuth/interface/CaloJetMaker.h"

// #include "JetMETCorrections/Objects/interface/JetCorrector.h"

typedef math::XYZTLorentzVector LorentzVector_;


CaloJetMaker::CaloJetMaker(const edm::ParameterSet& iConfig) {

    produces<std::vector<float> > ("calojetspt").setBranchAlias("calojets_pt");
    produces<std::vector<float> > ("calojetseta").setBranchAlias("calojets_eta");
    produces<std::vector<float> > ("calojetsphi").setBranchAlias("calojets_phi");
    // produces<std::vector<float> > ("calojetscorL1FastL2L3").setBranchAlias("calojets_corL1FastL2L3"  ); 

    caloJetsInputTag = iConfig.getParameter<edm::InputTag>("caloJetsInputTag_");
    // caloJetsCorrectorL1FastL2L3 = iConfig.getParameter<std::string>("caloJetsCorrectorL1FastL2L3");
}


CaloJetMaker::~CaloJetMaker() {}

void  CaloJetMaker::beginJob() {
}

void CaloJetMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void CaloJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
    std::auto_ptr<std::vector<float> > calojets_pt   (new std::vector<float>);
    std::auto_ptr<std::vector<float> > calojets_eta  (new std::vector<float>);
    std::auto_ptr<std::vector<float> > calojets_phi  (new std::vector<float>);
    // std::auto_ptr<std::vector<float> > calojets_corL1FastL2L3  (new std::vector<float>);

    edm::Handle<edm::View<reco::CaloJet> > calojet_h;
    iEvent.getByLabel(caloJetsInputTag, calojet_h);

    // const JetCorrector* correctorL1FastL2L3 = JetCorrector::getJetCorrector(caloJetsCorrectorL1FastL2L3, iSetup);
    for(edm::View<reco::CaloJet>::const_iterator calojet_it = calojet_h->begin(); calojet_it != calojet_h->end(); calojet_it++){

      if(calojet_it->pt() < 5.0) continue;

      // float corL1FastL2L3 = correctorL1FastL2L3->correction(*calojet_it,   iEvent, iSetup);

      calojets_pt  ->push_back(calojet_it->pt());
      calojets_eta ->push_back(calojet_it->eta());
      calojets_phi ->push_back(calojet_it->phi());
      // calojets_corL1FastL2L3 ->push_back(corL1FastL2L3);

    } 

    iEvent.put(calojets_pt,   "calojetspt" );
    iEvent.put(calojets_eta,  "calojetseta" );
    iEvent.put(calojets_phi,  "calojetsphi" );
    // iEvent.put(calojets_corL1FastL2L3,  "calojetscorL1FastL2L3" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloJetMaker);
