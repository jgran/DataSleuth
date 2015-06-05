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

typedef math::XYZTLorentzVector LorentzVector_;


CaloJetMaker::CaloJetMaker(const edm::ParameterSet& iConfig) {

    produces<std::vector<float> > ("calojetspt").setBranchAlias("calojets_pt");
    produces<std::vector<float> > ("calojetseta").setBranchAlias("calojets_eta");
    produces<std::vector<float> > ("calojetsphi").setBranchAlias("calojets_phi");

    caloJetsInputTag = iConfig.getParameter<edm::InputTag>("caloJetsInputTag_");
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

    edm::Handle<edm::View<reco::CaloJet> > calojet_h;
    iEvent.getByLabel(caloJetsInputTag, calojet_h);

    for(edm::View<reco::CaloJet>::const_iterator calojet_it = calojet_h->begin(); calojet_it != calojet_h->end(); calojet_it++){

      if(calojet_it->pt() < 5.0) continue;

      calojets_pt  ->push_back(calojet_it->pt());
      calojets_eta ->push_back(calojet_it->eta());
      calojets_phi ->push_back(calojet_it->phi());

    } 

    iEvent.put(calojets_pt,   "calojetspt" );
    iEvent.put(calojets_eta,  "calojetseta" );
    iEvent.put(calojets_phi,  "calojetsphi" );
}

//define this as a plug-in
DEFINE_FWK_MODULE(CaloJetMaker);
