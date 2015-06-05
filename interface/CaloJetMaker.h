#ifndef CaloJetMaker_H
#define CaloJetMaker_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class CaloJetMaker : public edm::EDProducer {
public:
    explicit CaloJetMaker (const edm::ParameterSet&);
    ~CaloJetMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::InputTag caloJetsInputTag;
    
};


#endif
