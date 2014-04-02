// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeProducer
// 
/**\class CascadeProducer CascadeProducer.cc MyProducers/CascadeProducer/src/CascadeProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: CascadeProducer.cc,v 1.12 2010/02/20 21:02:02 wmtan Exp $
//
//


// system include files
#include <memory>

#include "RecoVertex/V0Producer/interface/CascadeProducer.h"

// Constructor
CascadeProducer::CascadeProducer(const edm::ParameterSet& iConfig) :
  theParams(iConfig) {

  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >("Cascade");
}

// (Empty) Destructor
CascadeProducer::~CascadeProducer() {

}


//
// Methods
//

// Producer Method
void CascadeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   CascadeFitter theVees(theParams, iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     casCandidates( new reco::VertexCompositeCandidateCollection );
   casCandidates->reserve( theVees.getCass().size() );

   std::copy( theVees.getCass().begin(),
              theVees.getCass().end(),
              std::back_inserter(*casCandidates) );

   // Write the collections to the Event
   iEvent.put( casCandidates, std::string("Cascade") );
}


void CascadeProducer::beginJob() {
}


void CascadeProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(CascadeProducer);
