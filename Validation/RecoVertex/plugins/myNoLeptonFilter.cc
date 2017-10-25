// -*- C++ -*-
//
// Package:    Validation/RecoVertex
// Class:      myNoLeptonFilter
// 
/**\class myNoLeptonFilter myNoLeptonFilter.cc Validation/RecoVertex/plugins/myNoLeptonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Julie Malcles
//         Created:  Thu, 29 Oct 2015 09:40:52 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"

//for AOD: 
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

using namespace std;
using namespace edm;
//using namespace math;

class myNoLeptonFilter : public edm::EDFilter {
   public:
      explicit myNoLeptonFilter(const edm::ParameterSet&);
      ~myNoLeptonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
  EDGetTokenT<View<reco::Track> > allTrackToken_;
  EDGetTokenT<View<reco::Track> > noLepTrackToken_;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
myNoLeptonFilter::myNoLeptonFilter(const edm::ParameterSet& iConfig):
  allTrackToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag> ( "allTrackTag" ) ) ),
  noLepTrackToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag> ( "noLepTrackTag" ) ) )
{
   //now do what ever initialization is needed

}


myNoLeptonFilter::~myNoLeptonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
myNoLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //cout<<" FILTER"<< endl;
   using namespace edm;
  Handle<View<reco::Track> > allTracks;
  iEvent.getByToken( allTrackToken_, allTracks );
  //cout<<" FILTER got tracks"<< endl;
  Handle<View<reco::Track> > noLepTracks;
  iEvent.getByToken( noLepTrackToken_, noLepTracks );
  //cout<<" FILTER got tracks no mu"<< endl;
  if(allTracks->size()-noLepTracks->size()!=2) return false;
  else cout<<" PASSING"<< endl;
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
myNoLeptonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
myNoLeptonFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
myNoLeptonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
myNoLeptonFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
myNoLeptonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
myNoLeptonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
myNoLeptonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(myNoLeptonFilter);
