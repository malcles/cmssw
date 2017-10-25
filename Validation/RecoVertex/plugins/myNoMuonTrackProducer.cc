// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"

//for AOD: 
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//for miniAOD:
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"


using namespace std;
using namespace edm;
using namespace reco;
using namespace math;


//
// class declaration
//

class myNoMuonTrackProducer : public edm::EDProducer 
{
public:
  explicit myNoMuonTrackProducer(const edm::ParameterSet&);
  ~myNoMuonTrackProducer();
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  EDGetTokenT<View<reco::Candidate> > recoCandidateToken_;
  EDGetTokenT<View<reco::Muon> > muonToken_;
  EDGetTokenT<View<reco::Track> > generalTrackToken_;

  // ----------member data ---------------------------
  
  //double ptMin_;

  bool doRemoveMuons = true;
  bool verbose = true;
  bool isData = false;
  bool muonEventMC = false;
  int  nMuonInAccMC=0;
  int  nZInAccMC =0;

  int nMuonReco=0;
  int nPatMuonReco=0;
  int nZReco=0;

  int muonCount = 0;
  int notMuonCount = 0;
  int notMuonBestTrackCount = 0;

};


//
// constructors and destructor
//
myNoMuonTrackProducer::myNoMuonTrackProducer(const ParameterSet& iConfig):
  recoCandidateToken_( consumes<View<reco::Candidate> >( iConfig.getParameter<InputTag> ( "recoCandidatesTag" ) ) ),
  muonToken_( consumes<View<reco::Muon> >( iConfig.getParameter<InputTag>( "muonTag" ) ) ),
  generalTrackToken_( consumes<View<reco::Track> >( iConfig.getParameter<InputTag> ( "generalTrackTag" ) ) )
    //  patMuonToken_( consumes<View<pat::Muon> >( iConfig.getParameter<InputTag>( "patMuonTag" ) ) )

{
  produces< reco::TrackCollection >();
  doRemoveMuons = iConfig.getUntrackedParameter<bool>( "doRemoveMuons", true );
  verbose = iConfig.getUntrackedParameter<bool>( "verbose", true );
  isData = iConfig.getUntrackedParameter<bool>( "isData", false );
  cout<<" doRemoveMuons="<<doRemoveMuons<<" isData="<<isData<<" verbose="<<verbose<< endl;
}


myNoMuonTrackProducer::~myNoMuonTrackProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
myNoMuonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  muonCount = 0;
  notMuonCount = 0;
  notMuonBestTrackCount = 0;

  //cout << " [DEBUG]: myNoMuonTrackProducer : Entering myNoMuonTrackProducer " << endl;

  Handle<View<reco::Candidate> > recoCandidates;
  iEvent.getByToken( recoCandidateToken_, recoCandidates );
  //  cout<<" got reco cands"<<endl;

  Handle<View<reco::Muon> >  muons;
  iEvent.getByToken( muonToken_, muons );
  //  cout<<" got muons"<<endl;

  //  Handle<View<pat::Muon> >  patMuons;
  // iEvent.getByToken( patMuonToken_, patMuons );
  //  cout<<" got pat muons"<<endl;

  Handle<View<reco::Track> > generalTracks;
  iEvent.getByToken( generalTrackToken_, generalTracks );
  cout<<" got tracks, isData= "<<isData << endl;
    
  std::unique_ptr<reco::TrackCollection> MuonLessTracks(new reco::TrackCollection);


  //cout <<" [DEBUG]: recoCandidates->size() = " << recoCandidates->size() << " muons size:"<< muons->size()<<" generalTracks size:" << generalTracks->size()<< endl;

  muonEventMC = false;
  nMuonInAccMC=0;
  
  
  
  
  nMuonReco=0;
  for(unsigned int j = 0 ; j < muons->size(); j++){
    Ptr<reco::Muon> the_muon = muons->ptrAt(j);
    //cout << " Entering Muons loop " << endl;
    //cout<<" muon #"<<j<<" pt="<<the_muon->pt()<<" eta="<<the_muon->eta()<< endl;
    if( the_muon->isPFMuon() && (the_muon->isGlobalMuon() || the_muon->isTrackerMuon()) && the_muon->pt()>10 ) { 
      if(verbose) cout<<"JM SELmuon #"<<j<<" pt="<<the_muon->pt()<<" eta="<<the_muon->eta()<< endl;
      nMuonReco++;
    }
  }

  Ptr<reco::Muon> pat_muon1;
  Ptr<reco::Muon> pat_muon2;
  double mass=0;

 if(verbose) cout << "JM nMuonsReco="<< nMuonReco<<endl;

  if( nMuonReco>=2 ){

    float smallerDM=999.;
    float Zmass=91.19;
    Ptr<reco::Muon> pat_muontmp1;
    Ptr<reco::Muon> pat_muontmp2;
  
    for(unsigned int j = 0 ; j < muons->size(); j++){
	
      pat_muontmp1= muons->ptrAt(j);
      if(! ( pat_muontmp1->isPFMuon() && ( pat_muontmp1->isGlobalMuon() ||  pat_muontmp1->isTrackerMuon()) &&  pat_muontmp1->pt()>10 ) ) continue;

      for(unsigned int i = j+1 ; i < muons->size(); i++){

	pat_muontmp2= muons->ptrAt(i);
	if(! ( pat_muontmp2->isPFMuon() && ( pat_muontmp2->isGlobalMuon() ||  pat_muontmp2->isTrackerMuon()) &&  pat_muontmp2->pt()>10 ) ) continue;
	
	XYZTLorentzVector p4Sum;
	p4Sum += pat_muontmp1->p4();
	p4Sum += pat_muontmp2->p4();

	double masstmp = p4Sum.M();
	double DM=fabs(masstmp-Zmass);
	if(DM<smallerDM){
	  pat_muon1=pat_muontmp1;
	  pat_muon2=pat_muontmp2;
	  mass=masstmp;
	  smallerDM=DM;
	}
      }
    }    
    if(verbose) cout<<"JM dimuon invariant mass:"<<mass<< endl;
  }

  bool pass=false;
  //if( mass>50 && mass<130 && nMuonReco>=2 ){
  if( mass>50 && mass<1000 && nMuonReco>=2 ){ // Changed for Pasquale
      if(verbose) cout<< "JM DIMUON CANDIDATE" << endl;
      pass=true;
    }
 
    if(pass){
      if(verbose)  cout<< "JM Muon 1:"<<pat_muon1->pt()<<" "<<pat_muon1->eta()<<" "<<pat_muon1->phi()<<endl;
      if(verbose)  cout<< "JM Muon 2:"<<pat_muon2->pt()<<" "<<pat_muon2->eta()<<" "<<pat_muon2->phi()<<endl;
    }
    
    int nRemoved=0;

    for(unsigned int j = 0 ; j < generalTracks->size(); j++){
      
      Ptr<reco::Track> the_track = generalTracks->ptrAt(j);
    
      bool doMatch=false;
      if(pass){

	if(the_track->pt()>10){
	  if(verbose) cout<<"JM theTrack#"<<j<<" pt="<<the_track->pt()<<" eta="<<the_track->eta()<<" "<<the_track->phi()<< endl; 
	  
	  float dPt1=fabs(the_track->pt()-pat_muon1->pt());
	  float dPt2=fabs(the_track->pt()-pat_muon2->pt());
	  float dR1=deltaR( the_track->eta(),  the_track->phi(),  pat_muon1->eta(), pat_muon1->phi() );
	  float dR2=deltaR( the_track->eta(),  the_track->phi(),  pat_muon2->eta(), pat_muon2->phi() );
	  
	  if(verbose) cout<<"JM dR1="<<dR1<<" dR2="<<dR2<<" dPt1="<<dPt1<<" dPt2="<<dPt2<< endl;
	  
	  if(doRemoveMuons && ((dPt1<0.01 && dR1<0.005) || (dPt2<0.01 && dR2<0.005))){
	    doMatch=true;
	    nRemoved++;
	  }	
	}
      }
      if(!doMatch) MuonLessTracks->push_back(*the_track);            
    }
    
    if(verbose && pass) cout<<" JM nRemoved="<<nRemoved<< endl;
    
    if(verbose) cout <<" JM myNoMuonTrackProducer : MuonLessTracks->size() = " << MuonLessTracks->size() <<" generalTracks->size() = "<<generalTracks->size()<< endl;
    
    iEvent.put(std::move(MuonLessTracks));  
}

// ------------ method called once each job just before starting event loop  ------------
void 
myNoMuonTrackProducer::beginJob(const edm::EventSetup&)
{
  edm::LogWarning("myNoMuonTrackProducer") << "begin job";
  LogDebug("myNoMuonTrackProducer") << "begin job";
}

// ------------ method called once each job just after ending the event loop  ------------
void 
myNoMuonTrackProducer::endJob() {

  edm::LogWarning("myNoMuonTrackProducer") << "end job";
  LogDebug("myNoMuonTrackProducer") << "end job";
}

//define this as a plug-in
DEFINE_FWK_MODULE(myNoMuonTrackProducer);
