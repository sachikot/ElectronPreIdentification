// -*- C++ -*-
//
// Package:    ElePF/ElePFAnalyzer
// Class:      ElePFAnalyzer
// 
/**\class ElePFAnalyzer ElePFAnalyzer.cc ElePF/ElePFAnalyzer/plugins/ElePFAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sachiko Toda
//         Created:  Sun, 07 Dec 2014 04:40:22 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// added by Sachiko
#include "RecoParticleFlow/PFTracking/plugins/GoodSeedProducer.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoParticleFlow/PFTracking/interface/ElectronSeedMerger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include <TTree.h>
#include <vector>
#include <string>
#include <iostream>
#include <TMath.h>

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class ElePFAnalyzer : public edm::EDAnalyzer {
public:
  explicit ElePFAnalyzer(const edm::ParameterSet&);
  ~ElePFAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const& run, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void FillSeedInformation(const reco::PFRecTrackCollection & PfRTColl, 
  				   reco::GsfTrack, 
  				   bool,
  				   const edm::Event&,
  				   const edm::EventSetup&);

  bool isaV( int pdgId);
  bool isaJP( int pdgId);
  bool isaBhadron(int pdgId);
  bool isaDhadron(int pdgId);
  int parent(const reco::Candidate& mom);
  //move mother code for Gen-particle here
  virtual int  tracer(const reco::Candidate& genCand);


  // ----------member data ---------------------------
  ESHandle<MagneticField> magField;
  ESHandle<TrackerGeometry> pDD;
  edm::Handle<reco::BeamSpot> theBeamSpot;
  math::XYZVector B_;
  PFTrackTransformer *pfTransformer_;
  PFTrackTransformer* pfTkTransformer_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticles_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracks_;

  TTree *EventInfo;
  Int_t nEle;
  Int_t nPion;
  Int_t nSeed;
  std::vector<Float_t> gen_ele_pt;
  std::vector<Float_t> gen_ele_eta;
  std::vector<Float_t> gen_ele_phi;
  std::vector<Int_t>   gen_ele_isFromMom;

  std::vector<Float_t> gen_pion_pt;
  std::vector<Float_t> gen_pion_eta;
  std::vector<Float_t> gen_pion_phi;
  std::vector<Int_t>   gen_pion_isFromMom;

  std::vector<Bool_t> isTrackerDriven;
  std::vector<Bool_t> isEcalDriven;
  std::vector<Float_t> gsf_seed_pt;
  std::vector<Float_t> gsf_seed_eta;
  std::vector<Float_t> gsf_seed_phi;
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
ElePFAnalyzer::ElePFAnalyzer(const edm::ParameterSet& iConfig):
  pfTransformer_(0),
  pfTkTransformer_(0),
  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  gsfTracks_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("gsfTracks")))
{

  edm::Service<TFileService> fs;
  EventInfo = fs->make<TTree> ("EventTree", "Event Tree");
  EventInfo->Branch("nEle", &nEle, "nEle/I");
  EventInfo->Branch("nPion", &nPion, "nPion/I");
  EventInfo->Branch("nSeed", &nSeed, "nSeed/I");
  EventInfo->Branch("gen_ele_pt", &gen_ele_pt);
  EventInfo->Branch("gen_ele_eta", &gen_ele_eta);
  EventInfo->Branch("gen_ele_phi", &gen_ele_phi);
  EventInfo->Branch("gen_ele_isFromMom", &gen_ele_isFromMom);

  EventInfo->Branch("gen_pion_pt", &gen_pion_pt);
  EventInfo->Branch("gen_pion_eta", &gen_pion_eta);
  EventInfo->Branch("gen_pion_phi", &gen_pion_phi);
  EventInfo->Branch("gen_pion_isFromMom", &gen_pion_isFromMom);

  EventInfo->Branch("TrackerDrivenSeed", &isTrackerDriven);
  EventInfo->Branch("EcalDrivenSeed", &isEcalDriven);
  EventInfo->Branch("gsf_seed_pt", &gsf_seed_pt);
  EventInfo->Branch("gsf_seed_eta", &gsf_seed_eta);
  EventInfo->Branch("gsf_seed_phi", &gsf_seed_phi);
}

ElePFAnalyzer::~ElePFAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete pfTransformer_;
  delete pfTkTransformer_;

  pfTransformer_ = nullptr;
  pfTkTransformer_ = nullptr;
}


//
// member functions
//--------------------------------------------------
//0 - not from W, Z, B, or C
//1 - W or Z
//2 - B
//3 - C

bool ElePFAnalyzer::isaV(int pidAbs){
  int res = false;
  if(pidAbs==24 || pidAbs==23 || pidAbs==22)res=true;
  //if (pidAbs ==23) res= true; // for Zee events
  return res;
}

bool ElePFAnalyzer::isaBhadron(int pidAbs){

  bool isB = false;
  if(pidAbs>500){
              
    pidAbs/= 100;
                
    if(pidAbs<60 && pidAbs>50)pidAbs/= 10;
    int mod10 = pidAbs % 5;
               
    if(mod10 == 0) {                       
      isB = true;                      
    }
  }
  else  {
    //              std::cout<<"PID too low to be a B meson"<<std::endl;
  }
  return isB;
}

bool ElePFAnalyzer::isaDhadron(int pidAbs){

  bool isD = false;
  if(pidAbs>400){
    pidAbs/= 100;
    if(pidAbs<50 && pidAbs>40)pidAbs/= 10;
    int mod10 = pidAbs % 4;
                
    if(mod10 == 0 && pidAbs<10) {
      isD = true;                       
    }
  }
  else{
                
    //std::cout<<"PID too low to be a B meson"<<std::endl;
  }
  return isD;
}

int ElePFAnalyzer::parent(const reco::Candidate& mom){
  int pdgId = mom.pdgId();

  if ( isaV      (fabs(pdgId)) ) return 1;
  if ( isaBhadron(fabs(pdgId)) ) return 2;
  if ( isaDhadron(fabs(pdgId)) ) return 3;
      
  return 0;
}

int ElePFAnalyzer::tracer(const reco::Candidate& genCand) {

  int parentId = 0;
  if ( genCand.numberOfMothers() != 0 ){

    const reco::Candidate& mom0 = *genCand.mother(0);
    parentId = parent(mom0);
    
    if ( parentId == 0 && !(mom0.pdgId() == 111 || mom0.pdgId()==221) ){
      if (mom0.numberOfMothers() != 0 ){
	const reco::Candidate& mom1 = *mom0.mother(0);
	parentId = parent(mom1);
      
	if(parentId == 0){
	  if (mom1.numberOfMothers() != 0 ){
	    const reco::Candidate& mom2 = *mom1.mother(0);
	    parentId = parent(mom2);
	  
	    if (parentId == 0 ){
	      if (mom2.numberOfMothers() != 0 ){
		const reco::Candidate& mom3 = *mom2.mother(0);
		parentId = parent(mom3);
	      
		if(parentId == 0){
		  if (mom3.numberOfMothers() != 0 ){
		    const reco::Candidate& mom4 = *mom3.mother(0);
		    parentId = parent(mom4);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return parentId;
}

void ElePFAnalyzer::FillSeedInformation(const reco::PFRecTrackCollection & PfRTkColl,
					reco::GsfTrack gsftracks,
					bool otherColl,
					const edm::Event& iEvent,
					const edm::EventSetup& iSetup){

  const reco::BeamSpot& bs = *(theBeamSpot.product());
  TSCBLBuilderNoMaterial tscblBuilder;
  const GeomDet *detc = 0;

  TrajectorySeed::const_iterator hit_end = gsftracks.seedRef()->recHits().second;
  //You don't need to loop over for hits. The good direction is given by the latest hit
  detc = pDD->idToDet( (*(--hit_end)).geographicalId() );

  TrajectoryStateOnSurface t = trajectoryStateTransform::transientState((*gsftracks.seedRef()).startingState(), &(detc->surface()), &(*magField));
  const FreeTrajectoryState & stateForProjectionToBeamLine = *t.freeState();
  TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine, bs);
  GlobalPoint v = tscbl.trackStateAtPCA().position();
  math::XYZPoint pos( v.x(), v.y(), v.z() );
  GlobalVector p = tscbl.trackStateAtPCA().momentum();
  math::XYZVector mom( p.x(), p.y(), p.z() );

  float seedPt = sqrt(p.x()*p.x() + p.y()*p.y());
  gsf_seed_pt.push_back(seedPt);
  gsf_seed_eta.push_back(mom.eta());
  gsf_seed_phi.push_back(mom.phi());
}// end of fill seed info

// ------------ method called for each event  ------------
void
ElePFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  gen_ele_pt.clear();
  gen_ele_eta.clear();
  gen_ele_phi.clear();
  gen_ele_isFromMom.clear();

  gen_pion_pt.clear();
  gen_pion_eta.clear();
  gen_pion_phi.clear();
  gen_pion_isFromMom.clear();
  isTrackerDriven.clear();
  isEcalDriven.clear();
  gsf_seed_pt.clear();
  gsf_seed_eta.clear();
  gsf_seed_phi.clear();

  // Handles, collections etc. ===========================================
  Handle<GenParticleCollection> genParticle;
  iEvent.getByToken(genParticles_, genParticle);

  Handle<reco::GsfTrackCollection> gsfTrackColl;
  iEvent.getByToken(gsfTracks_, gsfTrackColl);
  GsfTrackCollection gsftracks = *(gsfTrackColl.product());

  Handle<reco::PFRecTrackCollection> thePfRecTrackCollection;
  iEvent.getByLabel("pfTrack",thePfRecTrackCollection);
  const PFRecTrackCollection& PfRTkColl = *(thePfRecTrackCollection.product());

  nEle = 0;
  nPion = 0;
  for( size_t i = 0; i < genParticle->size(); ++i ){
    const GenParticle& genCand = (*genParticle)[i];

    bool electron = abs(genCand.pdgId()) == 11;
    bool pion     = abs(genCand.pdgId()) == 211;
   
    int parentId = tracer(genCand);
   
    if (electron == 1){
      gen_ele_pt.push_back(genCand.pt());
      gen_ele_phi.push_back(genCand.phi());
      gen_ele_eta.push_back(genCand.eta());
      gen_ele_isFromMom.push_back(parentId);
      ++nEle;
    } 
    if(pion == 1) {
      gen_pion_pt.push_back(genCand.pt());
      gen_pion_phi.push_back(genCand.phi());
      gen_pion_eta.push_back(genCand.eta());
      gen_pion_isFromMom.push_back(parentId);
      ++nPion;
    }
  } // loop over for gen paticle

  iSetup.get<IdealMagneticFieldRecord>().get(magField);
  iSetup.get<TrackerDigiGeometryRecord>().get(pDD);

  InputTag beamSpot_(string("offlineBeamSpot"));
  iEvent.getByLabel(beamSpot_, theBeamSpot);

  nSeed = 0;
  for (unsigned int igsf=0; igsf<gsftracks.size();igsf++) {
    GsfTrackRef trackRef(gsfTrackColl, igsf);

    //reco-seed information
    FillSeedInformation(PfRTkColl, gsftracks[igsf], false, iEvent, iSetup);

    bool isTrackerDrivenForThisSeed = true;
    bool isEcalDrivenForThisSeed = true;
    ElectronSeedRef SeedRef= trackRef->extra()->seedRef().castTo<ElectronSeedRef>();

    if (SeedRef->caloCluster().isNull()){
      isEcalDrivenForThisSeed = false;
    }
    if (SeedRef->ctfTrack().isNull()){
      isTrackerDrivenForThisSeed = false;
    }

    isTrackerDriven.push_back(isTrackerDrivenForThisSeed);
    isEcalDriven.push_back(isEcalDrivenForThisSeed);
    ++nSeed;
  }

  EventInfo->Fill();
}


  // ------------ method called once each job just before starting event loop  ------------
  void 
    ElePFAnalyzer::beginJob()
  {
  }

  // ------------ method called once each job just after ending the event loop  ------------
  void 
    ElePFAnalyzer::endJob() 
  {
  }

  // ------------ method called when starting to processes a run  ------------

void 
ElePFAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& es)
{
  //Magnetic Field
  ESHandle<MagneticField> magneticField;
  es.get<IdealMagneticFieldRecord>().get(magneticField);
  B_=magneticField->inTesla(GlobalPoint(0,0,0));
  
  pfTransformer_= new PFTrackTransformer(B_);
  pfTransformer_->OnlyProp();
  
  pfTkTransformer_= new PFTrackTransformer(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));
  pfTkTransformer_->OnlyProp();  
}


  // ------------ method called when ending the processing of a run  ------------

void 
ElePFAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
  delete pfTransformer_;
  pfTransformer_ = nullptr;
  delete pfTkTransformer_;
  pfTkTransformer_ = nullptr;
}

  // ------------ method called when starting to processes a luminosity block  ------------
  /*
    void 
    ElePFAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }
  */

  // ------------ method called when ending the processing of a luminosity block  ------------
  /*
    void 
    ElePFAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
    {
    }
  */

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void
    ElePFAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(ElePFAnalyzer);
