// -*- C++ -*-
//
// Package:    ElePFAnalyzer
// Class:      ElePFAnalyzer
// 
/**\class ElePFAnalyzer ElePFAnalyzer.cc ElePF/ElePFAnalyzer/plugins/ElePFAnalyzer.cc
*/
//
// Original Authors:  Sachiko Toda & Yurii Maravin
// $Id$
//
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedTrackerVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoParticleFlow/PFTracking/interface/ElectronSeedMerger.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "RecoParticleFlow/PFTracking/interface/GoodSeedProducer.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h" 
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//add for the input variables
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFResolutionMap.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"  
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFEGammaAlgo.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "RecoParticleFlow/PFProducer/interface/PFElectronExtraEqual.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TVector2.h"
#include "Math/VectorUtil.h"
#include <algorithm>

using namespace edm;
using namespace std;
using namespace reco;

typedef std::pair<const reco::PFBlockElement*,bool> PFFlaggedElement;

class ElePFAnalyzer : public edm::EDAnalyzer {
public:
  explicit ElePFAnalyzer(const edm::ParameterSet&);
  ~ElePFAnalyzer();

private:
  virtual void clearBDTTreeVariables();
  // returns parentId
  //move mother code for Gen-particle here
  virtual int  tracer(const reco::Candidate& genCand);
  // returns a track per given genCand or false
  virtual bool findMatch(const reco::Candidate& genCand, 
			 Handle<TrackingParticleCollection> TPCollectionH,
			 reco::SimToRecoCollection q,
			 reco::Track& track,
			 int&   indexTrack,
			 float& dR,
			 float& sharedHits);
			
  // checks if two of the tracks are innermost false if track1 is near, true if track2 is near, returns true for sameLayer if they are in the same layer
  virtual bool isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer);

  // returns a colection of pairs(track, sharedHits) of a given pdgId in the event
  virtual void findRecoTracks(Handle<TrackingParticleCollection> TPCollectionH,
			      reco::SimToRecoCollection q, 
			      std::vector<std::pair<const reco::Track, float> >& tracksAndSharedHits, 
			      int pdgId);

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(const edm::Run & run, const edm::EventSetup &) override;
  virtual void endRun(const edm::Run & run, const edm::EventSetup &) override;

  virtual void FillSeedInformation(const reco::PFRecTrackCollection & PfRTColl, 
				   reco::ElectronSeed, 
				   bool, 
				   const edm::Event&, 
				   const edm::EventSetup&);


  bool isaV( int pdgId);
  bool isaBhadron(int pdgId);
  bool isaDhadron(int pdgId);
  int parent(const reco::Candidate& mom);

  void fillPreIdRefValueMap( edm::Handle<reco::TrackCollection> tkhandle,
			     const edm::OrphanHandle<reco::PreIdCollection>&,
			     edm::ValueMap<reco::PreIdRef>::Filler & filler);
  int getBin( float );


  // ----------member data ---------------------------
  ESHandle<MagneticField> magField;
  ESHandle<TrackerGeometry> pDD;
  edm::Handle<reco::BeamSpot> theBeamSpot;

  //  edm::InputTag PreIdMapLabel_;
  edm::InputTag gsfTrackLabel_;
  edm::InputTag pfNuclear_;
  edm::InputTag pfTrackLabel_;

  // YM additions
  Handle<GenParticleCollection> genParticles_;
  edm::InputTag label_tp_;
  edm::InputTag tracksTag;
  edm::InputTag simtracksTag;
  // YM end of my additions

  std::string   fitterName_;
  std::string   smootherName_;
  ///Fitter
  edm::ESHandle<TrajectoryFitter> fitter_;
  ///Smoother
  edm::ESHandle<TrajectorySmoother> smoother_;

  /// Map used to create the TrackRef, PreIdRef value map
  std::map<reco::TrackRef,unsigned> refMap_;

  edm::InputTag pfCLusTagECLabel_;
  edm::InputTag pfCLusTagHCLabel_;
  edm::InputTag pfCLusTagPSLabel_;
  std::vector<edm::InputTag> tracksContainers_;
  ///TRACK QUALITY
  bool useQuality_;
  reco::TrackBase::TrackQuality trackQuality_;
  Double_t      clusThreshold_;
  Double_t      minEp_;
  Double_t      maxEp_;
  Double_t      minPt_;
  Double_t      maxPt_;
  Double_t      maxEta_;
  Bool_t        usePreshower_;

  ///Vector of clusters of the PreShower
  std::vector<reco::PFCluster> ps1Clus;
  std::vector<reco::PFCluster> ps2Clus;
  Bool_t useNuclear_;

  ///ISOLATION REQUEST AS DONE IN THE TAU GROUP
  double HcalIsolWindow_;
  double EcalStripSumE_minClusEnergy_;
  double EcalStripSumE_deltaEta_;
  double EcalStripSumE_deltaPhiOverQ_minValue_;
  double EcalStripSumE_deltaPhiOverQ_maxValue_;
  double minEoverP_;
  double maxHoverP_;

  std::string filename_;
  ///B field
  math::XYZVector B_;
  PFTrackTransformer *pfTransformer_;
  edm::ParameterSet conf_;
  PFResolutionMap* resMapEtaECAL_;
  PFResolutionMap* resMapPhiECAL_;
  PFTrackTransformer* pfTkTransformer_;
};

ElePFAnalyzer::ElePFAnalyzer(const edm::ParameterSet& iConfig):
  pfTransformer_(0),
  conf_(iConfig),
  resMapEtaECAL_(0),
  resMapPhiECAL_(0),
  pfTkTransformer_(0)
{
  pfTrackLabel_     = iConfig.getParameter<edm::InputTag>( "PFRecTrackLabel" );
  pfNuclear_        = iConfig.getParameter<edm::InputTag>( "PFNuclear" );
  gsfTrackLabel_    = iConfig.getParameter<edm::InputTag>( "GsfTrackModuleLabel" );

  //YM get the input tags for the collections to read from the event
  label_tp_         = iConfig.getParameter<edm::InputTag>("label_tp");
  tracksTag         = iConfig.getParameter<edm::InputTag>("tracksTag");
  simtracksTag      = iConfig.getParameter<edm::InputTag>("simtracksTag");

  tracksContainers_ = iConfig.getParameter< vector < InputTag > >("TkColList");
  useQuality_       = iConfig.getParameter<bool>("UseQuality");
  trackQuality_     = TrackBase::qualityByName(iConfig.getParameter<std::string>("TrackQuality"));

  fitterName_       = iConfig.getParameter<string>("Fitter");
  smootherName_     = iConfig.getParameter<string>("Smoother");  
  pfCLusTagECLabel_ = iConfig.getParameter<InputTag>("PFEcalClusterLabel");
  pfCLusTagHCLabel_ = iConfig.getParameter<InputTag>("PFHcalClusterLabel");  
  pfCLusTagPSLabel_ = iConfig.getParameter<InputTag>("PFPSClusterLabel");
  clusThreshold_    = iConfig.getParameter<double>("ClusterThreshold");
  minEp_            = iConfig.getParameter<double>("MinEOverP");
  maxEp_            = iConfig.getParameter<double>("MaxEOverP");
  minPt_            = iConfig.getParameter<double>("MinPt");
  maxPt_            = iConfig.getParameter<double>("MaxPt");
  maxEta_           = iConfig.getParameter<double>("MaxEta");
  usePreshower_     = iConfig.getParameter<bool>("UsePreShower");
  filename_         = iConfig.getParameter<string>( "filename" );

}

ElePFAnalyzer::~ElePFAnalyzer()
{
  delete pfTransformer_;
  delete pfTkTransformer_;
  delete resMapEtaECAL_;
  delete resMapPhiECAL_;

  pfTransformer_ = nullptr;
  pfTkTransformer_ = nullptr;
  resMapEtaECAL_ = nullptr;
  resMapPhiECAL_ = nullptr;

}

bool ElePFAnalyzer::isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer) {
  // copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoParticleFlow/PFTracking/plugins/PFElecTkProducer.cc?revision=1.23&view=markup
  reco::HitPattern hitPattern1 = track1.hitPattern();
  reco::HitPattern hitPattern2 = track2.hitPattern();

  // retrieve the first valid hit
  int hitCounter1 = 0;
  trackingRecHit_iterator hitsIt1;
  for(hitsIt1 = track1.recHitsBegin(); hitsIt1 != track1.recHitsEnd(); ++hitsIt1, ++hitCounter1)
    { if (((**hitsIt1).isValid())) break; }
  int hitCounter2 = 0;
  trackingRecHit_iterator hitsIt2;
  for(hitsIt2 = track2.recHitsBegin(); hitsIt2 != track2.recHitsEnd(); ++hitsIt2, ++hitCounter2)
    { if (((**hitsIt2).isValid())) break; }
  uint32_t hit1 = hitPattern1.getHitPattern(hitCounter1);
  uint32_t hit2 = hitPattern2.getHitPattern(hitCounter2);
  if ( hitPattern1.getSubStructure(hit1) != hitPattern2.getSubStructure(hit2) ) 
    return hitPattern2.getSubStructure(hit2) < hitPattern1.getSubStructure(hit1);
  else if ( hitPattern1.getLayer(hit1) != hitPattern2.getLayer(hit2) ) 
    return hitPattern2.getLayer(hit2) < hitPattern1.getLayer(hit1);
  else {
    sameLayer = true;
    return false;
  }
}

// returns true if found, and the reco::Track track is returned
void ElePFAnalyzer::findRecoTracks(Handle<TrackingParticleCollection> TPCollectionH,
				  reco::SimToRecoCollection q,
				  std::vector<std::pair<const reco::Track, float> >& tracksAndSharedHits, 
				  int pdgId) {

  tracksAndSharedHits.clear();
  // get the collection
  const TrackingParticleCollection tPC = *(TPCollectionH.product());
  if ( tPC.size() == 0 ) return;

  unsigned int nTPs = tPC.size();
  // loop over tracking particles
  for (unsigned int iTP=0; iTP < nTPs; ++iTP){
    const TrackingParticle& lTP(tPC[iTP]);
    
    edm::Ref<TrackingParticleCollection> tp(TPCollectionH, iTP);
    if ( lTP.pt() > 2.0 && abs(lTP.pdgId()) == pdgId && fabs(lTP.eta()) < 2.6 ) {
      try{ 
	std::vector<std::pair<RefToBase<Track>, double> > trackV = q[tp];
	if ( trackV.size() == 1 ) {
	  std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
	  RefToBase<Track> tr = it->first;
	  Track track1 = *tr;
	  float sharedHits = it->second;
	  std::pair<const Track, float> newElement = make_pair(track1, sharedHits);
	  tracksAndSharedHits.push_back(newElement);
	}
	if ( trackV.size() > 1 ) {
	  std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
	  RefToBase<Track> tr = it->first;
	  Track track1 = *tr;
	  float nHitsRatio1 = it->second;
	  ++it; tr = it->first;
	  Track track2 = *tr;
	  float nHitsRatio2 = it->second;

	  bool sameLayer = false;
	  bool result = isInnerMost(track1, track2, sameLayer);
	  Track finalTrack;
	  float sharedHits = 0;
	  if ( !sameLayer ) {
	    if ( result ) {
	      //cout << "track2 is inner most, with pT = " << track2.pt() << '\t' << nHitsRatio2 << endl;
	      finalTrack = track2;
	      sharedHits = nHitsRatio2;
	    }
	    else {
	      //cout << "track1 is inner most, with pT = " << track1.pt() << '\t' << nHitsRatio1 << endl;
	      finalTrack = track1;
	      sharedHits = nHitsRatio1;
	    }
	  } else {
	    //cout << "track1 and track2 are in the same layer!" << endl;
	    if ( nHitsRatio2 > nHitsRatio1 ) {
	      // use track2
	      finalTrack = track2;
	      sharedHits = nHitsRatio2;
	    }
	    else {
	      // use track1
	      finalTrack = track1;
	      sharedHits = nHitsRatio1;
	    }
	  }
	  std::pair<const Track, float> newElement = make_pair(finalTrack, sharedHits);
	  tracksAndSharedHits.push_back(newElement);
	}

	// only for electrons
	if ( abs(pdgId) == 11 ) {
	  double dRMAX = 100.0;
	  int iel = -1;

	  for( size_t i = 0; i < genParticles_->size(); ++i ) {
	    const GenParticle& genCand = (*genParticles_)[i];
	    if ( abs(genCand.pdgId()) != 11 ) continue; 
	    double eta0 = lTP.eta();
	    double phi0 = lTP.phi();
	    double eta1 = genCand.eta();
	    double phi1 = genCand.phi();
	    double dphi = acos(cos(phi0 - phi1));
	    double deta = eta0 - eta1;
	    double dr = sqrt(deta*deta + dphi*dphi);
	    if ( dr < dRMAX ) {
	      iel = i;
	      dRMAX = dr;
	    }
	    if ( iel != -1 && dRMAX < 0.05 ) {
	      // check if this is from the B flavor
	      //cout << "Match with genCand is found at " << dRMAX << '\t' << lTP.pt() << " = " << (*genParticles_)[iel].pt() << '\t' << tracer((*genParticles_)[iel]) << endl;
	      ;
	    }
	    else
	      ;
	    //cout << "No match is found in the genCand block for TrackingParticle " << lTP.pt() << endl;
	  }
	}
      } catch (Exception event) {
	// no track is found to be matched, do nothing and move on to the next tracking particle
	;
      }
    }
  }
}
// ------------ method called for each event  ------------
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

bool ElePFAnalyzer::findMatch(const reco::Candidate& genCand, 
			      Handle<TrackingParticleCollection> TPCollectionH,
			      reco::SimToRecoCollection q,
			      reco::Track& track,
			      int& indexTrack,
			      float& dR,
			      float& sharedHits) {

  // get the collection
  const TrackingParticleCollection tPC = *(TPCollectionH.product());
  if ( tPC.size() == 0 ) return false;

  unsigned int nTPs = tPC.size();
  // loop over tracking particles
  int counter = 0;
  for(unsigned int iTP = 0; iTP < nTPs; ++iTP) {
    const TrackingParticle& lTP(tPC[iTP]);
    float deta = lTP.eta() - genCand.eta();
    float dphi = acos(cos(lTP.phi() - genCand.phi()));
    dR = sqrt(deta*deta + dphi*dphi);

    if ( fabs(lTP.pt() - genCand.pt()) > 0.1 || abs(lTP.pdgId()) != abs(genCand.pdgId()) || dR > 0.1 ) continue;
    ++counter;
    //if ( counter > 1 ) cout << "Found more than one tracking particle!!! __LINE__" << '\t' << counter << endl;
    edm::Ref<TrackingParticleCollection> tp(TPCollectionH, iTP);
    try {
      std::vector<std::pair<RefToBase<Track>, double> > trackV = q[tp];
      if ( trackV.size() == 1 ) {
	std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
	RefToBase<Track> tr = it->first;
	indexTrack = tr.key();
	track = *tr;
	sharedHits = it->second;
	return true;
      } 
      if ( trackV.size() > 1 ) {
	std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
	RefToBase<Track> tr = it->first;
	Track track1 = *tr;
	int indexTrack1 = tr.key();
	float nHitsRatio1 = it->second;
	++it; tr = it->first;
	int indexTrack2 = tr.key();
	Track track2 = *tr;
	float nHitsRatio2 = it->second;
	bool sameLayer = false;
	bool result = isInnerMost(track1, track2, sameLayer);
	Track finalTrack;
	if ( !sameLayer ) {
	  if ( result ) {
	    finalTrack = track2;
	    sharedHits = nHitsRatio2;
	    indexTrack = indexTrack2;
	  }
	  else {
	    //cout << "track1 is inner most, with pT = " << track1.pt() << '\t' << nHitsRatio1 << endl;
	    finalTrack = track1;
	    sharedHits = nHitsRatio1;
	    indexTrack = indexTrack1;
	  }
	} else {
	  //cout << "track1 and track2 are in the same layer!" << endl;
	  if ( nHitsRatio2 > nHitsRatio1 ) {
	    // use track2
	    finalTrack = track2;
	    sharedHits = nHitsRatio2;
	    indexTrack = indexTrack2;
	  }
	  else {
	    // use track1
	    finalTrack = track1;
	    sharedHits = nHitsRatio1;
	    indexTrack = indexTrack1;
	  }
	}
	track = finalTrack;
	return true;
      }
    } catch (Exception event) {
      // no track is found to be matched, do nothing
      return false;
    }
  }
  return false;
}

int ElePFAnalyzer::tracer(const reco::Candidate& genCand) {
  int parentId = 0;
  if ( genCand.numberOfMothers() != 0 ){
    const reco::Candidate& mom0 = *genCand.mother(0);
    parentId = parent(mom0);
    
    if ( parentId == 0 && !(mom0.pdgId() == 111 || mom0.pdgId()==221) ){
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
  return parentId;
}

void ElePFAnalyzer::FillSeedInformation(const reco::PFRecTrackCollection & PfRTkColl,
				       reco::ElectronSeed EleSeed,
				       bool otherColl,
				       const edm::Event& iEvent,
				       const EventSetup& iSetup){
   
  const reco::BeamSpot& bs = *(theBeamSpot.product());
  TSCBLBuilderNoMaterial tscblBuilder;
  const GeomDet *detc = 0;

  TrajectorySeed::const_iterator hit_end = EleSeed.recHits().second;

  //You don't need to loop over for hits. The good direction is given by the latest hit
  detc = pDD->idToDet( (*(--hit_end)).geographicalId() );

  TrajectoryStateOnSurface t = trajectoryStateTransform::transientState((EleSeed).startingState(), &(detc->surface()), &(*magField));
  const FreeTrajectoryState & stateForProjectionToBeamLine = *t.freeState();
  TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine, bs);
  GlobalPoint v = tscbl.trackStateAtPCA().position();
  math::XYZPoint pos( v.x(), v.y(), v.z() );
  GlobalVector p = tscbl.trackStateAtPCA().momentum();
  math::XYZVector mom( p.x(), p.y(), p.z() );

  float seedPt = sqrt(p.x()*p.x() + p.y()*p.y());

  seed_pt.push_back(seedPt);
  seed_eta.push_back(mom.eta());
  seed_phi.push_back(mom.phi());

}// end of fill seed info


int ElePFAnalyzer::parent(const reco::Candidate& mom){
  int pdgId = mom.pdgId();

  if ( isaV      (fabs(pdgId)) ) return 1;
  if ( isaBhadron(fabs(pdgId)) ) return 2;
  if ( isaDhadron(fabs(pdgId)) ) return 3;
      
  return 0;
}


void ElePFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  // Handles, collections etc. ===========================================
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", theHitsAssociator);
  TrackAssociatorBase* associatorByHits = (TrackAssociatorBase *)theHitsAssociator.product();

  iSetup.get<IdealMagneticFieldRecord>().get(magField);
  iSetup.get<TrackerDigiGeometryRecord>().get(pDD);

  InputTag beamSpot_(string("offlineBeamSpot"));
  iEvent.getByLabel(beamSpot_, theBeamSpot);

  // YM
  Handle<reco::TrackCollection> theTrackCollection;
  iEvent.getByLabel(tracksTag, theTrackCollection);
  reco::TrackCollection  Tk=*(theTrackCollection.product());

  Handle<View<Track> > trackCollectionH;
  iEvent.getByLabel(tracksTag, trackCollectionH);
  const View<Track> tC = *(trackCollectionH.product());

  Handle<vector<Trajectory> > tjCollection;
  iEvent.getByLabel(tracksTag, tjCollection);
  vector<Trajectory> Tj=*(tjCollection.product());

  Handle<SimTrackContainer> simTrackCollection;
  iEvent.getByLabel(simtracksTag, simTrackCollection);
  const SimTrackContainer simTC = *(simTrackCollection.product());

  Handle<TrackingParticleCollection> TPCollectionH;
  iEvent.getByLabel(label_tp_, TPCollectionH);
  const TrackingParticleCollection tPC = *(TPCollectionH.product());

  if ( tPC.size() == 0 ) {
    edm::LogInfo("ElePFAnalyzer") 
      << "TP Collection for efficiency studies has size = 0! Skipping Event.";
    return;
  }
  iEvent.getByLabel("genParticles", genParticles_);
  
  Handle<reco::GsfTrackCollection> gsftrackscoll;
  iEvent.getByLabel(gsfTrackLabel_,gsftrackscoll);

  Handle<reco::PFRecTrackCollection> thePfRecTrackCollection;
  iEvent.getByLabel("pfTrack",thePfRecTrackCollection);
  const PFRecTrackCollection& PfRTkColl = *(thePfRecTrackCollection.product());

  //ECAL clusters	      
  Handle<PFClusterCollection> theECPfClustCollection;
  iEvent.getByLabel(pfCLusTagECLabel_,theECPfClustCollection);
  vector<PFCluster> basClus;
  vector<PFCluster>::const_iterator iklus;
  for (iklus=theECPfClustCollection.product()->begin();
       iklus!=theECPfClustCollection.product()->end();
       iklus++){
    if((*iklus).correctedEnergy()>clusThreshold_) basClus.push_back(*iklus);
  }

  iSetup.get<TrajectoryFitter::Record>().get(fitterName_, fitter_);
  iSetup.get<TrajectoryFitter::Record>().get(smootherName_, smoother_);

  // clear temporary maps
  refMap_.clear();

  //HCAL clusters
  Handle<PFClusterCollection> theHCPfClustCollection;
  iEvent.getByLabel(pfCLusTagHCLabel_,theHCPfClustCollection);
  
  //PS clusters
  Handle<PFClusterCollection> thePSPfClustCollection;
  iEvent.getByLabel(pfCLusTagPSLabel_,thePSPfClustCollection);

  ps1Clus.clear();
  ps2Clus.clear();

  for (iklus=thePSPfClustCollection.product()->begin();
       iklus!=thePSPfClustCollection.product()->end();
       iklus++){
    if ((*iklus).layer()== PFLayer::PS1) ps1Clus.push_back(*iklus);
    if ((*iklus).layer()== PFLayer::PS2) ps2Clus.push_back(*iklus);
  }

  reco::SimToRecoCollection q = 
    associatorByHits->associateSimToReco(trackCollectionH, 
  					 TPCollectionH, 
  					 &iEvent,
  					 &iSetup);


  Handle<reco::ElectronSeedCollection> seedelectrons;
  iEvent.getByLabel("electronMergedSeeds", seedelectrons);

  ElectronSeedCollection eleSeeds = *( seedelectrons.product() );
  nSeed = 0;
  for ( unsigned int iseed = 0; iseed < eleSeeds.size(); ++iseed ){
    ElectronSeedRef SeedRef( seedelectrons, iseed );

    //reco-seed information
    FillSeedInformation(PfRTkColl, eleSeeds[iseed], false, iEvent, iSetup);

    bool isTrackerDrivenForThisSeed = true;
    bool isEcalDrivenForThisSeed = true;

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

  // ======================================================================
  // Main loop: over gen particles in the event
  nEle = 0;
  nPion = 0;
  for( size_t i = 0; i < genParticles_->size(); ++i ){
    const GenParticle& genCand = (*genParticles_)[i];

    bool electron = abs(genCand.pdgId()) == 11;
    bool pion     = abs(genCand.pdgId()) == 211;

    if ( genCand.pt() < 2.0 || (!(electron || pion)) ) continue; 

    // check if it is pion - we do not want any electrons around
    if ( pion ) {
      float mindR2 = 1000;
      for( vector<unsigned int>::iterator it = iGenElectrons.begin(); it != iGenElectrons.end(); ++it ) {
	const GenParticle& electron = (*genParticles_)[*it];
	float dphi = acos(cos(genCand.phi() - electron.phi()));
	float deta = genCand.eta() - electron.eta();
	float dr2 = dphi*dphi + deta*deta;
	if ( dr2 < mindR2 ) mindR2 = dr2;
      }
      if ( sqrt(mindR2) < 0.3 ) continue;
    }

    // consider electrons and pions for pT > 2 GeV
    int parentId = tracer(genCand);

    // we fill the tree per each gen-level particle
    // so clear the variables
    clearBDTTreeVariables();

    // fill gen level info for a particle:
    gen_pt_ = genCand.pt();
    gen_phi_ = genCand.phi();
    gen_eta_ = genCand.eta();
    gen_pdgId_ = genCand.pdgId();
    gen_parentId_ = parentId;

    if (electron == 1){
      gen_ele_pt.push_back(genCand.pt());
      gen_ele_phi.push_back(genCand.phi());
      gen_ele_eta.push_back(genCand.eta());
      gen_ele_status.push_back(genCand.status());
      gen_ele_id.push_back(genCand.pdgId());
      gen_ele_isFromMom.push_back(parentId);
      ++nEle;
    }
    if (pion == 1){
      gen_pion_pt.push_back(genCand.pt());
      gen_pion_phi.push_back(genCand.phi());
      gen_pion_eta.push_back(genCand.eta());
      gen_pion_status.push_back(genCand.status());
      gen_pion_id.push_back(genCand.pdgId());
      gen_pion_isFromMom.push_back(parentId);
      ++nPion;
    }

    // find the nearest electron
    gen_dRElectron_ = 1000;
    for( size_t imc = 0; imc < genParticles_->size(); ++imc ){
      if ( imc == i ) continue; // do not consider itself :)
      const GenParticle& genCandToCheck = (*genParticles_)[imc];
      if ( abs(genCandToCheck.pdgId()) == 11 ) {
	double deta = genCandToCheck.eta() - gen_eta_;
	double dphi = acos(cos(genCandToCheck.phi() - gen_phi_));
	double dr2 = deta*deta + dphi*dphi;
	if (dr2 < gen_dRElectron_ ) gen_dRElectron_ = dr2;
      }
    }
    if ( gen_dRElectron_ != 1000 ) gen_dRElectron_ = sqrt(gen_dRElectron_);

    // find the nearest electron
    sim_dRElectron_ = 1000;
    for( size_t imc = 0; imc < simTC.size(); ++imc ){
      SimTrack simtrack = simTC[imc];
      if ( abs(simtrack.type()) == 11 && simtrack.momentum().pt() > 1.0 ) {
	double deta = simtrack.momentum().eta() - gen_eta_;
	double dphi = acos(cos(simtrack.momentum().phi() - gen_phi_));
	double dr2 = deta*deta + dphi*dphi;
	if (dr2 < sim_dRElectron_ ) sim_dRElectron_ = dr2;
      }
    }
    if ( sim_dRElectron_ != 1000 ) sim_dRElectron_ = sqrt(sim_dRElectron_);

    Track track;
    float sharedHits = 0;
    float dR = 100;
    int indexTrack;
    Track trRefToBase;
    bool result = findMatch(genCand, TPCollectionH, q, 
			    track, indexTrack, dR, sharedHits);
    if ( result ) gen_match_ = 1;
    else          gen_match_ = 0;

    //cout << genCand.pdgId() << '\t' << genCand.pt(); 
    if ( result ) {
      //cout << '\t' << "matched!" << '\t' 
      //<< track.pt() << '\t' << dR << '\t' << sharedHits << endl;
      float px_in = track.innerMomentum().x();
      float py_in = track.innerMomentum().y();
      float pt_in = sqrt(px_in*px_in + py_in*py_in);
      float dPt = (pt_in - track.outerPt())/pt_in;

      trk_pt_       = track.pt();
      trk_pTOB_     = Tj[indexTrack].lastMeasurement().updatedState().globalMomentum().mag();
      trk_eta_      = track.eta();
      trk_phi_      = track.phi();
      trk_chi2_     = track.chi2();
      trk_ndof_     = track.ndof();
      trk_nchi2_    = track.normalizedChi2();
      trk_dpt_      = dPt;
      trk_nhits_    = track.found();
      trk_nmatched_ = sharedHits;
      trk_quality_  = track.quality(trackQuality_);

      TrackRef trackRef(theTrackCollection, indexTrack);
      reco::PFRecTrack pfrectrack( trackRef->charge(), 
				   reco::PFRecTrack::KF, 
				   indexTrack, trackRef);

      Trajectory FakeTraj;
      bool valid = pfTkTransformer_->addPoints(pfrectrack, *trackRef, FakeTraj);
      pfrectrack.calculatePositionREP();

      // momentum and position of the track for old style matching
      float pfmass=  0.0005; 
      float pfoutenergy=sqrt((pfmass*pfmass)+Tk[indexTrack].outerMomentum().Mag2());
      XYZTLorentzVector mom =XYZTLorentzVector(Tk[indexTrack].outerMomentum().x(),
					       Tk[indexTrack].outerMomentum().y(),
					       Tk[indexTrack].outerMomentum().z(),
					       pfoutenergy);

      XYZTLorentzVector pos =   XYZTLorentzVector(Tk[indexTrack].outerPosition().x(),
						  Tk[indexTrack].outerPosition().y(),
						  Tk[indexTrack].outerPosition().z(),
						  0.);
      
      BaseParticlePropagator theOutParticle = 
	BaseParticlePropagator( RawParticle(mom,pos),
				0,0,B_.z());
      theOutParticle.setCharge(Tk[indexTrack].charge());
      
      theOutParticle.propagateToEcalEntrance(false);
      
      float toteta=1000;
      float totphi=1000;
      float dr=1000;
      float EE = 0;
      float feta=0;
      math::XYZPointF ElecTrkEcalPos(0,0,0);
      PFClusterRef clusterRef;
      math::XYZPoint meanShowerSaved;
      if(theOutParticle.getSuccess()!=0){
	ElecTrkEcalPos=math::XYZPointF(theOutParticle.vertex().x(),
				       theOutParticle.vertex().y(),
				       theOutParticle.vertex().z());
	bool isBelowPS=(fabs(theOutParticle.vertex().eta())>1.65) ? true :false;	

	unsigned clusCounter=0;	
	float max_ee = 0;
	for(vector<PFCluster>::const_iterator aClus = basClus.begin();
	    aClus != basClus.end(); aClus++,++clusCounter) {
	  
	  double ecalShowerDepth
	    = PFCluster::getDepthCorrection(aClus->correctedEnergy(),
					    isBelowPS,
					    false);
	  
	  math::XYZPoint meanShower=math::XYZPoint(theOutParticle.vertex())+
	    math::XYZTLorentzVector(theOutParticle.momentum()).Vect().Unit()*ecalShowerDepth;	
	  
	  float etarec=meanShower.eta();
	  float phirec=meanShower.phi();
	  float tmp_ep=aClus->correctedEnergy()/trk_pTOB_;
	  float tmp_phi=fabs(aClus->position().phi()-phirec);
	  if (tmp_phi>TMath::TwoPi()) tmp_phi-= TMath::TwoPi();
	  float tmp_dr=sqrt(pow(tmp_phi,2)+
			    pow((aClus->position().eta()-etarec),2));
	  

	  if ((tmp_dr<dr)&&(tmp_ep>minEp_)&&(tmp_ep<maxEp_)){
	    dr=tmp_dr;
	  
	    if(dr < 0.2){

	      if(aClus->correctedEnergy() > max_ee){
		max_ee = aClus->correctedEnergy();

		toteta=aClus->position().eta()-etarec;
		totphi=tmp_phi;
		EE = aClus->correctedEnergy();
		feta= aClus->position().eta();
		meanShowerSaved = meanShower;

		clusterRef = PFClusterRef(theECPfClustCollection,clusCounter);
	      }
	    }
	  } // loop over dr
	} // loop over basic clusters
      } // out particle success
      trk_dr_ = dr;

      double ecaletares = resMapEtaECAL_->GetBinContent(resMapEtaECAL_->FindBin(feta,EE));
      double ecalphires = resMapPhiECAL_->GetBinContent(resMapPhiECAL_->FindBin(feta,EE));

      // geometrical compatibility
      float chieta = (toteta!=1000) ? toteta/ecaletares : toteta;
      float chiphi = (totphi!=1000) ? totphi/ecalphires : totphi;
      float chichi = sqrt(chieta*chieta + chiphi*chiphi); 
      trk_ecalChi_ = chichi;
      trk_ecalDeta_ = fabs(toteta);
      trk_ecalDphi_ = fabs(totphi);

      //       // Check for PS clusters
      //       std::vector<std::pair<PFCluster, float> > psClustersAndDist;
      //       psClustersAndDist.clear();
      //       float min_dist = -1.0; // minimal distance to consider in PS matching
      //       if ( fabs(feta) > 1.68 ) { // in PS area
      // 	for ( PFClusterCollection::const_iterator ips = ps1Clus.begin();
      // 	      ips != ps1Clus.end(); ++ips ) {
      // 	  float dist = testPreshowerDistance(*clusterRef, *ips);
      // 	  if ( dist == -1.0 || ( min_dist != -1.0 && dist > min_dist ) ) continue;
      // 	  if ( dist < min_dist || min_dist == -1.0 ) {
      // 	    psClustersAndDist.push_back(std::make_pair(*ips, dist));
      // 	  }
      // 	}

      // 	for ( PFClusterCollection::const_iterator ips = ps2Clus.begin();
      // 	      ips != ps2Clus.end(); ++ips ) {
      // 	  float dist = testPreshowerDistance(*clusterRef, *ips);
      // 	  if ( dist == -1.0 || ( min_dist != -1.0 && dist > min_dist ) ) continue;
      // 	  if ( dist < min_dist || min_dist == -1.0 ) {
      // 	    psClustersAndDist.push_back(std::make_pair(*ips, dist));
      // 	  }
      // 	}
	
      // 	for(std::vector<std::pair<PFCluster, float> >::const_iterator ips = psClustersAndDist.begin();
      // 	    ips != psClustersAndDist.end(); ++ips) {
      // 	  ++ecal_ps_;
      // 	  PFCluster ps = ips->first;
      // 	  const PFLayer::Layer pslayer = ps.layer();
      // 	  const double psenergy = ps.energy();
      // 	  ecal_ps1e_ += (PFLayer::PS1 == pslayer)*psenergy;
      // 	  ecal_ps2e_ += (PFLayer::PS2 == pslayer)*psenergy;
      // 	}
      //       }//loop over for eta cut |eta| > 1.68

      // EoP ----------
      ecal_e_ = EE;
      trk_ep_ = EE/trk_pTOB_;

      Trajectory::ConstRecHitContainer tmp;
      TrajectorySeed Seed = (*trackRef->seedRef());
      
      Trajectory::ConstRecHitContainer hits=Tj[indexTrack].recHits();
      for (int ih=hits.size()-1; ih>=0; ih--)  tmp.push_back(hits[ih]);
      vector<Trajectory> FitTjs=(fitter_.product())->fit(Seed,
							 tmp,
							 Tj[indexTrack].lastMeasurement().updatedState());
      
      if(FitTjs.size()>0){
	if(FitTjs[0].isValid()){
	  vector<Trajectory> SmooTjs=(smoother_.product())->trajectories(FitTjs[0]);
	  if(SmooTjs.size()>0){
	    if(SmooTjs[0].isValid()){
	      //Track refitted with electron hypothesis
	      float pt_out = SmooTjs[0].firstMeasurement().updatedState().globalMomentum().perp();
	      float pt_in = SmooTjs[0].lastMeasurement().updatedState().globalMomentum().perp();
	      float dpt = (pt_in>0) ? fabs(pt_out-pt_in)/pt_in : 0.;
	      trk_dptGSF_ = dpt;
	      trk_chiRatio_ = SmooTjs[0].chiSquared()/Tj[indexTrack].chiSquared();
	      trk_chiReduced_ = trk_chiRatio_*track.normalizedChi2();
	    }
	  }//SmooTjs.size
	}//loop over for FitTjs.isValid()
      }//loop over for FitTjs.size()

      // check if the track is matched with the ECAL
      const reco::PFTrajectoryPoint& atECAL = 
	pfrectrack.extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);

      if ( valid && atECAL.isValid()) {
	
	float minDist = 1000;
	PFCluster matchedCluster;
	for(vector<PFCluster>::const_iterator aClus = basClus.begin();
	    aClus != basClus.end(); ++aClus) {
	  
	  PFCluster clust = *aClus;
	  clust.calculatePositionREP();
	  float dist = LinkByRecHit::testTrackAndClusterByRecHit(pfrectrack, clust);
	  
	  if ( dist > 0. && dist < minDist ){
	    minDist = dist;
	    matchedCluster = clust;
	  }
	}//loop over for aClus
      }//loop over for atECAL
    }//loop over result
    else {
      //cout << endl;
      ;
    }
  }//loop over genparticle
}


    // ------------ method called once each job just before starting event loop  ------------
    void ElePFAnalyzer::beginJob()
    {
      ;
    }

    // ------------ method called once each job just after ending the event loop  ------------
    void ElePFAnalyzer::endJob() 
    {
      file->cd();
      seeds->Write();
      bdtTree->Write();
      EventInfo->Write();
      file->Write();
    }

    // ------------ method called when starting to processes a run  ------------

    void ElePFAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup & es)
    {
  
      //Magnetic Field
      ESHandle<MagneticField> magneticField;
      es.get<IdealMagneticFieldRecord>().get(magneticField);
      B_=magneticField->inTesla(GlobalPoint(0,0,0));
  
      pfTransformer_= new PFTrackTransformer(B_);
      pfTransformer_->OnlyProp();
  
      pfTkTransformer_= new PFTrackTransformer(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));
      pfTkTransformer_->OnlyProp();  

      //Resolution maps
      FileInPath ecalEtaMap(conf_.getParameter<string>("EtaMap"));
      FileInPath ecalPhiMap(conf_.getParameter<string>("PhiMap"));
      resMapEtaECAL_ = new PFResolutionMap("ECAL_eta",ecalEtaMap.fullPath().c_str());
      resMapPhiECAL_ = new PFResolutionMap("ECAL_phi",ecalPhiMap.fullPath().c_str());

  
      //read PS threshold
      FileInPath parPSFile(conf_.getParameter<string>("PSThresholdFile"));
      ifstream ifsPS(parPSFile.fullPath().c_str());
      //      for (int iy=0;iy<12;iy++) ifsPS >> thrPS[iy];
  
    }

    // ------------ method called when ending the processing of a run  ------------

    void  ElePFAnalyzer::endRun(const edm::Run & run, const edm::EventSetup & es)
    {
      delete pfTransformer_;
      pfTransformer_ = nullptr;
      delete pfTkTransformer_;
      pfTkTransformer_ = nullptr;

    }
int ElePFAnalyzer::getBin(float pt){
  int ip=0;
  if (pt<6) ip=0;
  else {  if (pt<12) ip=1;
    else ip=2;
  }
  return ip;
}  

//     double ElePFAnalyzer::testPreshowerDistance(PFCluster eeclus,
// 						PFCluster psclus) {

//       const reco::PFCluster::REPPoint& pspos = psclus.positionREP();
//       const reco::PFCluster::REPPoint& eepos = eeclus.positionREP();
//       // same side of the detector?
//       if ( eeclus.z()*psclus.z() < 0 ) return -1.0;

//       const double dphi = std::abs(TVector2::Phi_mpi_pi(eepos.phi() - pspos.phi()));
//       if ( dphi > 0.6 ) return -1.0;
//       const double deta = std::abs(eepos.eta() - pspos.eta());
//       if ( deta > 0.3 ) return -1.0;
//       return LinkByRecHit::testECALAndPSByRecHit(eeclus, psclus, false);
//     }


    void ElePFAnalyzer::fillPreIdRefValueMap( Handle<TrackCollection> tracks,
					      const edm::OrphanHandle<reco::PreIdCollection>& preidhandle,
					      edm::ValueMap<reco::PreIdRef>::Filler & filler){
      std::vector<reco::PreIdRef> values;

      unsigned ntracks=tracks->size();
      for(unsigned itrack=0;itrack<ntracks;++itrack)
	{
	  reco::TrackRef theTrackRef(tracks,itrack);
	  std::map<reco::TrackRef,unsigned>::const_iterator itcheck=refMap_.find(theTrackRef);
	  if(itcheck==refMap_.end()) 
	    {
	      // the track has been early discarded
	      values.push_back(reco::PreIdRef());
	    }
	  else
	    {
	      edm::Ref<reco::PreIdCollection> preIdRef(preidhandle,itcheck->second);
	      values.push_back(preIdRef);
	      // std::cout << " Checking Refs " << (theTrackRef==preIdRef->trackRef()) << std::endl;
	    }
	}
      filler.insert(tracks,values.begin(),values. end());
    }


    //define this as a plug-in
    DEFINE_FWK_MODULE(ElePFAnalyzer);
  

