// system include files 
#include <memory>
#include <vector>
#include <map>
#include <cmath>
#include <cstring>
#include <algorithm>

// Framework Core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Muon
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//Tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

//Tracking Particles, Match to simulated tracks
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"

//SimTracks
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
//SimVertex
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

//Generation level particles
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//pdgTable
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

//CSC Digis
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

/*#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
*/

//Geometry
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

//ROOT
#include "TTree.h"
#include "TFile.h"

using namespace std;

class CSCPriEff : public edm::EDFilter {
   public:
      explicit CSCPriEff(const edm::ParameterSet&);
      ~CSCPriEff();
      enum ParticleType	{LightMeson=1,CharmedMeson=2,ccbarMeson=3,BottomMeson=4,bbarMeson=5,LightBaryon=6,CharmedBaryon=7,BottomBaryon=8,DiQuarks=9,Lepton=10,W=11,Z=12,Other=13};
      enum TheMuonType {PromptMuFromWZ=10,PromptMuFromLightMeson=11,PromptMuFromHeavyMeson=12,PromptMuFromLightBaryon=13,PromptMuFromHeavyBaryon=14,NotPromptMufromWZ=15,PromptMuFromOthers=16,PunchThrough=20,PunchThroughAndDecayInFlight=21,DecayInFlightFromLightMeson=31,DecayInFlightFromHeavyMeson=32,DecayInFlightFromLightBaryon=33,DecayInFlightFromHeavyBaryon=34,NoMuSysHit=40,Others=1,NothingMatched=0};
      //Details:
      //NotPromptMufromWZ=15,PromptMuFromOthers=16 are error cases
      //PunchThrough: it is a hadron and doesn't decay in the detector region
      //PunchThroughAndDecayInFlight: it is a hadron and decays outside of HCAL (including HO)
      //NoMuSysHit: Candidates can be hardon, electron.
   private:
      virtual void beginJob();
      virtual Bool_t filter(edm::Event&, const edm::EventSetup&);
      virtual reco::Muon::ArbitrationType MuonArbitrationTypeFromString( const std::string &);
      virtual void HepMCParentTree(HepMC::GenParticle *);
      virtual void SimTrackDaughtersTree(const SimTrack * thisTrk, TrackingParticleRef tpr);
      inline Bool_t SimHitToSegment(const PSimHit& hit);
      inline ParticleType ParticleCata(int);
      inline void ClearVecs();
      inline Bool_t IstheSameDChain(vector<int> &,vector<int> &);
      virtual TheMuonType classification(vector<Int_t> &);
      const TrackingRecHit* getHitPtr(trackingRecHit_iterator iter) const {return &**iter;}
      inline Int_t FindSimTrackRecordingPosition( Int_t SavedSimTrkPos ) {
	Int_t count=0;
	vector<Bool_t>::const_iterator IsParSavedAsHep_iter = IsParSavedAsHep.begin();
	for (; IsParSavedAsHep_iter != IsParSavedAsHep.end(); IsParSavedAsHep_iter++ ) {
	  if (!*IsParSavedAsHep_iter) {
	    if (count==SavedSimTrkPos) break;
	    else count++;
	  }
	}
	return IsParSavedAsHep_iter-IsParSavedAsHep.begin();
      }
      inline Int_t FindHepMCRecordingPosition( Int_t SavedHepParPos ) {
	Int_t count=0;
	vector<Bool_t>::const_iterator IsParSavedAsHep_iter = IsParSavedAsHep.begin();
	for (; IsParSavedAsHep_iter != IsParSavedAsHep.end(); IsParSavedAsHep_iter++ )
	  if (*IsParSavedAsHep_iter) {
	    if (count==SavedHepParPos) break;
	    else count++;
	  }
	return IsParSavedAsHep_iter-IsParSavedAsHep.begin();
      }
      inline Bool_t IsTheSameSegment(UInt_t Seg1,UInt_t Seg2) {
	if ( (*XSegment)[Seg1]!=(*XSegment)[Seg2] || (*YSegment)[Seg1]!=(*YSegment)[Seg2] ) return false;
	if ( (*Ring)[Seg1] != (*Ring)[Seg2] || (*Chamber)[Seg1] != (*Chamber)[Seg2] ) return false;
	return true;
      }
      virtual void endJob() ;

      
//-------------member data------------
   //TTree&TFile
   TFile *file;
   string FileName;
   TTree *Tracks_Tree,*Summarization_Tree,*ErrorMsg_Tree;
   struct General_EvtInfo
   {
     ULong64_t RUN,EVENT,LS,ORBIT,BX;
   } Info;
   //Tracker Muon
   vector<Float_t> *pt,*eta,*phi,*Vertex_x,*Vertex_y,*Vertex_z,*isoR03sumPt,*isoR03emEt,*isoR03hadEt,*isoR03hoEt,*isoR03nJets,*isoR03nTracks,*isoR05sumPt,*isoR05emEt,*isoR05hadEt,*isoR05hoEt,*isoR05nJets,*isoR05nTracks,*isoemVetoEt,*isohadVetoEt,*isohoVetoEt,*normalizedChi2,*TrkRelChi2,*CaloE_emMax,*CaloE_emS9,*CaloE_emS25,*CaloE_hadMax,*CaloE_hadS9,*Calo_emPos_R,*Calo_emPos_eta,*Calo_emPos_phi,*Calo_hadPos_R,*Calo_hadPos_eta,*Calo_hadPos_phi,*dEdx,*dEdxError,*TrkKink,*GlbKink;
   vector<Bool_t> *chargeMinus,*isGlobalMu,*isTrackerMu;
   vector<Int_t> *dEdx_numberOfSaturatedMeasurements,*dEdx_numberOfMeasurements;
   //Tracks
   edm::InputTag tracksTag;
   string dEdxTag;
   vector<UInt_t> *InnerTrack_nValidTrackerHits,*InnerTrack_nValidPixelHits,*InnerTrack_nLostTrackerHits,*InnerTrack_nLostPixelHits,*InnerTrack_ndof,*GlobalTrack_ndof;
   vector<Float_t> *InnerTrack_chi2;
   UInt_t minTrackHits;
   //Muon Selectors
   unsigned int num_Cuts;
   reco::Muon::ArbitrationType MuonArbitrationType;
   vector<muon::SelectionType> OfficialMuonSelectors;
   vector<Bool_t> * SelectorPassed[30],*MySelector;
   //CSC Digis
   edm::InputTag CSCDigisTag;
   //Muon Chamber and Segment Match  0-3=CSC station 1-4; save position: Mu*4+station-1
   vector<Float_t> *TrackDistToChamberEdge,*TrackDistToChamberEdgeErr,*XTrack,*YTrack,*XErrTrack,*YErrTrack,*XSegment,*YSegment,*XErrSegment,*YErrSegment,*DRTrackToSegment,*DRErrTrackToSegment,*XSimSegment,*YSimSegment;
   vector<Bool_t> *IsSegmentOwnedExclusively,*IsSegmentBestInStationByDR,*IsSegmentBelongsToTrackByDR,*IsSegmentBelongsToTrackByCleaning;
   vector<UInt_t> *StationMask,*RequiredStationMask;
   vector<Byte_t> *Chamber,*ChamberSimSegment,*MuonIndex,*MuonIndexSimSegment,*NumberOfLCTsInChamber,*NumberOfHitsInSegment;
   vector<Char_t> *Ring,*RingSimSegment;
   struct CSCSimHitInfo {
     Float_t X,Y;
     Byte_t Layer;
   };
   map<Int_t,CSCSimHitInfo> FirstLayer,LastLayer;//the first and last layer SimHits in one chamber. Int_t is the ChamberID
   double maxChamberDist,maxChamberDistPull;
   //PrimaryVertex
   string PrimaryVerticesTag;
   vector<Float_t> *vx,*vxError,*vy,*vyError,*vz,*vzError;
   //Tracking particles ( combination of simulated tracks )
   edm::InputTag TPInputTag;
   vector<Float_t> *TrkParticles_pt,*TrkParticles_eta,*TrkParticles_phi;
   vector<Double_t> *SharedHitsRatio,*MCMatchChi2;
   vector<Int_t> *TrkParticles_pdgId,*TrkParticles_charge;
   //To Get the Decay Chains
      //saved variables
   Float_t HepMCFourVec[4];
   vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*Gen_vx,*Gen_vy,*Gen_vz,*Gen_vt;
   vector<Int_t> *Gen_pdgId;
   vector<Long64_t> *MuonType;//all possible chains; -1 is at the beginning of each chain; -2is at the beginning of each reco-muons;
   vector<Bool_t> *IsParInHep,*IsParHasMuonHits;
   vector<Int_t> *DChains,*theSameWithMuon;
   //PSimHits
   //temparory variables
   string HepMCTag;
   bool ChainRecord;
   vector< vector<Int_t> > SimChains;
   vector< vector<Int_t> > HepMCChains;
   vector<Int_t> DChain;
   vector<const SimTrack *> SavedSimTrk;
   vector<HepMC::GenParticle *> SavedHepPar;
   vector<Bool_t> IsParSavedAsHep;
   vector<SimVertex> SVC;

   //Summarization
   struct Summary
   {
     ULong64_t Total_Events,Total_TrackerMuons,Total_GlobalMuon,Total_GlobalnotTrackerMuon;
   }Summarization;
   //Errors and Warnings reports
   struct ErrorMsg
   {
     ULong64_t ErrorCode,Run,Event;
   }Error;
   //Others
   static bool FirstEntry;
   ULong64_t NumMisMatch;
};
// ------------ method called once each job just before starting event loop  ------------

Bool_t CSCPriEff::FirstEntry=true;

void 
CSCPriEff::beginJob(){
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CSCPriEff::endJob() {
}
