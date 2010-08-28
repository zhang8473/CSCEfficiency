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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
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

/*
//CSC rechit
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

//Geometry
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include <Geometry/CSCGeometry/interface/CSCLayer.h>
#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
*/

//ROOT
#include "TTree.h"
#include "TFile.h"

using namespace std;

class CSCPriEff : public edm::EDFilter {
   public:
      explicit CSCPriEff(const edm::ParameterSet&);
      ~CSCPriEff();
      enum ParticleType	{LightMeson=1,CharmedMeson=2,ccbarMeson=3,BottomMeson=4,bbarMeson=5,LightBaryon=6,CharmedBaryon=7,BottomBaryon=8,DiQuarks=9,Lepton=10,W=11,Z=12,Other=13};

   private:
      virtual void beginJob();
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void HepMCParentTree(HepMC::GenParticle *);
      virtual void SimTrackDaughtersTree(SimTrack *);
      inline ParticleType ParticleCata(int);
      inline void ClearVecs();
      inline bool IstheSameDChain(vector<int> &,vector<int> &);
      virtual Byte_t classification(vector<Int_t> &);
      const TrackingRecHit* getHitPtr(trackingRecHit_iterator iter) const {return &**iter;}
      inline Int_t FindSimTrackRecordingPosition( Int_t ParToSimPos ) {
	Int_t count=0;
	vector<Bool_t>::const_iterator IsParInHep_iter = IsParInHep->begin();
	for (; IsParInHep_iter != IsParInHep->end(); IsParInHep_iter++ ) {
	  if (!*IsParInHep_iter) {
	    if (count==ParToSimPos) break;
	    else count++;
	  }
	}
	return IsParInHep_iter-IsParInHep->begin();
      }
      inline Int_t FindHepMCRecordingPosition( Int_t ParToHepPos ) {
	Int_t count=0;
	vector<Bool_t>::const_iterator IsParInHep_iter = IsParInHep->begin();
	for (; IsParInHep_iter != IsParInHep->end(); IsParInHep_iter++ )
	  if (*IsParInHep_iter) {
	    if (count==ParToHepPos) break;
	    else count++;
	  }
	return IsParInHep_iter-IsParInHep->begin();
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
   vector<Float_t> *pt,*eta,*phi,*Vertex_x,*Vertex_y,*Vertex_z,*isoR03sumPt,*isoR03emEt,*isoR03hadEt,*isoR03hoEt,*isoR03nJets,*isoR03nTracks,*isoR05sumPt,*isoR05emEt,*isoR05hadEt,*isoR05hoEt,*isoR05nJets,*isoR05nTracks,*isoemVetoEt,*isohadVetoEt,*isohoVetoEt,*normalizedChi2,*TrkRecoChi2,*CaloE_emMax,*CaloE_emS9,*CaloE_emS25,*CaloE_hadMax,*CaloE_hadS9,*Calo_emPos_R,*Calo_emPos_eta,*Calo_emPos_phi,*Calo_hadPos_R,*Calo_hadPos_eta,*Calo_hadPos_phi,*dEdx,*dEdxError,*TrkKink,*GlbKink;
   vector<Bool_t> *chargeMinus,*isGlobalMu,*isTrackerMu;
   vector<Int_t> *dEdx_numberOfSaturatedMeasurements,*dEdx_numberOfMeasurements;
   //Tracks
   edm::InputTag tracksTag;
   string dEdxTag;
   vector<UInt_t> *nValidTrackerHits,*nValidPixelHits,*nLostTrackerHits,*nLostPixelHits,*nBadMuonHits;
   UInt_t minTrackHits;
   //Muon Chamber Match
   vector<Float_t> *TrackDist[4],*TrackDistErr[4],*DXTrackToSegment[4],*DYTrackToSegment[4],*DXErrTrackToSegment[4],*DYErrTrackToSegment[4],*DRTrackToSegment[4],*DRErrTrackToSegment[4];
   vector<Bool_t> *IsSegmentBelongsToTrackByDR[4],*IsSegmentBelongsToTrackByCleaning[4];
   vector<UInt_t> *StationMask,*RequiredStationMask;
   double maxChamberDist,maxChamberDistPull;
   //PrimaryVertex
   string PrimaryVerticesTag;
   vector<Float_t> *vx,*vxError,*vy,*vyError,*vz,*vzError;
   //Tracking particles ( combination of simulated tracks )
   vector<Float_t> *TrkParticles_pt,*TrkParticles_eta,*TrkParticles_phi;
   vector<Double_t> *MatchQuality;
   vector<Int_t> *TrkParticles_pdgId,*TrkParticles_charge;
   //Muon Selectors
   unsigned int num_Cuts;
   vector<string> Cuts;
   vector<Bool_t> * SelectorPassed[30],*MySelector;
   //To Get the Decay Chains
      //saved variables
   Float_t HepMCFourVec[4];
   vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*Gen_vx,*Gen_vy,*Gen_vz,*Gen_vt;
   vector<Int_t> *Gen_pdgId;
   vector<Long64_t> *MuonType;//all possible chains; -1 is at the beginning of each chain; -2is at the beginning of each reco-muons;
   vector<Bool_t> *IsParInHep;
   vector<Int_t> *DChains,*theSameWithMuon;
      //temparory variables
   string HepMCTag;
   bool ChainRecord;
   vector< vector<Int_t> > SimChains;
   vector< vector<Int_t> > HepMCChains;
   vector<Int_t> DChain;
   vector<SimTrack *> ParToSim,MaskOut;
   vector<HepMC::GenParticle *> ParToHep;
   map< SimTrack *, vector<SimTrack *> > Daughters;
   vector<SimVertex> SVC;

   //Summarization
   struct Summary
   {
     ULong64_t Total_Events,Total_TrackerMuons,Total_GlobalMuon,Total_GlobalnotTrackerMuon;
   }Summarization;
   //Errors and Warnings reports
   struct ErrorMsg
   {
     Byte_t ErrorCode;
     ULong64_t Run,Event;
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
