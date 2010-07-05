// system include files 
#include <memory>
#include <vector>
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

//Tracking Particles, Match to simulated tracks
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

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

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void HepMCDaughtersTree(HepMC::GenParticle *genPar);
      virtual void HepMCParentTree(HepMC::GenParticle *);
      virtual void SimTrackDaughtersTree(edm::SimVertexContainer &, edm::SimTrackContainer &, vector< vector<int> > &, unsigned int);
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
   Double_t minP,maxP,MaxRelpT,MaxDR,tmp_MaxDR2,MaxDR2,maxNormChi2;
   UInt_t minTrackHits;
   //Tracker Muon
   vector<Float_t> *pt,*eta,*phi,*Vertex_x,*Vertex_y,*Vertex_z,*isoR03sumPt,*isoR03emEt,*isoR03hadEt,*isoR03hoEt,*isoR03nJets,*isoR03nTracks,*isoR05sumPt,*isoR05emEt,*isoR05hadEt,*isoR05hoEt,*isoR05nJets,*isoR05nTracks,*isoemVetoEt,*isohadVetoEt,*isohoVetoEt,*normalizedChi2,*TrkRecoChi2,*CaloE_emMax,*CaloE_emS9,*CaloE_emS25,*CaloE_hadMax,*CaloE_hadS9,*Calo_emPos_R,*Calo_emPos_eta,*Calo_emPos_phi,*Calo_hadPos_R,*Calo_hadPos_eta,*Calo_hadPos_phi,*dEdx,*dEdxError,*TrkKink,*GlbKink;
   vector<Bool_t> *chargeMinus,*isGlobalMu;
   vector<Int_t> *dEdx_numberOfSaturatedMeasurements,*dEdx_numberOfMeasurements,*whichGenPar,*whichGenId;
   //Simulated Particles
   vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*Gen_vx,*Gen_vy,*Gen_vz,*Gen_vt;
   vector<Int_t> *Gen_pdgId;
   vector<Byte_t> *Gen_charge;
   vector<Int_t> DChain, *DChains; //Simulated chains; -1 is at the beginning of each chain; -2is at the beginning of each reco-muons;
   int particlenum;
   //PrimaryVertex
   string PrimaryVerticesTag;
   vector<Float_t> *vx,*vxError,*vy,*vyError,*vz,*vzError;
   //Tracking particles ( combination of simulated tracks )
   vector<Float_t> *TrkParticles_pt,*TrkParticles_eta,*TrkParticles_phi,*TrkParticles_chi2;
   vector<Int_t> *TrkParticles_pdgId;
   string MatchAlgorithm;
   //Tracks
   string tracksTag,dEdxTag;
   //HepMC Collection
   string HepMCTag;
   //Summarization
   struct Summary
   {
     ULong64_t Total_Events,Total_RecoMuons,Total_GenMuons,Total_GenElectrons,Total_GenPhotons;
   }Summarization;
   //Errors and Warnings reports
   struct ErrorMsg
   {
     Byte_t ErrorCode;
     ULong64_t Run,Event;
   }Error;
   //Others
   static Bool_t FirstEntry;
   Bool_t LocalRun;
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
