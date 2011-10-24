#define TrackingParticles
//#define DimuonAnalysis
// system include files 
#include <memory>
#include <vector>
#include <map>
#include <cmath>
#include <cstring>
#include <algorithm>

// Framework Core
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Standard Filters
#include "TriggerResultsFilter.h"
#include "GoodVertexFilter.h"
#include "FilterOutScraping.h"

//Simulation Infomation
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Muon
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

//Tracks
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
//Tracking Tools
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//Tracking Particles, Match to simulated tracks
#ifdef TrackingParticles
  #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
  #include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
  #include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
  #include "SimTracker/Records/interface/TrackAssociatorRecord.h"
  #include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
  #include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
  //SimTracks
  #include "SimDataFormats/Track/interface/SimTrackContainer.h"
  //SimVertex
  #include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#else
  #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
  #include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#endif
//HLT
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//Vertex
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

//CSC Digis
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

//Geometry
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"

//MSG Logger
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//ROOT
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

#define maxFilterObjects 30
#define MaxHLTObjDeviation2 0.1 // the maximum deviation between pt and eta of reco objs and HLT objs to be accepted as a HLT objs. this value is the veto cone of muon squared
using namespace edm;
using namespace std;
using namespace reco;

class MuonSysPriEff : public EDAnalyzer {
 public:
  explicit MuonSysPriEff(const edm::ParameterSet&);
  ~MuonSysPriEff();
  enum ParticleType {LightMeson=1,CharmedMeson=2,ccbarMeson=3,BottomMeson=4,bbarMeson=5,LightBaryon=6,CharmedBaryon=7,BottomBaryon=8,DiQuarks=9,Lepton=10,W=11,Z=12,Other=13};
  enum TheMuonType {PromptMuFromWZ=10,PromptMuFromLightMeson=11,PromptMuFromHeavyMeson=12,PromptMuFromLightBaryon=13,PromptMuFromHeavyBaryon=14,NotPromptMufromWZ=15,PromptMuFromOthers=16,PunchThrough=20,PunchThroughAndDecayInFlight=21,DecayInFlightFromLightMeson=31,DecayInFlightFromHeavyMeson=32,DecayInFlightFromLightBaryon=33,DecayInFlightFromHeavyBaryon=34,NoMuSysHit=40,Others=1,NothingMatched=0};
 private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual reco::Muon::ArbitrationType MuonArbitrationTypeFromString( const std::string &);
  inline  void ClearVecs_RECO();
  inline  Float_t YDistToHVDeadZone(Float_t yLocal, Int_t StationRing);
#ifdef TrackingParticles
  virtual void GetDecayChains(TrackingParticleRef tpr,vector<int> *DChains, vector <TheMuonType> &type, HepMC::GenEvent *HepGenEvent);
  virtual void HepMCParentTree(HepMC::GenParticle *);
  virtual void SimTrackDaughtersTree(const SimTrack * thisTrk, TrackingParticleRef tpr);
  inline Bool_t IstheSameDChain(vector<int> &,vector<int> &);
  inline Bool_t SimHitsToSegment(const PSimHit& hit);
  inline ParticleType ParticleCata(int);
  inline TheMuonType classification(vector<Int_t> &);
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
  inline Int_t FindHepMCRecordingPosition( Int_t SavedHepParPos ) {
    Int_t count=0;
    vector<Bool_t>::const_iterator IsParInHep_iter = IsParInHep->begin();
    for (; IsParInHep_iter != IsParInHep->end(); IsParInHep_iter++ )
      if (*IsParInHep_iter) {
	if (count==SavedHepParPos) break;
	else count++;
      }
    return IsParInHep_iter-IsParInHep->begin();
  }
#endif
  //-------------control parameters------------
  Bool_t FirstEvent,IsMC,HasDigi,ThrowBadEvents;
  TriggerResultsFilter *HLTFilter;
  GoodVertexFilter *VertexFilter;
  FilterOutScraping *ScrapingFilter;
  //-------------member data------------
  //TTree&TFile
  TFile *file;
  string FileName;
  TTree *Muons_Tree,*RunInfo_Tree;
  //Event Information
  struct General_EvtInfo
  {
    ULong64_t RUN,EVENT,LS,ORBIT,BX;
  } Info;
  Bool_t isRealData,isHLTTriggerred,isGoodVertex,isNoScrapping;
  Double_t GenEventWeight;
  //PU information
  Int_t numberOfPUVertices;
  //Muon
  Float_t MuonPtCut,MuonEtaMin,MuonEtaMax,MinNumberOfMuons,DiMuonInvarMassCut;
  vector<Float_t> *pt,*eta,*phi,*Vertex_x,*Vertex_y,*Vertex_z,*isoR03sumPt,*isoR05sumPt,*normalizedChi2,*TrkRelChi2,*dEdx,*dEdxError,*TrkKink,*GlbKink;
  vector<Bool_t> *chargeMinus,*isGlobalMu,*isTrackerMu,*IsPUTrack;
  vector<Int_t> *dEdx_numberOfSaturatedMeasurements,*dEdx_numberOfMeasurements,*numberOfMatchedSegments,*numberOfMatchedStations;
#ifdef DimuonAnalysis
  //Dimuon saved position: num_mu*mu1-mu1*(mu1+1)/2+(mu2-mu1)-1 (mu2>mu1) (mu1 starts at 0)
  vector<Float_t> *DiMuonInvariantMass,*CosThetaStar,*AngleBetweenDiMuon;
#endif
  //Tracks
  edm::InputTag tracksTag;
  string dEdxTag;
  vector<UInt_t> *InnerTrack_nValidTrackerHits,*InnerTrack_nValidPixelHits,*InnerTrack_nLostTrackerHits,*InnerTrack_nLostPixelHits,*InnerTrack_ndof,*GlobalTrack_ndof,*numberOfValidMuonHits;
  vector<Float_t> *InnerTrack_chi2,*GlobalTrack_chi2,*DXYwtBS,*DXYwtPV,*DXYErrwtPV;
  //CSC Digis
  edm::InputTag CSCDigisTag;
  //Muon Selectors
  unsigned int num_Cuts;
  reco::Muon::ArbitrationType MuonArbitrationType;
  vector<muon::SelectionType> OfficialMuonSelectors;
  vector<Bool_t> * SelectorPassed[30];
  //Muon Chamber and Segment Match
  vector<Float_t> *TrackDistToChamberEdge,*TrackDistToChamberEdgeErr,*XTrack,*YTrack,*XErrTrack,*YErrTrack,*XSegment,*YSegment,*XLCT,*YLCT,*XErrSegment,*YErrSegment,*dYToDeadZone;
  vector<Bool_t> *IsCSCChamber,*IsSegmentOwnedExclusively,*IsSegmentBelongsToTrackByDR,*IsSegmentBelongsToTrackByCleaning;
  vector<Byte_t> *SectorChamber,*MuonIndex,*NumberOfLCTsInChamber,*NumberOfHitsInSegment;
  vector<Char_t> *StationRing;
  double maxChamberDist,maxChamberDistPull;
  //PrimaryVertex
  string PrimaryVerticesTag;
  vector<Float_t> *vx,*vxError,*vy,*vyError,*vz,*vzError;
  //Run Information and HLT
  UInt_t HLTSize;
  HLTConfigProvider HLTConfig;
  string *HLTTableName;
  string HLTProc;
  vector<string> *HLTNamesSet;
  vector<Bool_t> *HLTacceptance;
  edm::InputTag TriggerResultsTag;
  ULong64_t NumberOfEvents;//total number of events in this run

  //HLTObjects
  UInt_t num_HLTsSaveObjs;
  Int_t HLTFilterNamesAcceptenceIndex[maxFilterObjects];
  vector<Float_t> *HLTObj_pt[maxFilterObjects],*HLTObj_eta[maxFilterObjects],*HLTObj_phi[maxFilterObjects];
  vector<Bool_t> *isHLTObj[maxFilterObjects];
  vector<string> HLTObj_HLTNames;
  vector<edm::InputTag> HLT_ModuleNames;
  edm::InputTag TriggerEventTag;

#ifdef TrackingParticles // combination of simulated tracks
  edm::InputTag TPInputTag;
  UInt_t minTrackHits;
  vector<Float_t> *TrkParticles_pt,*TrkParticles_eta,*TrkParticles_phi;
  vector<Int_t> *TrkParticles_pdgId,*TrkParticles_charge;
  vector<Double_t> *SharedHitsRatio,*MCMatchChi2;
  vector<Bool_t> *IsParInHep,*IsParHasMuonHits;
  vector<Int_t> *TTTruthDChains,*SegTruthDChains,*theSameWithMuon;//all possible chains; -1 is at the beginning of each chain; -2is at the beginning of each reco-muons;
  vector<Long64_t> *TTTruthMuType,*SegTruthMuType;
  //temparory variables
  bool ChainRecord;
  vector< vector<Int_t> > SimChains;
  vector< vector<Int_t> > HepMCChains;
  vector<Int_t> DChain;
  vector<const SimTrack *> SavedSimTrk;
  vector<HepMC::GenParticle *> SavedHepPar;
  vector<TrackingParticleRef> SavedTP;
  map< SimTrack *, vector<SimTrack *> > Daughters;
  vector<SimVertex> SVC;
  //SimHits
  struct CSCChamberSimHitsInfo {
     Int_t ChamberID;
     Float_t FirstHit_X,FirstHit_Y,LastHit_X,LastHit_Y;
     Byte_t HitsMask;//1-6 bit stand for 6 layers
     CSCChamberSimHitsInfo(Int_t ChamberID,Float_t a,Float_t b,Float_t c,Float_t d,Byte_t e):ChamberID(ChamberID),FirstHit_X(a),FirstHit_Y(b),LastHit_X(c),LastHit_Y(d),HitsMask(e) {}//creator
   };
   vector<CSCChamberSimHitsInfo> ChamberSimHits;//the first and last layer SimHits in one chamber. Int_t is the ChamberID
  //Others
  ULong64_t NumMisMatch;
#endif

  Float_t HepMCFourVec[4];
  vector<Float_t> *Gen_pt,*Gen_eta,*Gen_phi,*Gen_vx,*Gen_vy,*Gen_vz,*Gen_vt;
  vector<Int_t> *Gen_pdgId;
  string HepMCTag;
  inline void ClearVecs_HepMC();

  inline Bool_t IsTheSameSegment(UInt_t Seg1,UInt_t Seg2) {
    if ( (*IsCSCChamber)[Seg1]!=(*IsCSCChamber)[Seg2] ) return false;
    if ( (*XSegment)[Seg1]!=(*XSegment)[Seg2] || (*YSegment)[Seg1]!=(*YSegment)[Seg2] ) return false;
    if ( (*StationRing)[Seg1] != (*StationRing)[Seg2] || (*SectorChamber)[Seg1] != (*SectorChamber)[Seg2] ) return false;
    return true;
  }
};
