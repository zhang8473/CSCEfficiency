#ifndef RPhysicsStudies_TPTrackMuonSys_H
#define RPhysicsStudies_TPTrackMuonSys_H

/** \class TPTrackMuonSys
 *  Author  : Chaouki Boulahouache, Jinzhong Zhang
 *  Intitute: Rice University, Northeastern University
 *  
 *  Package:  TPTrackMuonSys
 */

//#define GetCSCHitsBefore500
//#define GetPUFromEarlierThan_4_1_2
//#define GetPUFromEarlierThan_4_4_0
#define TrackingParticles
//#define m_debug
//#define jz_debug
// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <fstream>
#include <TCanvas.h>
#include <TH2F.h>
// user include files
// Framework...
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

/*Math*/
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Point3D.h"
/*Geometry*/
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCIndexer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"

/*Bad Chambers*/
#include "CalibMuon/CSCCalibration/interface/CSCConditions.h"
#include "CondFormats/CSCObjects/interface/CSCBadChambers.h"
#include "CondFormats/DataRecord/interface/CSCBadChambersRcd.h"

/*CSC-LCT*/
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "DataFormats/MuonDetId/interface/CSCTriggerNumbering.h"

//Global Triggers
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DataRecord/interface/L1MuTriggerScalesRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1MuTriggerScales.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Beamspot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// Tracks...
#include "DataFormats/TrackReco/interface/Track.h"  // reconstructed tracks that are stored in the AOD and RECO.
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"//isolation
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
//Track-detector associator
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

// Vertices...
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Muons...
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//MC truth information
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//Tracking Particles, Match to simulations/generation level particles
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

#include <vector>
#include <string.h>
#include <iostream>

#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

#define MAXNTRACKS 1500

//Chambers positions                   1/1,4      1/2         1/3       2/1,2    3/1,2      4/1,2
//const Float_t MEZOdd[6]      ={615.95, 711.76, 694.16, 839.96, 924.04, 1013.5}; // 839.96
//const Float_t MEZEven[6]     ={586.65, 684.36, 694.16, 815.16, 948.84, 1038.3}; // 839.96
const Float_t MEZ[6]         ={601.3,  696.11, 696.11, 827.56, 936.44, 1025.9}; // 839.96

class TPTrackMuonSys : public edm::EDAnalyzer {
 public:
  enum ParticleType {LightMeson=1,CharmedMeson=2,ccbarMeson=3,BottomMeson=4,bbarMeson=5,LightBaryon=6,CharmedBaryon=7,BottomBaryon=8,DiQuarks=9,Lepton=10,W=11,Z=12,Muon=13,JPsi=14,Other=15};
  enum TheTrackType {PromptMuFromW=10,PromptMuFromZ=11,PromptMuFromJPsi=12,PromptMuFromLightMeson=13,PromptMuFromHeavyMeson=14,PromptMuFromLightBaryon=15,PromptMuFromHeavyBaryon=16,NotPromptMufromWZ=17,PromptMuFromOthers=18,PunchThrough=20,PunchThroughAndDecayInFlight=21,DecayInFlightFromLightMeson=31,DecayInFlightFromHeavyMeson=32,DecayInFlightFromLightBaryon=33,DecayInFlightFromHeavyBaryon=34,NoMuSysHit=40,Others=01,NothingMatched=00};
  /// Constructor
  TPTrackMuonSys(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~TPTrackMuonSys();

 protected:
 private:
  virtual void beginJob() ;
  virtual void beginRun(const Run& r, const EventSetup& iSet);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  TrajectoryStateOnSurface cylExtrapTrkSam(reco::TrackRef track, double rho);
  TrajectoryStateOnSurface surfExtrapTrkSam(reco::TrackRef track, double z);
  FreeTrajectoryState freeTrajStateMuon(reco::TrackRef track);

  struct MCParticleInfo
  {
    Bool_t IsThisFromSimTrk,IsParticleFromGenerator,DoesParticleHaveMuonHit,IsParticleBornInsideOfARegion,IsPileup;
    Int_t pdgId;
  };
      
  MCParticleInfo MCParticleInfo_Creator( const SimTrack * thisTrk, TrackingParticleRef tpr );
  MCParticleInfo MCParticleInfo_Creator( HepMC::GenParticle * GenParticle ) {
    MCParticleInfo TBA;
    TBA.IsThisFromSimTrk=false;
    TBA.IsParticleFromGenerator=true;    TBA.DoesParticleHaveMuonHit=false;
    TBA.IsParticleBornInsideOfARegion=true;    TBA.IsPileup=false;
    TBA.pdgId=GenParticle->pdg_id();
    return TBA;
  }

  vector<SimVertex> SVC;
  vector< MCParticleInfo > MCParticlesList;
  vector<Int_t> DChain;
  vector< vector<Int_t> > SimChains;
  vector< vector<Int_t> > HepMCChains;
  vector<const SimTrack *> SavedSimTrk;
  vector<HepMC::GenParticle *> SavedHepPar;

  inline Int_t FindSimTrackInMCParticlesList( Int_t PosInSimTrackVec ) {
    Int_t count=0;
    vector<MCParticleInfo>::const_iterator MCParticlesList_iter = MCParticlesList.begin();
    for (; MCParticlesList_iter != MCParticlesList.end(); MCParticlesList_iter++ )
      if ( MCParticlesList_iter->IsThisFromSimTrk ) {
	if (count==PosInSimTrackVec) break;
	else count++;
      }
    return MCParticlesList_iter-MCParticlesList.begin();
  }

  inline Int_t FindHepMCInMCParticlesList( Int_t PosInHepMCVec ) {
    Int_t count=0;
    vector<MCParticleInfo>::const_iterator MCParticlesList_iter = MCParticlesList.begin();
    for (; MCParticlesList_iter != MCParticlesList.end(); MCParticlesList_iter++ )
      if ( !MCParticlesList_iter->IsThisFromSimTrk ) {
	if (count==PosInHepMCVec) break;
	else count++;
      }
    return MCParticlesList_iter-MCParticlesList.begin();
  }
      
  Bool_t GetDecayChains(TrackingParticleRef tpr, HepMC::GenEvent *HepGenEvent, ULong64_t &type, Int_t &truth_thesamewith, vector<vector< vector<Int_t> > > & SimChains);
  void SimTrackDaughtersTree(const SimTrack * thisTrk, TrackingParticleRef tpr);
  void HepMCParentTree(HepMC::GenParticle *genPar);
  inline ParticleType ParticleCata(Int_t pid);
  Bool_t SaveAndClassify(vector<Int_t> &Chain, vector<TheTrackType> &types, Int_t & truth_thesamewith, vector<vector< vector<Int_t> > > &  SimChains);
  TheTrackType Classify(vector<Int_t> &Chain);
  Bool_t IstheSameDChain(const vector<Int_t> &Chain1,const vector<Int_t> &Chain2);

  CSCConditions theDbConditions;
  HLTConfigProvider  hltConfigProvider_;
  
  void chamberCandidates(Int_t station, Float_t feta, Float_t phi, std::vector <int> &coupleOfChambers);
  Int_t ringCandidate(Int_t station, Float_t feta, Float_t phi);
  void fillChamberPosition();
  Short_t thisChamberCandidate(Short_t station, Short_t ring, Float_t phi);

  ///// Functions needed...
  void getCSCSegWkeyHalfStrip(const std::vector<CSCRecHit2D> &theseRecHits, Float_t &cStrp, Float_t &ckWG);
  
  Float_t YDistToHVDeadZone(Float_t yLocal, Int_t StationAndRing);

  vector<Float_t> GetEdgeAndDistToGap(reco::TrackRef trackRef, CSCDetId & detid);

  reco::MuonCollection::const_iterator matchTTwithMT(reco::TrackCollection::const_iterator &itrack);

  bool matchTTwithCSCRecHit(bool trackDir,
			    Int_t j, 
			    reco::TrackRef trackRef, 
			    edm::Handle<CSCRecHit2DCollection> recHits, 
			    //std::vector<CSCRecHit2DCollection> recHitOut, 
			    std::vector<CSCRecHit2D> &recHitOut, 
			    std::vector<Int_t > &deltaRecHitX,
			    std::vector<Int_t > &deltaRecHitY);
  

  bool matchTTwithRPCEChit(bool trackDir, 
			   Int_t j, 
			   reco::TrackRef trackRef, 
			   edm::Handle<RPCRecHitCollection> rpcRecHits, 
			   RPCRecHitCollection::const_iterator &rpcHitOut);

  inline Float_t TrajectoryDistToSeg( TrajectoryStateOnSurface *TrajSuf, CSCSegmentCollection::const_iterator segIt);
 
  TrajectoryStateOnSurface* matchTTwithCSCSeg( reco::TrackRef trackRef, edm::Handle<CSCSegmentCollection> cscSegments, 
					       CSCSegmentCollection::const_iterator &cscSegOut, CSCDetId & idCSC );

  Bool_t matchTTwithCSCSeg(Bool_t trackDir, Int_t j, reco::TrackRef trackRef, edm::Handle<CSCSegmentCollection> cscSegments, 
			   CSCSegmentCollection::const_iterator &cscSegOut );


  Bool_t matchCSCSegWithLCT(edm::Handle<CSCCorrelatedLCTDigiCollection> mpclcts, 
			    CSCDetId & idCSC, 
			    Int_t TT,
			    Float_t TrkPhi, Float_t TrkEta,
			    Float_t c1, Float_t w1,
			    CSCCorrelatedLCTDigiCollection::const_iterator &mpcItOut,
			    CSCCorrelatedLCTDigiCollection::const_iterator &mpcHsWkOut,
			    Bool_t *xMatch,
			    Float_t *mDAngle,
			    Float_t *diffTrkEta,
			    Float_t *diffTrkPhi,
			    Float_t *delHStrp,
			    Float_t *delWkey
			    );

  LocalPoint * matchTTwithLCTs(Float_t xPos, Float_t yPos, Short_t ec, Short_t st, Short_t &rg, Short_t cham, 
			       edm::Handle<CSCCorrelatedLCTDigiCollection> mpclcts, Float_t &dRTrkLCT, Int_t &lctBX );

  Int_t getNLayerMatchedCSCSeg(CSCSegmentCollection::const_iterator &cscSegMatch,
			       edm::Handle<CSCRecHit2DCollection> recHits,
			       Float_t *delRecSegX,
			       Float_t *delRecSegY,
			       Int_t &nGHits);

  CSCSectorReceiverLUT *srLUTs_[2][6][5];
  const L1MuTriggerScales* theTriggerScales;

  edm::Handle<CSCSegmentCollection> cscSegments;
  edm::Handle<reco::MuonCollection> muons;
  edm::ESHandle<MagneticField> theBField;
  //edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  // DT Geometry
  //edm::ESHandle<DTGeometry> dtGeom;
  edm::ESHandle<RPCGeometry> rpcGeo;
  // CSC Geometry
  edm::ESHandle<CSCGeometry> cscGeom;

  // Extrapolator to cylinder
  edm::ESHandle<Propagator> propagatorAlong;
  edm::ESHandle<Propagator> propagatorOpposite;
  edm::ParameterSet trackExtractorPSet_;
  // counters
  Int_t nEventsAnalyzed;
  Int_t treeCount;

  Float_t StationOnePhi[36];
  Float_t StationTwoOnePhi[18], StationTwoTwoPhi[36];
  Float_t StationThreeOnePhi[18], StationThreeTwoPhi[36];
  Float_t StationFourOnePhi[18], StationFourTwoPhi[36];


  unsigned long long etime; 
  Int_t bunchX, orbitNumb, expType, LumiBlock;
  //
  // The root file for the histograms.
  //
  TFile *theFile;

  //
  // General input parameters
  //
  std::string m_rootFileName, m_refRootFile, m_hltSingleMuTriName, m_hltDoubleMuTriName;
  double m_cellThreshold;
  // Steering parameters...
  edm::InputTag m_beamSpot;
  Bool_t m_isMC, m_doTrigger, m_doTrack, m_doMuon, m_doElectron, m_doEIDAndIso, m_doPhoton, m_doGenPart;
  Bool_t m_doTau, m_doJet, m_doBJet, m_doMET, m_doCaloTiming, m_doIslandHybrid;
  edm::InputTag m_gTracksTag, m_fitTracksTag, m_electronTag, m_muonTag, m_photonTag, m_mcTag;
  edm::InputTag m_jet0Tag, m_jet1Tag, m_jet2Tag, m_jet3Tag, m_jet4Tag, m_jet5Tag, m_jet6Tag;
  edm::InputTag m_trujet0Tag, m_trujet1Tag, m_trujet2Tag, m_trujet3Tag, m_trujet4Tag, m_trujet5Tag, m_trujet6Tag, m_vertexSrc, m_dEdxDiscrimTag;
  edm::InputTag m_HepMCTag;
  edm::InputTag trackProducer;
  edm::InputTag m_hlt, m_hltTrgEv, m_L1extraTag;
  vector<Int_t> m_HLTMuTrgBit;
  Int_t m_HLTDiMuTrgBit;
  vector<string> m_HLTMuTrgNames;
  string m_HLTDiMuTrgName;

  //badchambers
  edm::ESHandle<CSCBadChambers> pBad;
  CSCBadChambers* badChambers_;

  /// Ntuple variables
  TTree *fractNtuple, *RunInfo;
  //HLT
  vector<string> *HLTMuNames,*HLTMuObjModuleNames;
  string *HLTTableName, *HLTDiMuName, *HLTDiMuObjModuleName;
  vector<Int_t> *badChambersIndices; 
  vector<Bool_t> *HLTMuAcceptance;
  Bool_t HLTDiMuAcceptance,trgSingle;
  vector<Float_t> *minDRHLTMu;
  Float_t minDRHLTAllSingleMu,minDRHLTDiMu;

  // for general information
  Int_t isItMC, run_number, event_number;
  Float_t mcweight;
  UInt_t numberOfPUVertices; // the number of pileup interactions that have been added to the event from BX=0
  Float_t numberOfPUVerticesMixingTruth;// the "MixingTruth" mean number of pileup interactions for this event from which each bunch crossing has been sampled; same for all bunch crossings in an event (before possion smearing); in fact BX=-1, this value is zero, it's a bug I believe.
  UInt_t numberOfPUVerticesTot;  // all PX : BX=0 and out-of-time PU
  UInt_t numberOfPrimaryVertices ;
  //// for track information
  //
  Float_t tracks_ecalCrossedE, tracks_ecal3E, tracks_ecal5E, tracks_ecal3EMax, tracks_ecal5EMax;
  Float_t tracks_ecalTrueE, tracks_trkPosAtEcalEta, tracks_trkPosAtEcalPhi, tracks_ecalMaxPosEta, tracks_ecalMaxPosPhi;
  Float_t tracks_hcalCrossedE, tracks_hcal3E, tracks_hcal5E, tracks_hcal3EMax, tracks_hcal5EMax;
  Float_t tracks_hcalTrueE, tracks_trkPosAtHcalEta, tracks_trkPosAtHcalPhi, tracks_hcalMaxPosEta, tracks_hcalMaxPosPhi;
  Float_t tracks_hcalTrueECorrected;
  Float_t tracktruth_pt, tracktruth_e, tracktruth_p, tracktruth_id;
  Bool_t tracktruth_isPileup;
  ULong64_t tracktruth_type;
  Int_t tracktruth_thesamewith;
  TrackAssociatorParameters parameters_;
  TrackDetectorAssociator trackAssociator_;

  // Tracks
  Int_t tracks_algo, tracks_charge, tracks_numberOfValidHits, tracks_algN, tracks_numberOfLostHits, tracks_recHitsSize;
  Float_t tracks_id, tracks_chi2, tracks_e;
  Float_t tracks_px, tracks_py, tracks_pz, tracks_pt;
  Float_t tracks_eta, tracks_phi, tracks_dxy, tracks_d0, tracks_dsz, tracks_dz;
  Float_t tracks_vx, tracks_vy, tracks_vz, tracks_ndof;
  Float_t tracks_qoverp, tracks_lambda, tracks_quality;
  Float_t tracks_IsoR03Ratio,tracks_IsoR05Ratio;//ratio of track pT
  Float_t tracks_qoverpError, tracks_ptError, tracks_thetaError, tracks_lambdaError;
  Float_t tracks_etaError, tracks_phiError, tracks_dxyError, tracks_d0Error , tracks_dszError;
  Float_t tracks_dzError;
  Bool_t tracks_isCaloMuTrk,tracks_isTrackerMuTrk,trackVeto;

  Int_t nPosTrk, nNegTrk, nTotalTrks;
  ///
  Int_t myRegion;// 1-DT, 2-DT&CSC 3-CSC;

  Int_t mpc_endcap, mpc_ring, mpc_station, mpc_chamber;
  Int_t mpc_triggerSector, mpc_triggerCscId, mpc_trknmb, mpc_quality, mpc_bX, mpc_bX0, mpc_strip;
  Int_t mpc_pattern, mpc_bend, mpc_keyWG, mpc_syncErr, mpc_cscID, mpc_mpcLink;
  Float_t trkPhiScale, trkEtaScale, minDAngle, CSCSegMPCMatch, minPhi;
  
  Float_t muon_px, muon_py, muon_pz, muon_pt, muon_e, muon_eta, muon_phi, muon_charge, dRMU;
  Float_t muon_nChambers, muon_nChambersMSeg;
  Int_t muon_id, motherID;

  Float_t centerStrp, centerkWG;
  Int_t lctMatch1, lctMatch2, TTcase;
  Float_t lctPhiDiff1, lctEtaDiff1, lctDelHStrp1, lctDelWkey1; 
  Float_t lctPhiDiff2, lctEtaDiff2, lctDelHStrp2, lctDelWkey2; 

  Float_t MuTagPx, MuTagPy, MuTagPz, MuProbePx, MuProbePy, MuProbePz;
  Int_t MuTagHitsMuSys, MuTagHitsTrkSys;
  Float_t MuTagE, MuProbeE;
  Float_t MuTagPhiProj1, MuTagEtaProj1, MuTagPhiProj2, MuTagEtaProj2;
  Float_t MuTagIsoR03Ratio,MuTagIsoR05Ratio;//ratio of muon pT
  Float_t MuTagtracktruth_pt, MuTagtracktruth_p, MuTagtracktruth_id;
  Bool_t MuTagtracktruth_isPileup;
  ULong64_t MuTagtracktruth_type;
  Int_t MuTagtracktruth_thesamewith;

  /// For the ntuples..
  //
  Bool_t MuTagCaloL, MuTagCaloT, iSameVtx;

  Int_t mu_found, MuProbenHitsMuSys, MuProbenHitsTrkSys, MuProbenHitsPixSys;
  Int_t MuTagPromt, MuTagnSegTrkArb, MuProbeCharge;
  Float_t MuProbePt, MuProbeEta, MuProbePhi, MuTagPt, MuTagEta, MuTagPhi, invMass, deltaRTrkMu;      
  Float_t vtx_r, vtx_z, vtx_rError, vtx_zError, vtx_normChi2;
  Int_t vtx_size, igSameVtx;

  Int_t nTrkCountCSCSeg;
  /*CSC Chamber Candidates in each station*/
  Bool_t CSCEndCapPlus;
  Short_t CSCRg[4],CSCChCand[4];
  Bool_t CSCChBad[4];
  
  /*Extrapolated Tracks on CSC Chamber Candidates in each station*/
  Float_t CSCDyProjHVGap[4],CSCDyErrProjHVGap[4],CSCProjEdgeDist[4],CSCProjEdgeDistErr[4];//note: there is no so-called extrapolated track position on CSC chambers, because you need to know which layer.
  
  /*Segments characteristics*/
  Float_t CSCSegxLc[4],CSCSegyLc[4],CSCSegxErrLc[4],CSCSegyErrLc[4],CSCSegChisqProb[4],CSCdXdZTTSeg[4],CSCdYdZTTSeg[4];
  Int_t CSCnSegHits[4];

  /*Distance from the Extrapolated Tracks to CSC Segments, 99999. for no CSC segment found*/
  Float_t CSCDxTTSeg[4],CSCDxErrTTSeg[4],CSCDyTTSeg[4],CSCDyErrTTSeg[4],CSCDxyTTSeg[4],CSCDxyErrTTSeg[4];

  /*LCT characteristics*/
  Float_t CSCLCTxLc[4],CSCLCTyLc[4];
  Int_t CSCLCTbx[4];

  /*Distance from the Extrapolated Tracks to LCT, 99999. for no LCT found*/
  Float_t CSCDxTTLCT[4],CSCDxErrTTLCT[4],CSCDyTTLCT[4],CSCDyErrTTLCT[4],CSCDxyTTLCT[4],CSCDxyErrTTLCT[4];

  /*DetlaR between the extrapolated tracker track on muon system and the tagged muon*/
  Float_t dRTkMu[4];  
  /*Default decision of whether a segment or LCT is found*/
  Int_t segSt[4],lctSt[4];
};
#endif
