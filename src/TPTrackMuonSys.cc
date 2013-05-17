// -*- C++ -*-
//
// Package:    TPTrackMuonSys
// Class:      TPTrackMuonSys
// 
/**\class TPTrackMuonSys TPTrackMuonSys.cc RPhysicsStudies/TPTrackMuonSys/src/TPTrackMuonSys.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Chaouki Boulahouache,32 4-A05,+41227674763;
//           Jinzhong Zhang,32 4-C21,+41227671337,
//         Created:  Wed May 12 15:52:13 CEST 2010
// $Id$
//
// system include files
#include <memory>

// user include files
#include "../interface/TPTrackMuonSys.h"

TPTrackMuonSys::TPTrackMuonSys(const edm::ParameterSet& Conf) : theDbConditions( Conf ) {
  //now do what ever initialization is needed
  m_rootFileName         = Conf.getUntrackedParameter<std::string>("rootFileName","PhyTree.root");
  ///////////////////////////////
  // Various input parameters.
  //////////////////////////////

  // A debug flag, here I can set different debugging steps based on the m_debug value,
  m_gTracksTag = Conf.getUntrackedParameter<edm::InputTag>("gTracksTag"); // "generalTracks"
  m_dEdxDiscrimTag = Conf.getUntrackedParameter<edm::InputTag>("dedxTag"); //dedxHarmonic2
  trackProducer     = Conf.getParameter<edm::InputTag>("trackProducer");
  m_vertexSrc = Conf.getParameter<edm::InputTag>("vertexCollection");

  m_hlt               = Conf.getUntrackedParameter<edm::InputTag>("hltTag"); //TriggerResults::HLT
  m_hltTrgEv          = Conf.getUntrackedParameter<edm::InputTag>("hltEvTag"); //TriggerResults::HLT
  m_L1extraTag        = Conf.getUntrackedParameter<edm::InputTag>("L1extraTag");

  m_HepMCTag = Conf.getUntrackedParameter<string>("HepMCTag","generator");

  m_HLTMuTrgNames = Conf.getParameter< vector<string> >("HLTMuTrgNames");
  m_HLTDiMuTrgName = Conf.getParameter< string >("HLTDiMuTrgName");

  // TrackAssociator parameters
  edm::ParameterSet parameters = Conf.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  trackExtractorPSet_ = Conf.getParameter<edm::ParameterSet>("TrackExtractor");
  parameters_.loadParameters( parameters );
  trackAssociator_.useDefaultPropagator();

  ///  Now the MC specific information... if available...
  //
  m_isMC     = Conf.getUntrackedParameter<bool>("isMC");
  m_mcTag = Conf.getUntrackedParameter<edm::InputTag>("mcTag"); //
 
  // flags to switch on/off individual modules 
  // set counters to zero
  nEventsAnalyzed = 0;
  treeCount = 0;
  
  // Create the root file for the histograms
  theFile = new TFile(m_rootFileName.c_str(), "RECREATE");


  //////////////////// for LUTs
  /*
  bzero(srLUTs_, sizeof(srLUTs_));
  //Int_t endcap =1, sector =1;
  Bool_t TMB07=true;
  edm::ParameterSet srLUTset;
  srLUTset.addUntrackedParameter<bool>("ReadLUTs", false);
  srLUTset.addUntrackedParameter<bool>("Binary",   false);
  srLUTset.addUntrackedParameter<std::string>("LUTPath", "./");
  for(Int_t endcapItr = CSCDetId::minEndcapId(); endcapItr <= CSCDetId::maxEndcapId(); endcapItr++){
    for(Int_t sectorItr = CSCTriggerNumbering::minTriggerSectorId();sectorItr <= CSCTriggerNumbering::maxTriggerSectorId();sectorItr++){
      for(Int_t stationItr = 1; stationItr <= 4; stationItr++){
	if(stationItr == 1){
	  for(Int_t subsectorItr = 0; subsectorItr < 2; subsectorItr++){
	    srLUTs_[endcapItr-1][sectorItr-1][subsectorItr] = new CSCSectorReceiverLUT(endcapItr, sectorItr, subsectorItr+1, stationItr, srLUTset, TMB07); 
	  }
	} else {
	  srLUTs_[endcapItr-1][sectorItr-1][stationItr] = new CSCSectorReceiverLUT(endcapItr, sectorItr, 0, stationItr, srLUTset, TMB07); 
	} //if for station 1 or 234
      } // stationItr loop
    } // sectorItr loop
  } // endcapItr loop
  */
  // ------------------------------
  // --- Summary ntuple booking ---
  // ------------------------------

  fractNtuple         = new TTree("Fraction","/Fraction");
  fractNtuple->Branch("run_number",      &run_number, "run_number/I") ;
  fractNtuple->Branch("event_number",    &event_number, "event_number/I") ;
  fractNtuple->Branch("LumiBlock",       &LumiBlock,"LumiBlock/I") ;
  fractNtuple->Branch("bunchX",          &bunchX, "bunchX/I") ;
  fractNtuple->Branch("orbitNumb",       &orbitNumb,"orbitNumb/I") ;
  fractNtuple->Branch("mcweight",       &mcweight,"mcweight/F") ;
  fractNtuple->Branch("numberOfPUVertices",       &numberOfPUVertices, "numberOfPUVertices/i");
  fractNtuple->Branch("numberOfPUVerticesMixingTruth",       &numberOfPUVerticesMixingTruth, "numberOfPUVerticesMixingTruth/F");
  fractNtuple->Branch("numberOfPUVerticesTot",       &numberOfPUVerticesTot, "numberOfPUVerticesTot/i");
  fractNtuple->Branch("numberOfPrimaryVertices" , &numberOfPrimaryVertices , "numberOfPrimaryVertices/i"); 
  fractNtuple->Branch("trgSingle",       &trgSingle,"trgSingle/O") ;
  fractNtuple->Branch("nTrkCountCSCSeg",    &nTrkCountCSCSeg,   "nTrkCountCSCSeg/I") ;
  //  fractNtuple->Branch("nTrkCountRPCE",    &nTrkCountRPCE,   "nTrkCountRPCE/I") ;
  //fractNtuple->Branch("nTrkCountEC",       &nTrkCountEC,      "nTrkCountEC/I") ;
  fractNtuple->Branch("nPosTrk",       &nPosTrk,      "nPosTrk/I") ;
  fractNtuple->Branch("nNegTrk",       &nNegTrk,   "nNegTrk/I") ;
  fractNtuple->Branch("nTotalTrks",       &nTotalTrks,      "nTotalTrks/I") ;
  fractNtuple->Branch("trackVeto",       &trackVeto,      "trackVeto/O") ;
  fractNtuple->Branch("myRegion",           &myRegion,           "myRegion/I") ;
  fractNtuple->Branch("MuTagPt",         &MuTagPt,         "MuTagPt/F") ;
  fractNtuple->Branch("MuTagEta",        &MuTagEta,        "MuTagEta/F") ;
  fractNtuple->Branch("MuTagPhi",        &MuTagPhi,        "MuTagPhi/F") ;
  fractNtuple->Branch("MuTagIsoR03Ratio",        &MuTagIsoR03Ratio,        "MuTagIsoR03Ratio/F") ;
  fractNtuple->Branch("MuTagIsoR05Ratio",        &MuTagIsoR05Ratio,        "MuTagIsoR05Ratio/F") ;
  fractNtuple->Branch("MuTagPromt",      &MuTagPromt,      "MuTagPromt/I");
  fractNtuple->Branch("MuTagnSegTrkArb", &MuTagnSegTrkArb, "MuTagnSegTrkArb/I");
  fractNtuple->Branch("MuTagCaloL",      &MuTagCaloL,      "MuTagCaloL/O") ;
  fractNtuple->Branch("MuTagCaloT",      &MuTagCaloT,      "MuTagCaloT/O") ;
  fractNtuple->Branch("MuTagtracktruth_pt",            &MuTagtracktruth_pt,           "MuTagtracktruth_pt/F") ;
  fractNtuple->Branch("MuTagtracktruth_p",            &MuTagtracktruth_p,           "MuTagtracktruth_p/F") ;
  fractNtuple->Branch("MuTagtracktruth_id",            &MuTagtracktruth_id,           "MuTagtracktruth_id/F") ;
  fractNtuple->Branch("MuTagtracktruth_type",            &MuTagtracktruth_type,           "MuTagtracktruth_type/l") ;
  fractNtuple->Branch("MuTagtracktruth_isPileup",            &MuTagtracktruth_isPileup,         "MuTagtracktruth_isPileup/O") ;
  fractNtuple->Branch("MuTagtracktruth_thesamewith",            &MuTagtracktruth_thesamewith,           "MuTagtracktruth_thesamewith/I") ;
  fractNtuple->Branch("vtx_r",       &vtx_r,      "vtx_r/F") ;
  fractNtuple->Branch("vtx_z",       &vtx_z,      "vtx_z/F") ;
  fractNtuple->Branch("vtx_rError",       &vtx_rError,      "vtx_rError/F") ;
  fractNtuple->Branch("vtx_zError",       &vtx_zError,      "vtx_zError/F") ;
  fractNtuple->Branch("vtx_normChi2",       &vtx_normChi2,      "vtx_normChi2/F") ;
  fractNtuple->Branch("vtx_size",       &vtx_size,      "vtx_size/I") ;
  fractNtuple->Branch("iSameVtx",       &iSameVtx,      "iSameVtx/O") ;
  fractNtuple->Branch("invMass",            &invMass,            "invMass/F");
  fractNtuple->Branch("tracks_e",         &tracks_e,         "tracks_e/F");
  fractNtuple->Branch("tracks_pt",        &tracks_pt,        "tracks_pt/F");
  fractNtuple->Branch("tracks_eta",       &tracks_eta,       "tracks_eta/F");
  fractNtuple->Branch("tracks_phi",       &tracks_phi,       "tracks_phi/F");
  fractNtuple->Branch("tracks_charge",    &tracks_charge,    "tracks_charge/I");
  fractNtuple->Branch("tracks_id",        &tracks_id,        "tracks_id/I");
  fractNtuple->Branch("tracks_normchi2",      &tracks_chi2,      "tracks_normchi2/F");
  fractNtuple->Branch("tracks_dxy",       &tracks_dxy,       "tracks_dxy/F");
  fractNtuple->Branch("tracks_dz",        &tracks_dz,        "tracks_dz/F");
  fractNtuple->Branch("tracks_vx",        &tracks_vx,        "tracks_vx/F");
  fractNtuple->Branch("tracks_vy",        &tracks_vy,        "tracks_vy/F");
  fractNtuple->Branch("tracks_vz",        &tracks_vz,        "tracks_vz/F");
  fractNtuple->Branch("tracks_IsoR03Ratio",        &tracks_IsoR03Ratio,        "tracks_IsoR03Ratio/F");
  fractNtuple->Branch("tracks_IsoR05Ratio",        &tracks_IsoR05Ratio,        "tracks_IsoR05Ratio/F");
  fractNtuple->Branch("tracks_qoverp",    &tracks_qoverp,    "tracks_qoverp/F");
  fractNtuple->Branch("tracks_lambda",    &tracks_lambda,    "tracks_lambda/F");
  fractNtuple->Branch("tracks_recHitsSize", &tracks_recHitsSize,   "tracks_recHitsSize/I");
  fractNtuple->Branch("tracks_numberOfValidHits",        &tracks_numberOfValidHits,        "tracks_numberOfValidHits/I");
  fractNtuple->Branch("tracks_numberOfLostHits",        &tracks_numberOfLostHits,        "tracks_numberOfLostHits/I");
  fractNtuple->Branch("tracks_qoverpError",        &tracks_qoverpError,        "tracks_qoverpError/F");
  fractNtuple->Branch("tracks_ptError",        &tracks_ptError,        "tracks_ptError/F");
  fractNtuple->Branch("tracks_thetaError",        &tracks_thetaError,        "tracks_thetaError/F");
  fractNtuple->Branch("tracks_lambdaError",        &tracks_lambdaError,        "tracks_lambdaError/F");
  fractNtuple->Branch("tracks_etaError",        &tracks_etaError,        "tracks_etaError/F");
  fractNtuple->Branch("tracks_phiError",        &tracks_phiError,        "tracks_phiError/F");
  fractNtuple->Branch("tracks_dxyError",        &tracks_dxyError,        "tracks_dxyError/F");
  fractNtuple->Branch("tracks_d0Error",        &tracks_d0Error,        "tracks_d0Error/F");
  fractNtuple->Branch("tracks_dszError",        &tracks_dszError,        "tracks_dszError/F");
  fractNtuple->Branch("tracks_dzError",        &tracks_dzError,        "tracks_dzError/F");
  fractNtuple->Branch("tracks_isCaloMuTrk",        &tracks_isCaloMuTrk,        "tracks_isCaloMuTrk/O");
  fractNtuple->Branch("tracks_isTrackerMuTrk",        &tracks_isTrackerMuTrk,        "tracks_isTrackerMuTrk/O");
  // fractNtuple->Branch("tracks_dedx",            &tracks_dedx,           "tracks_dedx/F") ;
  // fractNtuple->Branch("tracks_dedxErr",            &tracks_dedxErr,           "tracks_dedxErr/F") ;
  //  fractNtuple->Branch("tracks_nSatMeas",            &tracks_nSatMeas,           "tracks_nSatMeas/F") ;
  //fractNtuple->Branch("tracks_nMeas",            &tracks_nMeas,           "tracks_nMeas/F") ;
  fractNtuple->Branch("tracktruth_pt",            &tracktruth_pt,           "tracktruth_pt/F") ;
  fractNtuple->Branch("tracktruth_p",            &tracktruth_p,           "tracktruth_p/F") ;
  fractNtuple->Branch("tracktruth_e",            &tracktruth_e,           "tracktruth_e/F") ;
  fractNtuple->Branch("tracktruth_id",            &tracktruth_id,           "tracktruth_id/F") ;
  fractNtuple->Branch("tracktruth_type",            &tracktruth_type,           "tracktruth_type/l") ;
  fractNtuple->Branch("tracktruth_isPileup",            &tracktruth_isPileup,         "tracktruth_isPileup/O") ;
  fractNtuple->Branch("tracktruth_thesamewith",            &tracktruth_thesamewith,           "tracktruth_thesamewith/I") ;

  fractNtuple->Branch("CSCEndCapPlus",            &CSCEndCapPlus,           "CSCEndCapPlus/O") ;

#define MakeBranchAllSt(Name,Type,Var) sprintf(BranchName,"%s%d",Name,st+1); \
  sprintf(BranchTitle,"%s%d/%s",Name,st+1,Type);			\
  fractNtuple->Branch(BranchName, &Var[st],BranchTitle)
  
  for (Short_t st=0;st<4;st++) {
    char BranchName[30],BranchTitle[30];
    /*CSC Chamber Candidates in each station*/
    MakeBranchAllSt("CSCRg","S",CSCRg);
    MakeBranchAllSt("CSCCh","S",CSCChCand);
    MakeBranchAllSt("CSCCBad","O",CSCChBad);

    /*Extrapolated Tracks on CSC Chamber Candidates in each station*/
    MakeBranchAllSt("CSCDyProjHVGap","F",CSCDyProjHVGap);
    MakeBranchAllSt("CSCDyErrProjHVGap","F",CSCDyErrProjHVGap);
    MakeBranchAllSt("CSCProjDistEdge","F",CSCProjEdgeDist);
    MakeBranchAllSt("CSCProjDistErrEdge","F",CSCProjEdgeDistErr);
 
    /*Segments characteristics*/
    MakeBranchAllSt("CSCSegxLc","F",CSCSegxLc);
    MakeBranchAllSt("CSCSegyLc","F",CSCSegyLc);
    MakeBranchAllSt("CSCSegxErrLc","F",CSCSegxErrLc);
    MakeBranchAllSt("CSCSegyErrLc","F",CSCSegyErrLc);
    MakeBranchAllSt("CSCdXdZTTSeg","F",CSCdXdZTTSeg);
    MakeBranchAllSt("CSCdYdZTTSeg","F",CSCdYdZTTSeg);
    MakeBranchAllSt("CSCSegChisqProb","F",CSCSegChisqProb);
    MakeBranchAllSt("CSCnSegHits","I",CSCnSegHits);

    /*Distance from the Extrapolated Tracks to CSC Segments, 99999. for no CSC segment found*/ 
    MakeBranchAllSt("CSCDxTTSeg","F",CSCDxTTSeg);
    MakeBranchAllSt("CSCDxErrTTSeg","F",CSCDxErrTTSeg);
    MakeBranchAllSt("CSCDyTTSeg","F",CSCDyTTSeg);
    MakeBranchAllSt("CSCDyErrTTSeg","F",CSCDyErrTTSeg);
    MakeBranchAllSt("CSCDxyTTSeg","F",CSCDxyTTSeg);
    MakeBranchAllSt("CSCDxyErrTTSeg","F",CSCDxyErrTTSeg);

    /*LCT characteristics*/
    MakeBranchAllSt("CSCLCTxLc","F",CSCLCTxLc);
    MakeBranchAllSt("CSCLCTyLc","F",CSCLCTyLc);
    MakeBranchAllSt("CSCLCTbx","I",CSCLCTbx);

    /*Distance from the Extrapolated Tracks to LCT, 99999. for no LCT found*/
    MakeBranchAllSt("CSCDxTTLCT","F",CSCDxTTLCT);
    MakeBranchAllSt("CSCDxErrTTLCT","F",CSCDxErrTTLCT);
    MakeBranchAllSt("CSCDyTTLCT","F",CSCDyTTLCT);
    MakeBranchAllSt("CSCDxErrTTLCT","F",CSCDxErrTTLCT);
    MakeBranchAllSt("CSCDxyTTLCT","F",CSCDxyTTLCT);
    MakeBranchAllSt("CSCDxyErrTTLCT","F",CSCDxyErrTTLCT);

    /*DetlaR between the extrapolated tracker track on muon system and the tagged muon*/
    MakeBranchAllSt("dRTkMu","F",dRTkMu);

    /*Default decision of whether a segment or LCT is found*/
    //    MakeBranchAllSt("segSt","I",segSt);
    //    MakeBranchAllSt("lctSt","I",lctSt);
  }

  HLTMuAcceptance=new vector<Bool_t>();
  fractNtuple->Branch("HLTMuAcceptance", &HLTMuAcceptance);
  fractNtuple->Branch("HLTDiMuAcceptance",     &HLTDiMuAcceptance,       "HLTDiMuAcceptance/F") ;
  minDRHLTMu=new vector<Float_t>();
  fractNtuple->Branch("minDRHLTMu",       &minDRHLTMu) ;
  fractNtuple->Branch("minDRHLTDiMu",     &minDRHLTDiMu,       "minDRHLTDiMu/F") ;
  fractNtuple->Branch("minDRHLTAllSingleMu",     &minDRHLTAllSingleMu,       "minDRHLTAllSingleMu/F") ;

  RunInfo = new TTree("RunInfo","RunInfo");
  RunInfo->Branch("RunNumber",&run_number,"RunNum/i");
  RunInfo->Branch("TableName",&HLTTableName);
  HLTMuNames=new vector<string>();
  RunInfo->Branch("HLTMuNames",&HLTMuNames);
  HLTMuObjModuleNames=new vector<string>();
  RunInfo->Branch("HLTMuObjModuleNames",&HLTMuObjModuleNames);
  RunInfo->Branch("HLTDiMuName",&HLTDiMuName);
  RunInfo->Branch("HLTDiMuObjModuleName",&HLTDiMuObjModuleName);
  RunInfo->Branch("BadCSCChamberIndexList",&badChambersIndices);
}

TPTrackMuonSys::~TPTrackMuonSys(){
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called to for each event  ------------
void 
TPTrackMuonSys::analyze(const edm::Event& event, const edm::EventSetup& setup){
  //cout << "HERE *************************" << endl;

  nEventsAnalyzed++;
  // access conditions data for this event 
  theDbConditions.initializeEvent( setup ); //initializeEvent(const edm::EventSetup & es);// fetch the maps from the database
 
  //////////////////////////////////////////////////////////////////////////
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  Bool_t beamExists = false;
  try{
    event.getByLabel(m_beamSpot,recoBeamSpotHandle);
    if ( recoBeamSpotHandle.isValid()) {
      beamSpot = *recoBeamSpotHandle;
      beamExists = true;
    }
  }catch (cms::Exception){
    edm::LogError("")<< "Error! Can't get m_beamSpot  by label. ";
  }
  event_number = event.id().event();
  LumiBlock    = (int)event.luminosityBlock();
  expType      = (int)event.experimentType();
  bunchX       = (int)event.bunchCrossing();
  orbitNumb    = (int)event.orbitNumber();
  Bool_t isRData      = event.isRealData(); //( L1Decision ? "passed" : "failed")

  //Get the Magnetic field from the setup
  setup.get<IdealMagneticFieldRecord>().get(theBField);
  // Get the GlobalTrackingGeometry from the setup
  // setup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  // Get the DT Geometry from the setup 
  // setup.get<MuonGeometryRecord>().get(dtGeom);

  setup.get<MuonGeometryRecord>().get(cscGeom);
//#define m_debug
#ifdef m_debug
  //check the chamber active region geometry and numbering --- something good to know
  if (nEventsAnalyzed==1) {
    Short_t rings[]={11,12,21,31,41,22};//,14,12,13,21,31,41,22};
    for (Byte_t iendcap=1;iendcap<3;iendcap++) 
      for (Byte_t iring=0;iring<sizeof(rings)/sizeof(Short_t);iring++) {
	CSCDetId id=CSCDetId(iendcap, rings[iring]/10, rings[iring]%10, 1, 1);
	printf("============ ME%c%d chamber outline=============\n",iendcap==1?'+':'-',rings[iring]);
	const CSCLayerGeometry *layerGeom = cscGeom->chamber(id)->layer(1)->geometry ();
	const Byte_t NStrips=layerGeom->numberOfStrips(),NWires=layerGeom->numberOfWireGroups();
	LocalPoint interSect_ = layerGeom->stripWireGroupIntersection(1, 5);
	printf("(strip 1, wire group 5) at (%.3f,%.3f)\n",interSect_.x(),interSect_.y());
	interSect_ = layerGeom->stripWireGroupIntersection(NStrips, 5);
	printf("(strip %d, wire group 5) at (%.3f,%.3f)\n",NStrips,interSect_.x(),interSect_.y());
	interSect_ = layerGeom->stripWireGroupIntersection(NStrips, NWires);
	printf("(strip %d, wire group %d) at (%.3f,%.3f)\n",NStrips, NWires,interSect_.x(),interSect_.y());
	interSect_ = layerGeom->stripWireGroupIntersection(1, NWires);
	printf("(strip 1, wire group %d) at (%.3f,%.3f)\n", NWires,interSect_.x(),interSect_.y());
	printf("   ======== middle x of wires(check)==========\n");
	const CSCWireTopology* wireTopology = layerGeom->wireTopology();
	for (Short_t wireGroup_id=10;wireGroup_id<=15;wireGroup_id++) {
	  Float_t wideWidth      = wireTopology->wideWidthOfPlane();
	  Float_t narrowWidth    = wireTopology->narrowWidthOfPlane();
	  Float_t length         = wireTopology->lengthOfPlane();
	  Float_t tangent = (wideWidth-narrowWidth)/(2.*length);
	  Float_t wireangle= wireTopology->wireAngle();
	  std::vector<float> wirecenters=wireTopology->wireValues(wireTopology->middleWireOfGroup(wireGroup_id));
	  printf("x center of wire %d is at %.3f=%.3f\n",wireGroup_id,wireTopology->wireValues(wireTopology->middleWireOfGroup(wireGroup_id))[0],wirecenters[2]*sin(wireangle)*tangent/2);
	}
	printf("y center of the wires is at %.3f\n",wireTopology->yOfWire(1)+0.5*wireTopology->lengthOfPlane());
      }
  }
#endif
  // Get the propagators
  setup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyRK", propagatorAlong   );
  setup.get<TrackingComponentsRecord>().get("SmartPropagatorAnyOpposite", propagatorOpposite);

  //vertices
  edm::Handle<reco::VertexCollection> pvHandle; 
  try{
    event.getByLabel(m_vertexSrc, pvHandle);
  }catch (cms::Exception){
    edm::LogError("")<< "Error! Can't get m_vertexSrc by label. ";
  }  
  const reco::VertexCollection & vertices = *pvHandle.product();	 
  numberOfPrimaryVertices =pvHandle ->size();

  // Get the RecHits collection :
  Handle<CSCRecHit2DCollection> recHits;
  event.getByLabel("csc2DRecHits",recHits);

  // get CSC segment collection
  event.getByLabel("cscSegments", cscSegments); 

  edm::Handle<CSCCorrelatedLCTDigiCollection> mpclcts;
  try{
    event.getByLabel("csctfunpacker","", mpclcts);
  }catch (cms::Exception){
    edm::LogError("")<< "Error! Can't get m_gTracksTag by label. ";
  }
   
  /*
  // get RPC and DT
  setup.get<MuonGeometryRecord>().get(rpcGeo);
   
  edm::Handle<RPCRecHitCollection> rpcRecHits;
  event.getByLabel("rpcRecHits",rpcRecHits);

  edm::Handle<DTRecSegment4DCollection> dtSegments;
  event.getByLabel( "dt4DSegments", dtSegments );
  */

  //  muon trigger scales
  edm::ESHandle<L1MuTriggerScales> trigscales_h;
  setup.get<L1MuTriggerScalesRcd> ().get(trigscales_h);
  theTriggerScales = trigscales_h.product();

  //if (! dtSegments.isValid()) throw cms::Exception("FatalError") << "Unable to find DTRecSegment4DCollection in event!\n";

  std::vector<DetId> chamberIds;
  std::vector<DetId>::const_iterator chamberIdIt;

  Handle<reco::TrackCollection> gTracks;
  // Bool_t trkBool[2]; trkBool[0]= true;
  try{
    event.getByLabel(m_gTracksTag, gTracks); 
  }catch (cms::Exception){
    edm::LogError("")<< "Error! Can't get m_gTracksTag by label. ";
    // trkBool[0] = false;
  }

  try{
    event.getByLabel("muons", muons); 
  }catch (cms::Exception){
    edm::LogError("")<< "Error! Can't get muons ";
  }
  Handle<ValueMap<DeDxData> > dEdxTrackHandle;
  event.getByLabel(m_dEdxDiscrimTag, dEdxTrackHandle);
  const ValueMap<DeDxData> dEdxTrack = *dEdxTrackHandle.product();

  // calo geometry
  edm::ESHandle<CaloGeometry> calogeo;
  setup.get<CaloGeometryRecord>().get(calogeo);
  if (! calogeo.isValid()) throw cms::Exception("FatalError") << "Unable to find CaloGeometryRecord in event!\n";
  
  //=======Generation and Simulation information========
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(m_gTracksTag,trackCollectionHV);
  //Simulated Vertices: the vertexId() is just the index
  Handle<SimVertexContainer> SVCollectionH;
  HepMC::GenEvent * HepGenEvent=NULL;
  
  //TrackingParticles
  Handle<TrackingParticleCollection> TPCollectionH ;
  RecoToSimCollection RecoToSimByHits;
  ESHandle<TrackAssociatorBase> AssociatorByHits;

  vector< vector< vector<Int_t> > > TracksSimChains, MuTracksSimChains;

  numberOfPUVertices=0; // the number of pileup interactions that have been added to the event from this bunch crossing
  numberOfPUVerticesMixingTruth=0;// the *true* mean number of pileup interactions for this event from which each bunch crossing has been sampled; same for all bunch crossings in an event
  numberOfPUVerticesTot=0;
  mcweight = 1.;
  Handle<edm::HepMCProduct> HepMCH;
  if (m_isMC) {
    event.getByLabel("g4SimHits", SVCollectionH);
    SVC = *SVCollectionH.product();
    event.getByLabel(m_HepMCTag, HepMCH);
    HepGenEvent=const_cast<HepMC::GenEvent *>( HepMCH->GetEvent() );
    Handle<GenEventInfoProduct> hEvtInfo;
    event.getByLabel("generator", hEvtInfo);
    mcweight = hEvtInfo->weight();
    event.getByLabel("mergedtruth",TPCollectionH);
    //event.getByType(TPCollectionH);
    setup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", AssociatorByHits);
    RecoToSimByHits = AssociatorByHits->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&setup);
    MCParticlesList.clear();
    SavedSimTrk.clear();
    SavedHepPar.clear();
    /*Start of getting pileup information*/
    edm::InputTag PileupSrc_("addPileupInfo");
    Handle<std::vector< PileupSummaryInfo > >  puInfo;
    if ( event.getByLabel(PileupSrc_, puInfo) ) {
#ifdef GetPUFromEarlierThan_4_1_2
      numberOfPUVertices = (*puInfo)[0].getPU_NumInteractions();
      numberOfPUVerticesMixingTruth=numberOfPUVertices;
      numberOfPUVerticesTot=numberOfPUVertices;
#else
      for(std::vector<PileupSummaryInfo>::const_iterator PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
	Int_t BX = PVI->getBunchCrossing();
	if(BX == 0) {
	  numberOfPUVertices = PVI->getPU_NumInteractions();
#ifdef GetPUFromEarlierThan_4_4_0
	  numberOfPUVerticesMixingTruth = -1.;//not available
#else
	  numberOfPUVerticesMixingTruth = PVI->getTrueNumInteractions();
#endif
	}
	numberOfPUVerticesTot+=numberOfPUVertices;
      }
#endif
    }//====== End of getting pileup information =============
  }//======= End of if m_isMC========

  //////////////////////////////////////////////   
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  setup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();
  Handle<L1GlobalTriggerReadoutRecord> gtBeforeMaskRecord;
  event.getByLabel("gtDigis", gtBeforeMaskRecord);
  trgSingle  = false;
  if (!gtBeforeMaskRecord.isValid()) {
    LogInfo("PhysicsTrees: ") << " Error: no L1GlobalTriggerReadoutRecord found with input tag --> gtDigis" ;
  } else {
    const TechnicalTriggerWord techDecisionWord = gtBeforeMaskRecord->technicalTriggerWord(0);   
    const DecisionWord gtDecisionWordBeforeMask = gtBeforeMaskRecord->decisionWord(0);
#ifdef m_debug       
    Int_t nl1s = 0;
#endif
    for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
#ifdef m_debug       
      if (nEventsAnalyzed <= 1) {
	cout <<" L1Menu Item m1 = " << nl1s 
	     <<" corresponds to AlgoName = " 
	     << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias() 
	     << " and the result= " << menu->gtAlgorithmResult((algo->second).algoAlias(), gtDecisionWordBeforeMask)
	     << endl;
      }
      nl1s = nl1s + 1;
#endif
      if((algo->second).algoName() == "L1_SingleMuOpen") trgSingle  = true;
      if((algo->second).algoName() == "L1_SingleMu5") trgSingle=false;
    }
  }

  /*-----------Start getting HLT results------------*/ 
  Handle<TriggerResults> triggerResults;
  Handle< trigger::TriggerEvent > handleTriggerEvent;
  HLTDiMuAcceptance=false;
  HLTMuAcceptance->clear();
  try {
    event.getByLabel(m_hlt, triggerResults);
    if ( triggerResults.product()->wasrun() ){
      LogInfo("")<<" At least one path out of " << triggerResults.product()->size() << " ran? " << triggerResults.product()->wasrun() << endl;	
      if ( triggerResults.product()->accept() ) {
	for (vector<Int_t>::const_iterator iter=m_HLTMuTrgBit.begin();iter<m_HLTMuTrgBit.end();iter++)
	  HLTMuAcceptance->push_back( triggerResults.product()->accept(*iter) );
	if ( m_HLTDiMuTrgBit>-1 ) HLTDiMuAcceptance=triggerResults.product()->accept( m_HLTDiMuTrgBit );
      } // end if at least one triggerResult accepted
    }// end if wasRun
  } catch (...) {// some old codes, it is possibly not going to happen in >CMSSW_4_4_X versions
    LogWarning("")<<"Not found HLT result, is HLT process name correct?"<< endl;	
    if(isRData == 0){
      if(event.getByLabel( m_hltTrgEv, handleTriggerEvent )){
	const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
	for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
	  std::string fullname = handleTriggerEvent->filterTag(ia).encode();
	  std::string name;
	  size_t p = fullname.find_first_of(':');
	  if ( p != std::string::npos) name = fullname.substr(0, p);
	  else name = fullname;
	  if ( &toc !=0 ) {
	    HLTMuAcceptance->assign(m_HLTMuTrgBit.size(),false);
	    for (vector<string>::const_iterator iter=HLTMuObjModuleNames->begin(); iter<HLTMuObjModuleNames->end();iter++ )
	      if ( name == *iter ) HLTMuAcceptance->at( iter-HLTMuObjModuleNames->begin() )=true;
	    if ( name == *HLTDiMuObjModuleName ) HLTDiMuAcceptance=true;
	  }
	}
      }
    }
  }
  Bool_t m_GotTrgObj=event.getByLabel( m_hltTrgEv, handleTriggerEvent );
  /*-----------End getting HLT results------------*/ 

  if(!gTracks.isValid() ) return;
  if(!muons.isValid() ) return;
  if (gTracks->size() > MAXNTRACKS) return;

#ifdef  m_debug
  cout << "Run: " << run_number << " Event: " << event_number << " LumiBlock: " <<  LumiBlock << " BX: " <<  bunchX 
       << " isRealDAT: " <<  isRData 
       << endl;
#endif

  Int_t trackNumber[2][4][4][36];
  Bool_t trkVeto[MAXNTRACKS];

  for(Int_t i1 = 0; i1 < 2; i1++)
    for(Int_t i2 = 0; i2 < 4; i2++)
      for(Int_t i3 = 0; i3 < 4; i3++)
	for(Int_t i4 = 0; i4 < 36; i4++)
	  trackNumber[i1][i2][i3][i4] = -1;
  
  for(Int_t i1 = 0; i1 < MAXNTRACKS; i1++)
    trkVeto[i1] = false;

  nPosTrk = 0;
  nNegTrk = 0;
  nTotalTrks = gTracks->size();
  for(reco::TrackCollection::const_iterator itTrack = gTracks->begin(); itTrack != gTracks->end(); itTrack++){//start loop tracks
    UInt_t itrk=itTrack - gTracks->begin();
    if(itTrack->charge() == 0) continue;
    if(itTrack->p() < 3.0)  continue;

    reco::TrackRef trackRef(gTracks, itrk );
    tracks_eta    = itTrack->eta(); tracks_phi = itTrack->phi();
    if(tracks_phi < 0 )tracks_phi = tracks_phi + 2*M_PI;
    tracks_chi2   = 9999.;
    if(itTrack->ndof() != 0)tracks_chi2   = itTrack->chi2()/itTrack->ndof();
    tracks_dxy  = itTrack->dxy(); tracks_dz   = itTrack->dz();
    if(beamExists){
      tracks_dxy  = itTrack->dxy(beamSpot.position()); 
      tracks_dz   = itTrack->dz(beamSpot.position());  
    }

    tracks_numberOfValidHits = itTrack->numberOfValidHits();
    if ( ! (tracks_numberOfValidHits>7  && fabs(tracks_dz)<24 && fabs(tracks_dxy)< 2.0 && tracks_chi2>0 && tracks_chi2<4) ) continue;
    Bool_t ec= ( tracks_eta > 0 );
    //     if(fabs(tracks_eta) < 1.2 ) trackDT = true;
    Bool_t GotCSCSegMatched = false;
    if(fabs(tracks_eta) > 0.9 )//if it is a CSCTrack
      for(UInt_t j =0; j < 6; j++){//start of matched segments check
	/* not necessary codes, it will give  the same result
	   Float_t zzPlaneME = 0.0;
	   if(trackCSC && MEZ[j] != 0){
	   zzPlaneME = MEZ[j];
	   if (!ec) zzPlaneME = -MEZ[j];
	   TrajectoryStateOnSurface tsos = surfExtrapTrkSam(trackRef, zzPlaneME);  
	   if (!tsos.isValid()) continue;
	   Float_t trkEta = tsos.globalPosition().eta(), trkPhi = tsos.globalPosition().phi();
	   Int_t rg = ringCandidate(j>2?j-1:1, trkEta, trkPhi);
	   if (rg==0) continue;
	   if( thisChamberCandidate(st, rg, trkPhi) %2 == 0 ){
	   zzPlaneME = MEZEven[j];
	   if(!ec)zzPlaneME = -MEZEven[j];
	   tsos = surfExtrapTrkSam(trackRef, zzPlaneME);  
	   }else{
	   zzPlaneME = MEZOdd[j];
	   if(!ec)zzPlaneME = -MEZOdd[j];
	   tsos = surfExtrapTrkSam(trackRef, zzPlaneME); 
	   }
	   if (!tsos.isValid())continue;*/
	CSCSegmentCollection::const_iterator cscSegOut;
	if( !matchTTwithCSCSeg(ec, j, trackRef, cscSegments,  cscSegOut) ) continue;
	GotCSCSegMatched=true;
	CSCDetId id  = (CSCDetId) cscSegOut->cscDetId();		   
	Byte_t endcapCSC = id.endcap() - 1,ringCSC = id.ring() - 1,stationCSC = id.station() - 1,chamberCSC = id.chamber() - 1;
 	if ( trackNumber[endcapCSC][stationCSC][ringCSC][chamberCSC] < 0 )
	  trackNumber[endcapCSC][stationCSC][ringCSC][chamberCSC] =  itrk;
	else {
	  trkVeto[ trackNumber[endcapCSC][stationCSC][ringCSC][chamberCSC] ]=true;
	  trkVeto[ itrk ]=true;
	}
      }//end of matched segments check
    
    if (GotCSCSegMatched) {
      if(ec){
	nPosTrk ++;
      }else{
	nNegTrk ++;
      }
    }
  }//end of first track loop

  std::string trackExtractorName = trackExtractorPSet_.getParameter<std::string>("ComponentName");
  reco::isodeposit::IsoDepositExtractor* muIsoExtractorTrack_ = IsoDepositExtractorFactory::get()->create( trackExtractorName, trackExtractorPSet_);
  for (reco::MuonCollection::const_iterator muIter1 = muons->begin(); muIter1 != muons->end(); ++muIter1) {  
    if (!muIter1->isTrackerMuon() ) continue;
    if(!muIter1->track().isNonnull()) continue;
    MuTagPt     = muIter1->track()->pt() ;
    if (MuTagPt < 5.) continue;
    
    Float_t mu1dxy = 1000., mu1dz = 1000.;

    if(beamExists){
      mu1dxy = muIter1->track()->dxy(beamSpot.position()); 
      mu1dz = muIter1->track()->dz(beamSpot.position()); 
    }else{
      mu1dxy = muIter1->track()->dxy(); 
      mu1dz = muIter1->track()->dz(); 
    }
    Float_t mu1Chi2 = muIter1->track()->normalizedChi2();
    MuTagHitsTrkSys  = muIter1->track()->hitPattern().numberOfValidTrackerHits();
    Bool_t goodTrack  = (fabs(mu1dxy) < 2.0 && fabs(mu1dz) < 24.0 && fabs(mu1Chi2) < 4.0 && MuTagHitsTrkSys > 7 ); //&& MuTagHitsMuSys > 3 
    if(!goodTrack)continue;
    
    MuTagnSegTrkArb  = muIter1->numberOfMatches(); // get number of chambers with matched segments
    if (MuTagnSegTrkArb < 2) continue;

    MuTagPx    = muIter1->track()->px() ;
    MuTagPy    = muIter1->track()->py() ;
    MuTagPz    = muIter1->track()->pz() ;
 
    MuTagE      = muIter1->track()->p() ;
    MuTagEta    = muIter1->track()->eta() ;
    MuTagPhi    = muIter1->track()->phi() ;

    Int_t MuTagcharge = muIter1->track()->charge() ;

    //    MuTagHitsPixSys  = muIter1->track()->hitPattern().numberOfValidPixelHits();
    //    MuTagHitsRPCSys  = muIter1->track()->hitPattern().numberOfValidMuonRPCHits();

    /// check for HLT matching
    minDRHLTDiMu=100.; minDRHLTAllSingleMu=100.;
    minDRHLTMu->assign(HLTMuObjModuleNames->size(),100.);
    Float_t PhiTemp1 = MuTagPhi<0?MuTagPhi + 2*M_PI:MuTagPhi;
    if ( m_GotTrgObj ){
      const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
      for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
	std::string fullname = handleTriggerEvent->filterTag(ia).encode();
	std::string name;
	size_t p = fullname.find_first_of(':');
	if ( p != std::string::npos) name = fullname.substr(0, p);
	else name = fullname;
	if ( &toc ==0) continue;
	const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
	for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
	  double l3phi  = toc[*ki].phi();
	  double l3eta  = toc[*ki].eta();
	  
	  if(l3phi < 0 )l3phi = l3phi + 2*M_PI;

	  Float_t deltaR_TTHLT   = sqrt((l3eta - MuTagEta) * (l3eta - MuTagEta) + (l3phi - PhiTemp1) *(l3phi - PhiTemp1));
	  if ( deltaR_TTHLT > 0.4 ) continue;
	  if ( name == *HLTDiMuObjModuleName ) {
	    if ( deltaR_TTHLT < minDRHLTDiMu ){
	      minDRHLTDiMu = deltaR_TTHLT;
	    }
	  }
	  Byte_t idx=0;
	  for (vector<string>::const_iterator iter=HLTMuObjModuleNames->begin(); iter<HLTMuObjModuleNames->end();iter++,idx++ )
	    if ( name == *iter ) {
	      if ( deltaR_TTHLT < minDRHLTMu->at(idx) ){
		minDRHLTMu->at(idx) = deltaR_TTHLT;
	      }
	      if ( deltaR_TTHLT < minDRHLTAllSingleMu ){
		minDRHLTAllSingleMu = deltaR_TTHLT;
	      }
	      break;
	    }
	}//end of loop ki
      }//end of loop ia
    }//end of if m_GotTrgObj
    MuTagIsoR03Ratio = MuTagPt>0?muIter1->isolationR03().sumPt/MuTagPt:9999.;
    MuTagIsoR05Ratio = MuTagPt>0?muIter1->isolationR05().sumPt/MuTagPt:9999.;
    MuTagCaloL= muon::isGoodMuon(*muIter1,muon::TM2DCompatibilityLoose);
    MuTagCaloT= muon::isGoodMuon(*muIter1,muon::TM2DCompatibilityTight);
     
    /*------------------------Start getting the Monte Carlo Truth--------------------------*/
    MuTagtracktruth_pt  = -9999.;
    MuTagtracktruth_p   = -9999.;
    MuTagtracktruth_id  = -9999.;
    MuTagtracktruth_type = 0;
    MuTagtracktruth_isPileup=false;
    MuTagtracktruth_thesamewith=-1;
    if ( m_isMC) {
      //TrackingParticles
      MuTracksSimChains.push_back(* new vector< vector<Int_t> >() );
      if ( TPCollectionH.isValid() ) {
	//SimToReco Tracks Association
	reco::TrackBaseRef trackBaseRef(muIter1->track());
	if ( RecoToSimByHits.find(trackBaseRef) != RecoToSimByHits.end() ) {
	  pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trackBaseRef].front();
	  //vector<pair<TrackingParticleRef, double> > TPCByChi2=RecoToSimByChi2[trk];
	  TrackingParticleRef tpr = BestMatch.first;
	  //Simulated Tracks
	  MuTagtracktruth_pt=tpr->pt();
	  MuTagtracktruth_p=tpr->p();
	  MuTagtracktruth_id=tpr->pdgId();
	  MuTagtracktruth_isPileup=GetDecayChains(tpr,HepGenEvent,MuTagtracktruth_type,MuTagtracktruth_thesamewith, MuTracksSimChains);
	}
      }//end of TPCollection.isValid
    }//end of m_isMC
    /*------------------------End of getting the Monte Carlo Truth--------------------------*/
   
    CLHEP::Hep3Vector Mu1Vector(MuTagPx, MuTagPy, MuTagPz);
      
    Float_t etaMuVec[4] ={-9999., -9999., -9999., -9999.};
    Float_t phiMuVec[4] ={-9999., -9999., -9999., -9999.};

    Float_t allCSC_mu = MEZ[0];
    if(fabs(MuTagEta) < 1.5)allCSC_mu = MEZ[1];
    if(MuTagEta < 0) allCSC_mu = -allCSC_mu;

    reco::TrackRef trackMuRef = muIter1->track();
    TrajectoryStateOnSurface tsos = surfExtrapTrkSam(trackMuRef, allCSC_mu);
    if(tsos.isValid()){
      etaMuVec[0] = tsos.globalPosition().eta();
      phiMuVec[0] = tsos.globalPosition().phi();       
    }
    allCSC_mu = MEZ[3];
    if(MuTagEta < 0) allCSC_mu = -allCSC_mu; 
    tsos = surfExtrapTrkSam(trackMuRef, allCSC_mu);
    if(tsos.isValid()){
      etaMuVec[1] = tsos.globalPosition().eta();
      phiMuVec[1] = tsos.globalPosition().phi();       
    }
     
    allCSC_mu = MEZ[4];
    if(MuTagEta < 0) allCSC_mu = -allCSC_mu; 
    tsos = surfExtrapTrkSam(trackMuRef, allCSC_mu);
    if(tsos.isValid()){
      etaMuVec[2] = tsos.globalPosition().eta();
      phiMuVec[2] = tsos.globalPosition().phi();       
    }
     
    allCSC_mu = MEZ[5];
    if(MuTagEta < 0) allCSC_mu = -allCSC_mu; 
    tsos = surfExtrapTrkSam(trackMuRef, allCSC_mu);
    if(tsos.isValid()){
      etaMuVec[3] = tsos.globalPosition().eta();
      phiMuVec[3] = tsos.globalPosition().phi();       
    }

    Bool_t firsttrackmatchingtoMuTag=true;

    for(reco::TrackCollection::const_iterator itTrack = gTracks->begin(); itTrack != gTracks->end(); itTrack++){
      UInt_t itrk = itTrack - gTracks->begin();
      if(itTrack->charge() == 0) continue;
      if(itTrack->p() < 3.0)  continue;
      trackVeto = trkVeto[itrk]; 
      reco::TrackRef trackRef(gTracks, itrk );

      if(itTrack->charge()*MuTagcharge != -1) continue;

      tracks_eta    = itTrack->eta();
      
      if(fabs(tracks_eta) > 0.9 ) {
	if(fabs(tracks_eta) < 1.2 ) myRegion = 2;
	else myRegion = 3;
      }
      else myRegion = 1;
      if (myRegion == 1) continue;//currently, only CSC

      tracks_phi = itTrack->phi();
      if(tracks_phi < 0 )tracks_phi = tracks_phi + 2*M_PI;

      tracks_chi2   = -9999.;
      if(itTrack->ndof() != 0)tracks_chi2   = itTrack->chi2()/itTrack->ndof();
      tracks_charge = itTrack->charge();

      tracks_dxy = itTrack->dxy(); tracks_dz = itTrack->dz();
      if(beamExists){
        tracks_dxy  = itTrack->dxy(beamSpot.position()); tracks_dz = itTrack->dz(beamSpot.position());
      }

      MuProbenHitsTrkSys  = itTrack->hitPattern().numberOfValidTrackerHits();
      MuProbenHitsPixSys  = itTrack->hitPattern().numberOfValidPixelHits();
       
      goodTrack = ( fabs(itTrack->eta())< 2.4 && fabs(tracks_dz)< 24.0 && 
		    fabs(tracks_dxy)< 2.0 && tracks_chi2> 0.0 &&  tracks_chi2< 4.0 && MuProbenHitsTrkSys > 7 );

      if (!goodTrack) continue;
      if(tracks_eta > 0) CSCEndCapPlus = true;
      else CSCEndCapPlus = false;
      
      tracks_pt    = itTrack->pt();
      tracks_e = itTrack->p(); 
      tracks_vx     = itTrack->vx();  // x coordinate of the reference point on track
      tracks_vy     = itTrack->vy();  // y coordinate of the reference point on track
      tracks_vz     = itTrack->vz();  // z coordinate of the reference point on track
      reco::IsoDeposit depTrk = muIsoExtractorTrack_->deposit(event, setup, *itTrack );
      tracks_IsoR03Ratio = tracks_pt>0?depTrk.depositWithin(0.3)/tracks_pt:9999.;
      tracks_IsoR05Ratio = tracks_pt>0?depTrk.depositWithin(0.5)/tracks_pt:9999.;
      tracks_qoverp = itTrack->qoverp(); // q/p 
      tracks_lambda = itTrack->lambda();
      tracks_recHitsSize= itTrack->recHitsSize();
      tracks_qoverpError = itTrack->qoverpError();//TrackAlgorithm tAlg    = itTrack->algo();
      tracks_ptError = itTrack->ptError();/// error on Pt (set to 1000 TeV if charge==0 for safety)
      tracks_thetaError = itTrack->thetaError();/// error on theta
      tracks_lambdaError = itTrack->lambdaError(); /// error on lambda
      tracks_etaError = itTrack->etaError(); tracks_phiError = itTrack->phiError(); tracks_dxyError = itTrack->dxyError();/// error on dxy
      tracks_d0Error = itTrack->d0Error(); tracks_dszError = itTrack->dszError();tracks_dzError = itTrack->dzError();/// error on dz
      tracks_numberOfValidHits = itTrack->numberOfValidHits(); ///unsigned short   number of valid hits found 
      tracks_numberOfLostHits = itTrack->numberOfLostHits(); ///unsigned short  number of cases where track crossed a layer without getting a hit.
      if(tracks_phi < 0)tracks_phi = tracks_phi + 2*M_PI;
      reco::MuonCollection::const_iterator matchedMu=matchTTwithMT(itTrack);
      if ( matchedMu!=muons->end() ) {
	tracks_isCaloMuTrk=matchedMu->isCaloMuon();
	tracks_isTrackerMuTrk=matchedMu->isTrackerMuon();
      }
      else {
	tracks_isCaloMuTrk=false;
	tracks_isTrackerMuTrk=false;
      }

      Bool_t trQuality = (fabs(MuTagEta) < 2.4 && fabs(tracks_eta) > 0.9 && fabs(tracks_eta) < 2.4 
			  && tracks_etaError < 0.003 && tracks_phiError < 0.003 
			  && tracks_ptError/tracks_pt < 0.05 && tracks_numberOfValidHits >= 10); // cuts removed from the SkimDPG.C file and put here...

      if(!trQuality)continue;

      Float_t mMu = 0.1134289256;
      invMass = pow( ( sqrt(pow(itTrack->p(),2)+ mMu*mMu) +  sqrt(pow(muIter1->track()->p(),2)+ mMu*mMu) ) ,2 ) -
	(
	 pow((itTrack->px() + muIter1->track()->px()),2) +
	 pow((itTrack->py() + muIter1->track()->py()),2) +
	 pow((itTrack->pz() + muIter1->track()->pz()),2)
	 );
      
      if(invMass < 0) continue;
      invMass = sqrt(invMass);
      Bool_t gotMass =  ((invMass > 2.5 &&  invMass < 3.6) || (invMass > 75.)) ;
      if(!gotMass)continue;

      /*------------------------Start getting the Monte Carlo Truth--------------------------*/
      tracktruth_pt  = -9999.;
      tracktruth_e   = -9999.;
      tracktruth_p   = -9999.;
      tracktruth_id  = -9999.;
      tracktruth_type = 0;
      tracktruth_isPileup=false;
      tracktruth_thesamewith=-1;
      
      if ( m_isMC) {
	//TrackingParticles
	if ( TPCollectionH.isValid() ) {
#ifdef jz_debug
	  //	cout<<"track "<<i<<" has "<<tracks_numberOfValidHits << " hits" <<endl;
#endif
	  if  (!firsttrackmatchingtoMuTag) MuTracksSimChains.push_back( MuTracksSimChains.back() );
	  TracksSimChains.push_back(* new vector< vector<Int_t> >() );
	  //SimToReco Tracks Association
	  reco::TrackBaseRef trackBaseRef(trackRef);
	  if ( RecoToSimByHits.find(trackBaseRef) != RecoToSimByHits.end() ) {
	    pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trackBaseRef].front();
	    //vector<pair<TrackingParticleRef, double> > TPCByChi2=RecoToSimByChi2[trk];
	    TrackingParticleRef tpr = BestMatch.first;
	    //Simulated Tracks
	    tracktruth_pt=tpr->pt();
	    tracktruth_e=tpr->energy();
	    tracktruth_p=tpr->p();
	    tracktruth_id=tpr->pdgId();
	    tracktruth_isPileup=GetDecayChains(tpr,HepGenEvent,tracktruth_type,tracktruth_thesamewith,TracksSimChains);
	  }
	}//end of TPCollection.isValid
      }//end of m_isMC
      /*------------------------End of getting the Monte Carlo Truth--------------------------*/

      CLHEP::Hep3Vector trackVector(itTrack->px(), itTrack->py(), itTrack->pz());
      ////////////  Check if they have the same good vertex...

      Bool_t gotVertexTrk1 = false, gotVertexTrk2 = false;

      iSameVtx = false;
      vtx_r = -9999.; vtx_z = -9999.;
      vtx_rError = -9999.; vtx_zError = -9999.;
      vtx_normChi2 = 100000.; vtx_size = -1;

      for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it) {
	if(!it->isValid())continue;
	if(it->isFake())continue;
	for(std::vector<TrackBaseRef>::const_iterator itr=it->tracks_begin() ; itr!=it->tracks_end() ; ++itr) {
	  if ( & (*muIter1->track().get() ) == &( * const_cast<reco::Track *>(itr->get()) ) )  {
	    //	    cerr<<"Got MuTag trk matched"<<endl;
	    gotVertexTrk1 = true;
	    if (gotVertexTrk2) break;
	  }
	  if ( & (* itTrack) == & ( * (itr->get()) ) ) {
	    //    cerr<<"Got TT trk matched"<<endl;
	    gotVertexTrk2 = true;
	    if (gotVertexTrk1) break;
	  }
	}// tracks within vertex

	if(gotVertexTrk1 && gotVertexTrk2){
	  vtx_size = it->tracksSize();
	  vtx_r = sqrt(pow(it->x(),2)+ pow(it->y(),2)); vtx_z = it->z();
	  if(it->x() != 0 && it->y() != 0) vtx_rError = vtx_r*sqrt(pow((it->xError()/it->x()),2) + pow((it->yError()/it->y()),2));
	  if(it->x() != 0 && it->y() == 0) vtx_rError = it->xError();
	  if(it->x() == 0 && it->y() != 0) vtx_rError = it->yError();
	  vtx_zError = it->zError();
	  if(it->ndof() != 0) vtx_normChi2 = it->chi2()/it->ndof();
	  iSameVtx = true;
	}
      }// for loop over vertices..

      /////////// finished checking the vertex..

      nTrkCountCSCSeg=0;
      for(Byte_t j =0; j < 4; j++){
	/*CSC Chamber Candidates in each station*/
	CSCRg[j]=-9999;
	CSCChCand[j]=-9999;
 	CSCChBad[j] = false;
	/*Extrapolated Tracks on CSC Chamber Candidates in each station*/
	CSCDyProjHVGap[j]=9999.;
	CSCDyErrProjHVGap[j]=-9999.;
	CSCProjEdgeDist[j]=-9999.;
	CSCProjEdgeDistErr[j]=-9999.;
	/*Segments characteristics*/
	CSCSegxLc[j]=-9999.;
	CSCSegyLc[j]=-9999.;
	CSCSegxErrLc[j]=-9999.;
	CSCSegyErrLc[j]=-9999.;
 	CSCSegChisqProb[j]=-9999.;
	CSCdXdZTTSeg[j]=-9999.;
	CSCdYdZTTSeg[j]=-9999.;
	CSCnSegHits[j]=-9999;
	/*Distance from the Extrapolated Tracks to CSC Segments, 9999. for no CSC segment found*/ 
	CSCDxTTSeg[j] = -9999.;
	CSCDxErrTTSeg[j] = -9999.;
	CSCDyTTSeg[j] = -9999.;
	CSCDyErrTTSeg[j] = -9999.;
        CSCDxyTTSeg[j] = -9999.;
 	CSCDxyErrTTSeg[j] = -9999.;
	/*LCT characteristics*/
	CSCLCTxLc[j] = -9999.;
	CSCLCTyLc[j] = -9999.;
	CSCLCTbx[j] = -9999;
	/*Distance from the Extrapolated Tracks to LCT, 9999. for no LCT found*/
	CSCDxTTLCT[j] = -9999.;
	CSCDxErrTTLCT[j] = -9999.;
	CSCDyTTLCT[j] = -9999.;
	CSCDyErrTTLCT[j] = -9999.;
        CSCDxyTTLCT[j] = -9999.;
 	CSCDxyErrTTLCT[j] = -9999.;
	/*DetlaR between the extrapolated tracker track on muon system and the tagged muon*/
	dRTkMu[j] = -9999.;
	/*Default decision of whether a segment or LCT is found*/
	segSt[j]=false;
	lctSt[j]=false;
      }
 
      for(UInt_t j =0; j < 6; j++) {

	if( MEZ[j] == 0)continue;

	Float_t zzPlaneME = MEZ[j];
	if(!CSCEndCapPlus) zzPlaneME = -MEZ[j];
	tsos = surfExtrapTrkSam(trackRef, zzPlaneME);  
	if (!tsos.isValid()) continue;

	Int_t st = j>2?j-2:0, ec=CSCEndCapPlus?1:2; // determine the station number...
	Float_t trkEta = tsos.globalPosition().eta(), trkPhi = tsos.globalPosition().phi();

	Short_t rg = ringCandidate(st+1, trkEta, trkPhi);
	if ( rg==-9999) continue;
	if (!(
	      (j==0&&(rg==1||rg==4)) ||
	      (j==1&&rg==2) ||
	      (j==2&&rg==3) ||
	      (j>2)
	      )) continue;//check if it is runing on the right ring for a certain Z
	else CSCRg[st]=rg;

	//for chamber overlap region
	/* This part of code is useless. Just directly extrapolate to chamber
	cerr<< "CHK"<<thisChamberCandidate(st+1, rg, trkPhi);
	if ( thisChamberCandidate(st+1, rg, trkPhi)%2 == 0 ) zzPlaneME = MEZEven[j];
	else zzPlaneME = MEZOdd[j];
	 
	if(!CSCEndCapPlus)zzPlaneME = -zzPlaneME;	 

	tsos = surfExtrapTrkSam(trackRef, zzPlaneME);
	if (!tsos.isValid()) continue;
	trkEta = tsos.globalPosition().eta(); trkPhi =  tsos.globalPosition().phi();
	rg = ringCandidate(st+1, trkEta, trkPhi);
	if ( rg==-9999 ) continue;
	*/

	//Calculate the solid angle between the muon and this extrapolated track
	if(trkPhi  < 0) trkPhi += 2*M_PI;
	/*dR between muon and track*/
	Float_t MPhi = phiMuVec[st];
	if (MPhi < 0 ) MPhi +=  2*M_PI;
	Float_t TPhi = tsos.globalPosition().phi();
	if (TPhi < 0 ) TPhi += 2*M_PI;
	dRTkMu[st] = deltaR( etaMuVec[st], MPhi, tsos.globalPosition().eta() , TPhi);
//	printf("Mu(%f,%f),Trk(%f,%f)-->dR(%f)",etaMuVec[st], MPhi, tsos.globalPosition().eta() , TPhi,dRTkMu[st]);
	//Find the chamber candidate. If two chamber overlaps, simply divide them from middle. e.g. if chamber1 is from phi=1-3, chamber 2 is from phi=2-4. Then phi=1-2.5-->chamber 1; phi=2.5-4-->chamber 2. Since the distribution over phi is uniform, it won't bias the result.
	CSCChCand[st] = thisChamberCandidate(st+1, rg, trkPhi);
	CSCDetId Layer0Id=CSCDetId(ec, st+1, rg,  CSCChCand[st], 0);//layer 0 is the mid point of the chamber. It is not a real layer.
	CSCChBad[st] = badChambers_->isInBadChamber( Layer0Id );
	//skip not-existing ME42 chambers
	if (CSCChBad[st]&&st==3&&rg==2) continue;
#ifdef jz_debug
	if (CSCChBad[st]) cerr<<(CSCEndCapPlus?"ME+":"ME-")<<st+1<<"/"<<rg<<"/"<<CSCChBad[st]<<" is a dead chamber."<<endl;
#endif
	/*Not necessary. Chamber search using phi is accurate enough.
	CSCDetId Layerid  = CSCDetId( ec, st+1, rg,  CSCChCand[st], 1 );
	vector<Float_t> EdgeAndDistToGap( GetEdgeAndDistToGap(trackRef,Layerid) );//values: 1-edge;2-err of edge;3-disttogap;4-err of dist to gap
	if (EdgeAndDistToGap[0]>0.||EdgeAndDistToGap[0]<-9990.) {
	  cerr<<"ch"<<CSCChCand[st]<<", oldedge="<<EdgeAndDistToGap[0]<<"cm --> ";
	  //try neighborhood chambers
	  CSCDetId Layerid_plusone  = CSCDetId( ec, st+1, rg,  CSCChCand[st]+1, 1 );
	  vector<Float_t> EdgeAndDistToGap_plusone( GetEdgeAndDistToGap(trackRef,Layerid_plusone) );
	  CSCDetId Layerid_minusone  = CSCDetId( ec, st+1, rg,  CSCChCand[st]-1, 1 );
	  vector<Float_t> EdgeAndDistToGap_minusone( GetEdgeAndDistToGap(trackRef,Layerid_minusone) );
	  if (EdgeAndDistToGap_plusone[0]<EdgeAndDistToGap[0]&&EdgeAndDistToGap_plusone[0]>-9990.) {
	    if (EdgeAndDistToGap_minusone[0]>EdgeAndDistToGap_plusone[0]) CSCChCand[st]+=1;
	    else { if (EdgeAndDistToGap_minusone[0]>-9990.) CSCChCand[st]-=1;}
	  }
	  else { if (EdgeAndDistToGap_minusone[0]<EdgeAndDistToGap[0]&&EdgeAndDistToGap_minusone[0]>-9990.) CSCChCand[st]-=1;}
	  cerr<<"ch"<<CSCChCand[st]<<", plusoneedge="<<EdgeAndDistToGap_plusone[0]<<"cm,"<<", minusoneedge="<<EdgeAndDistToGap_minusone[0]<<"cm."<<endl;
	  }*/
	for (Int_t ly=1;ly<7;ly++) {
	  CSCDetId Layerid = CSCDetId( ec, st+1, rg,  CSCChCand[st], ly );
	  vector<Float_t> EdgeAndDistToGap( GetEdgeAndDistToGap(trackRef,Layerid) );//values: 1-edge;2-err of edge;3-disttogap;4-err of dist to gap
	  if (EdgeAndDistToGap[0]>CSCProjEdgeDist[st]) {
	    CSCProjEdgeDist[st]=EdgeAndDistToGap[0];
	    CSCProjEdgeDistErr[st]=EdgeAndDistToGap[1];
	  }
	  if (EdgeAndDistToGap[2]<CSCDyProjHVGap[st]) {
	    CSCDyProjHVGap[st]=EdgeAndDistToGap[2];
	    CSCDyErrProjHVGap[st]=EdgeAndDistToGap[3];
	  }
	}
	//cerr<<"To Edge:"<<CSCProjEdgeDist[st]<<"; To HVGap:"<<CSCDyProjHVGap[st]<<endl;
	// Check the CSC segments in that region..
	CSCSegmentCollection::const_iterator cscSegOut;
	TrajectoryStateOnSurface *TrajToSeg = matchTTwithCSCSeg( trackRef, cscSegments,  cscSegOut, Layer0Id);//update tsos with segment Z position, the segment z has little difference with layer 0 Z after the second or the third decimal place (in cm).
	if ( TrajToSeg!=NULL ) {
	  /*
	  const GeomDet* tmp_gdet=cscGeom->idToDet(cscSegOut->cscDetId());
	  const CSCChamber* cscchamber = cscGeom->chamber(cscSegOut->cscDetId());
	  cerr<<"CSCSEGDetID:"<<cscSegOut->cscDetId()<<" -- posZ: "<<tmp_gdet->surface().position().z()
	      <<" -- CSCSeg posZ: "<< cscchamber->toGlobal(cscSegOut->localPosition()).z()<<endl;
	  */
	  /* Save the chamber ID */
	  CSCDetId id  = (CSCDetId) cscSegOut->cscDetId();
	  /* Save the segment postion and direction */
	  LocalPoint localSegPos = (*cscSegOut).localPosition();
	  //GlobalPoint globalSegPosition = cscchamber->toGlobal( localSegPos );
	  CSCSegxLc[st] = localSegPos.x(); 
	  CSCSegyLc[st] = localSegPos.y(); 
	  LocalError localSegErr = (*cscSegOut).localPositionError();
	  CSCSegxErrLc[st] = sqrt(localSegErr.xx());
	  CSCSegyErrLc[st] = sqrt(localSegErr.yy());
	  //cerr<<"CSCSegxLc:"<<localSegPos.x()<<"+-"<<CSCSegxErrLc[st]<<endl;
	  //cerr<<"CSCSegyLc:"<<localSegPos.y()<<"+-"<<CSCSegyErrLc[st]<<endl;

	  /* Save the segment quality */
	  CSCnSegHits[st] = cscSegOut->specificRecHits().size();
	  Int_t nDOFCSC = 2*CSCnSegHits[st]-4;
	  CSCSegChisqProb[st] = ChiSquaredProbability( double( (*cscSegOut).chi2() ), nDOFCSC );
	    
	  /* Save the difference between the ex-tracker track and the segment */
	  const GeomDet* gdet=cscGeom->idToDet(id);
	  LocalPoint localpCSC = gdet->surface().toLocal(TrajToSeg->freeState()->position());
	  CSCDxTTSeg[st] = CSCSegxLc[st] - localpCSC.x(); 
	  CSCDyTTSeg[st] = CSCSegyLc[st] - localpCSC.y(); 
	  CSCDxyTTSeg[st] = sqrt(pow(CSCDxTTSeg[st],2)+pow(CSCDyTTSeg[st],2));

	  LocalError localTTErr =TrajToSeg->localError().positionError();
	  CSCDxErrTTSeg[st] = localSegErr.xx() + localTTErr.xx(); //the x error of the distance between segment and track squared
	  CSCDyErrTTSeg[st] = localSegErr.yy() + localTTErr.yy(); //the y error of the distance between segment and track squared
	  CSCDxyErrTTSeg[st] = sqrt(pow(CSCDxTTSeg[st],2)*CSCDxErrTTSeg[st]+pow(CSCDyTTSeg[st],2)*CSCDyErrTTSeg[st])/CSCDxyTTSeg[st];
	  CSCDxErrTTSeg[st] = sqrt(CSCDxErrTTSeg[st]); CSCDyErrTTSeg[st] = sqrt(CSCDyErrTTSeg[st]);       
	    
	  LocalVector trackLocalDir = TrajToSeg->localDirection();
	  LocalVector segDir = (*cscSegOut).localDirection();
	  if ( trackLocalDir.z()!=0. && segDir.z()!=0.) {
	    Float_t dxdz_trk = trackLocalDir.x()/trackLocalDir.z();
	    Float_t dydz_trk = trackLocalDir.y()/trackLocalDir.z(); 
	    Float_t dxdz_seg = segDir.x()/segDir.z();
	    Float_t dydz_seg = segDir.y()/segDir.z();
	    if(fabs(id.station()) == 3 || fabs(id.station()) == 4) dydz_seg = -dydz_seg;
	    CSCdXdZTTSeg[st] = dxdz_trk - dxdz_seg;
	    CSCdYdZTTSeg[st] = dydz_trk - dydz_seg;
	  }
	  nTrkCountCSCSeg++;
	  delete TrajToSeg;
	}

       	////// Loop over MPC infromation to look for LCTs....
	CSCDetId Layer3id  = CSCDetId( ec, st+1, rg,  CSCChCand[st], 3 );//go to layer 3 that corresponds to the LCTPos
	const GeomDet* Layer3gdet=cscGeom->idToDet(Layer3id);
	tsos=surfExtrapTrkSam(trackRef, Layer3gdet->surface().position().z());
	if (tsos.isValid()) {
	  LocalPoint localL3pCSC = Layer3gdet->surface().toLocal(tsos.freeState()->position());
	  LocalPoint *LCTPos=matchTTwithLCTs( localL3pCSC.x(), localL3pCSC.y(), CSCEndCapPlus?1:2, st+1, CSCRg[st], CSCChCand[st], mpclcts, CSCDxyTTLCT[st], CSCLCTbx[st]);
	  if (LCTPos!=NULL) {
	    CSCLCTxLc[st]=LCTPos->x();
	    CSCLCTyLc[st]=LCTPos->y();
	    CSCDxTTLCT[st]=CSCLCTxLc[st]-localL3pCSC.x();
	    CSCDyTTLCT[st]=CSCLCTyLc[st]-localL3pCSC.y();
	    LocalError localTTErr =tsos.localError().positionError();
	    CSCDxErrTTLCT[st] = sqrt(localTTErr.xx()); CSCDyErrTTLCT[st] = sqrt(localTTErr.yy());
	    CSCDxyErrTTLCT[st] =sqrt(pow(CSCDxTTLCT[st],2)*localTTErr.xx() + pow(CSCDyTTLCT[st],2)*localTTErr.yy())/CSCDxyTTLCT[st];
	    delete LCTPos;
	    //cerr<<"segZ_TTy,LCTZ_TTy:"<<CSCSegyLc[st]-CSCDyTTSeg[st]<<","<<localL3pCSC.y()<<";"<<endl;
 	  }//end of if found LCT
	}//end of if layer 3 extrapolation is available
      } // for loop for the stations -- j
      /*
      for (Byte_t st=0;st<4;st++) {
	lctSt[st] = ( CSCDxyTTLCT[st] >0. && CSCDxyTTLCT[st] < 40 )?1:0;
	segSt[st] = ( CSCDxyTTSeg[st] >0. && CSCDxyTTSeg[st] < 40)?1:0;
	}*/
      //cerr<<"eta="<<tracks_eta<<":ME1"<<CSCRg[0]<<":TTLCTx,TTLCTy,TTLCTxy:"<<CSCDxTTLCT[0]<<","<<CSCDyTTLCT[0]<<","<<CSCDxyTTLCT[0]<<endl;
      fractNtuple->Fill();
      firsttrackmatchingtoMuTag=false;
    } // loop over tracks...
  } // loop over muon tracks
}

// ------------ method called once each job just before starting event loop  ------------
void TPTrackMuonSys::beginJob() {}
// ------------ method called once each job just after ending the event loop  ------------
void TPTrackMuonSys::endJob() {
  theFile->Write();
  theFile->Close();
  std::cout<<"Events in "<<nEventsAnalyzed<<std::endl;
}

//copied and modified from http://www.codeproject.com/KB/string/wildcmp.aspx
//It does not belong to any class
Bool_t wildcmp(const char *wild, const char *string) {
  // Written by Jack Handy - jakkhandy@hotmail.com
  const char *cp = NULL, *mp = NULL;

  while ((*string) && (*wild != '*')) {
    if ((*wild != *string) && (*wild != '?')) {
      return false;
    }
    wild++;
    string++;
  }

  while (*string) {
    if (*wild == '*') {
      if (!*++wild) {
        return true;
      }
      mp = wild;
      cp = string+1;
    } else if ((*wild == *string) || (*wild == '?')) {
      wild++;
      string++;
    } else {
      wild = mp;
      string = cp++;
    }
  }

  while (*wild == '*') {
    wild++;
  }
  return !*wild;
}

// ------------ method called in the beginning of each run  ------------
void TPTrackMuonSys::beginRun(const Run& r, const EventSetup& iSet)
{
  run_number = r.runAuxiliary().run();
  iSet.get<CSCBadChambersRcd>().get(pBad);
  badChambers_=const_cast<CSCBadChambers*>(pBad.product());
  badChambersIndices=new vector<Int_t>( badChambers_->container() );

  char plotname[100],plottitle[100];
  sprintf(plotname,"Run%d_BadChambers",run_number);
  sprintf(plottitle,"Known BadChambers in Run %d; chamber number",run_number);
  /* Draw a badchamber plot and save it to the ntuple*/
  TH2F *TH2F_BadChambers=new TH2F(plotname,plottitle,36,1,37,18,-9,9);
  const char *chambers[36]  = {"01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"};
  const char *rings[18] = {"ME-42","ME-41","ME-32","ME-31","ME-22","ME-21","ME-13","ME-12","ME-11","ME+11","ME+12","ME+13","ME+21","ME+22","ME+31","ME+32","ME+41","ME+42"};
  for (Byte_t i=0;i<36;i++) 
    TH2F_BadChambers->GetXaxis()->SetBinLabel(i+1,chambers[i]);
  for (Byte_t i=0;i<18;i++)
    TH2F_BadChambers->GetYaxis()->SetBinLabel(i+1,rings[i]);
  for( Short_t indexc = 1; indexc<=540; ++indexc ) {// chamber indices are in range 1-468 (CSCs 2008) or 469-540 (ME42)
    CSCDetId id = CSCIndexer().detIdFromChamberIndex( indexc ); 
    if ( !badChambers_->isInBadChamber( id ) ) continue;
    Byte_t ring=id.station()*10+id.ring();
    Float_t fillY;
    switch ( ring )
      {
      case 14:
	fillY=0.5;
	break;
      case 11:
	fillY=0.5;
	break;
      case 12:
	fillY=1.5;
	break;
      case 13:
	fillY=2.5;
	break;
      case 21:
	fillY=3.5;
	break;
      case 22:
	fillY=4.5;
	break;
      case 31:
	fillY=5.5;
	break;
      case 32:
	fillY=6.5;
	break;
      case 41:
	fillY=7.5;
	break;
      case 42:
	fillY=8.5;
	break;
      default:
	printf("Unexpected ring number: %d",ring);
	fillY=9.5;
      }
    if (id.endcap()==2) fillY*=-1;
    TH2F_BadChambers->Fill(id.chamber()+0.5,fillY);
#ifdef jz_debug
    cerr<<(id.endcap()==1?"ME+":"ME-")<<Int_t(ring)<<"/"<<id.chamber()<<"("<<indexc<<")";
#endif
  }
//  TCanvas *BadChambersView=new TCanvas("badch","badch",1200,1000);
//  BadChambersView->SetGrid();
  TH2F_BadChambers->SetStats(0);
  TH2F_BadChambers->SetMinimum(0);
  TH2F_BadChambers->SetMaximum(0.8);
  TH2F_BadChambers->Draw("colz");
  TH2F_BadChambers->SetLabelSize(0.035,"X");
#ifndef GetCSCHitsBefore500
  theFile->cd();
#endif
  TH2F_BadChambers->Write();
  // BadChambersView->Write("BadChambersPlot");
  /* End of drawing the badchamber plot*/

  Bool_t isConfigChanged = false;
  m_HLTMuTrgBit.clear();
  HLTMuNames->clear();
  HLTMuObjModuleNames->clear();
  m_HLTDiMuTrgBit=-1;
  HLTDiMuName=new string();
  HLTDiMuObjModuleName=new string();
  if ( hltConfigProvider_.init( r, iSet, m_hlt.process() , isConfigChanged ) ) {
    HLTTableName=new string( hltConfigProvider_.tableName() );
    const vector<string> & HLTNamesSet_=hltConfigProvider_.triggerNames();
    UInt_t idx=0;
    for ( vector<string>::const_iterator itertable=HLTNamesSet_.begin();itertable != HLTNamesSet_.end();itertable++,idx++ ) {
#ifdef jz_debug
      cout<<endl<<idx<<":"<<*itertable;
#endif
      Bool_t selected=false;
      for ( vector<string>::const_iterator iter=m_HLTMuTrgNames.begin();iter!=m_HLTMuTrgNames.end();iter++ ) {
	if ( wildcmp(iter->c_str(),itertable->c_str()) ) {
	  selected=true;
	  break;
	}
      }
      if (selected) {
	m_HLTMuTrgBit.push_back(idx);
	HLTMuNames->push_back(*itertable);
	HLTMuObjModuleNames->push_back( *(hltConfigProvider_.moduleLabels(idx).end()-2) );
#ifdef jz_debug
	cout<<", saved for singlemu";
#endif
      }
      if ( m_HLTDiMuTrgBit<0 && wildcmp(m_HLTDiMuTrgName.c_str(),itertable->c_str()) ) {
	m_HLTDiMuTrgBit=idx;
	HLTDiMuName->assign(*itertable);
	HLTDiMuObjModuleName->assign( *(hltConfigProvider_.moduleLabels(idx).end()-2) );
#ifdef jz_debug
	cout<<", saved for doublemu";
#endif
      }
    }//end of iter HLT table
    RunInfo->Fill();
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    LogWarning("DataLost") << " HLT config extraction failure with process name " <<m_hlt.process();
    // In this case, all access methods will return empty values!
  }
}

vector<Float_t> TPTrackMuonSys::GetEdgeAndDistToGap(reco::TrackRef trackRef, CSCDetId & detid) {
  vector<Float_t> result(4,9999.);
  result[3]=-9999;
  const GeomDet* gdet=cscGeom->idToDet( detid );
  TrajectoryStateOnSurface tsos=surfExtrapTrkSam(trackRef, gdet->surface().position().z());
  if (!tsos.isValid()) return result;
  LocalPoint localTTPos = gdet->surface().toLocal(tsos.freeState()->position());
  const CSCWireTopology* wireTopology = cscGeom->layer(detid)->geometry()->wireTopology();
  Float_t wideWidth      = wireTopology->wideWidthOfPlane();
  Float_t narrowWidth    = wireTopology->narrowWidthOfPlane();
  Float_t length         = wireTopology->lengthOfPlane();
  // If slanted, there is no y offset between local origin and symmetry center of wire plane
  Float_t yOfFirstWire   = fabs(wireTopology->wireAngle())>1.E-06 ? -0.5*length : wireTopology->yOfWire(1);
  // y offset between local origin and symmetry center of wire plane
  Float_t yCOWPOffset    = yOfFirstWire+0.5*length;
  // tangent of the incline angle from inside the trapezoid
  Float_t tangent = (wideWidth-narrowWidth)/(2.*length);
  // y position wrt bottom of trapezoid
  Float_t yPrime  = localTTPos.y()+fabs(yOfFirstWire);
  // half trapezoid width at y' is 0.5 * narrowWidth + x side of triangle with the above tangent and side y'
  Float_t halfWidthAtYPrime = 0.5*narrowWidth+yPrime*tangent;
  // x offset between local origin and symmetry center of wire plane is zero
  // x offset of ME11s is also zero. x center of wire groups is not at zero, because it is not parallel to x. The wire groups of ME11s have a complex geometry, see the code in m_debug.
  Float_t edgex = fabs(localTTPos.x()) - halfWidthAtYPrime;
  Float_t edgey = fabs(localTTPos.y()-yCOWPOffset) - 0.5*length;
  LocalError localTTErr = tsos.localError().positionError();
  //  cerr<<"x err = "<<sqrt(localTTErr.xx())<<" ; y err = "<<sqrt(localTTErr.yy())<<endl;
  if ( edgex > edgey ) {
    result[0] = edgex;
    result[1] = sqrt( localTTErr.xx() );
    //result[1] = sqrt(tsos.cartesianError().position().cxx()); 
  }
  else {
    result[0] = edgey;
    result[1] = sqrt( localTTErr.yy() );
    //result[1] = sqrt(tsos.cartesianError().position().cyy());
  }
  result[2]=YDistToHVDeadZone(localTTPos.y(), detid.station()*10+detid.ring());
  result[3]=sqrt( localTTErr.yy() );
  return result;//return values: 1-edge;2-err of edge;3-disttogap;4-err of dist to gap
}

////////////////////////////////////////////
reco::MuonCollection::const_iterator TPTrackMuonSys::matchTTwithMT(reco::TrackCollection::const_iterator &itrack)
{
  for (reco::MuonCollection::const_iterator muIter = muons->begin(); muIter != muons->end(); muIter++) {
    if (muIter->combinedMuon().isNonnull() ) 
      if ( &( * muIter->combinedMuon().get() ) == & (*itrack) ) return muIter;
    if (muIter->track().isNonnull() )
      if ( &(* muIter->track().get() ) == & (*itrack) ) return muIter;
    if (muIter->standAloneMuon().isNonnull() ) 
      if ( &( *muIter->standAloneMuon().get() ) == & (*itrack) ) return muIter;
  }
  return muons->end();
}

////// Match TT with RPCEndCapHit
Bool_t TPTrackMuonSys::matchTTwithCSCRecHit(Bool_t trackDir,
					    Int_t j, 
					    reco::TrackRef trackRef, 
					    edm::Handle<CSCRecHit2DCollection> recHits, 
					    //std::vector<CSCRecHit2DCollection> recHitOut,
					    std::vector<CSCRecHit2D> & recHitOut, 
					    std::vector<Int_t >& deltaRecHitX,
					    std::vector<Int_t >& deltaRecHitY)
{

  Bool_t aMatch = false;
  CSCRecHit2D tmp[6];

  Float_t minDX[6]={999., 999., 999., 999., 999., 999.};
  Float_t minDY[6]={999., 999., 999., 999., 999., 999.};
  Float_t minDR[6]={999., 999., 999., 999., 999., 999.};
  recHitOut.clear();

  Float_t rCut = 50.;

  for (CSCRecHit2DCollection::const_iterator recIt = recHits->begin(); recIt != recHits->end(); recIt++) {
    CSCDetId id = (CSCDetId)(*recIt).cscDetId();
    
    if((trackDir && id.endcap() == 2) || (!trackDir && id.endcap() == 1)) continue;
    if(
       (j == 0 && id.station() == 1 && ((id.ring() == 1) || (id.ring() == 4 ))) ||  // ME1/1,4
       (j == 1 && id.station() == 1 && (id.ring() == 2)) ||  // ME1/2
       (j == 2 && id.station() == 1 && (id.ring() == 3)) ||  // ME1/3
       (j == 3 && id.station() == 2) ||                      // ME2/1,2
       (j == 4 && id.station() == 3) ||                      // ME3/1,2
       (j == 5 && id.station() == 4)                         // ME4/1
       ){
      
      // Get pointer to the layer:
      const CSCLayer* csclayer = cscGeom->layer( id );
      if (!csclayer)continue;
      
      // Transform hit position from local chamber geometry to global CMS geom
      LocalPoint rhitlocal = (*recIt).localPosition();  
      GlobalPoint rhitglobal = csclayer->toGlobal(rhitlocal);
      Float_t RecHit_grecx   = rhitglobal.x();
      Float_t RecHit_grecy   = rhitglobal.y();
      TrajectoryStateOnSurface tsos = surfExtrapTrkSam(trackRef, rhitglobal.z());  	       
      if (!tsos.isValid()) continue;
      Float_t cscdeltaX = RecHit_grecx - tsos.globalPosition().x();
      Float_t cscdeltaY = RecHit_grecy - tsos.globalPosition().y();
      Float_t cscdeltaR = sqrt(pow(cscdeltaX,2)+pow(cscdeltaY,2));

      Int_t mlayer = id.layer() - 1;
      if(fabs(cscdeltaR) < fabs(minDR[mlayer]) &&  fabs(cscdeltaR) < rCut){
	aMatch = true;
	minDX[mlayer] = cscdeltaX; 
	minDY[mlayer] = cscdeltaY;
	minDR[mlayer] = cscdeltaR;
	tmp[mlayer] = (*recIt);
      }
    }
  }

  for (Int_t i = 0; i < 6; i++) {
    Float_t cscdeltaR = sqrt(pow(minDX[i],2)+pow(minDY[i],2));
    if(cscdeltaR  < rCut) recHitOut.push_back(tmp[i]);
  }
  for (Int_t i = 0; i < 6; i++) {
    deltaRecHitX.push_back(minDX[i]);
    deltaRecHitY.push_back(minDY[i]);
  }
  return aMatch;
}

////// Match TT with RPCEndCapHit
/*
Bool_t TPTrackMuonSys::matchTTwithRPCEChit(Bool_t trackDir, 
					   Int_t j, 
					   reco::TrackRef trackRef, 
					   edm::Handle<RPCRecHitCollection> rpcRecHits, 
					   RPCRecHitCollection::const_iterator &rpcHitOut) 
//				       Float_t deltaRPCX,  Float_t deltaRPCY)
{
  Bool_t aMatch = false;

  //Float_t deltaRPCX = 999.0, deltaRPCY = 999.0; 
  Float_t deltaRPCR = 999.0;// deltaRPCPhi = 999.0;

  for (RPCRecHitCollection::const_iterator rpcIt = rpcRecHits->begin(); rpcIt != rpcRecHits->end(); rpcIt++) {
    RPCDetId id = (RPCDetId)(*rpcIt).rpcId();
    const GeomDet* gdet=rpcGeo->idToDet(id);
    const BoundPlane & surface = gdet->surface();
    if(id.region() == 0)continue;
    if((trackDir && id.region() == -1) || (!trackDir && id.region() == 1) ) continue;
    if(
       (j >= 1 && j <= 2 && id.station() == 1) ||
       (j == 3 && id.station() == 2) ||
       (j == 4 && id.station() == 3) 
       ){
      
      LocalPoint localPos = (*rpcIt).localPosition();
      GlobalPoint globalPosition=surface.toGlobal(localPos); 
      TrajectoryStateOnSurface tsos = surfExtrapTrkSam(trackRef, globalPosition.z());  
      if (!tsos.isValid()) continue;
      
      Float_t RPCEdeltaXMin = globalPosition.x() - tsos.globalPosition().x();
      Float_t RPCEdeltaYMin = globalPosition.y() - tsos.globalPosition().y();
      //Float_t RPCEdeltaRMin = sqrt(pow(RPCEdeltaXMin,2)+pow(RPCEdeltaYMin,2));

      Float_t RPCEPhiProj = tsos.globalPosition().phi();
      Float_t RPCEPhiPos  = globalPosition.phi();
      
      ////// New for RPCEC
      if(RPCEPhiProj  < 0) RPCEPhiProj =  RPCEPhiProj + 2*M_PI;
      if(RPCEPhiPos  < 0) RPCEPhiPos =  RPCEPhiPos + 2*M_PI;
      
      //      Float_t RPCEdeltaPhiMin = fabs(RPCEPhiProj - RPCEPhiPos);//sqrt(pow(RPCEdeltaXMin,2)+pow(RPCEdeltaYMin,2));

      Float_t RPCEdeltaRMin = sqrt(pow(RPCEdeltaXMin,2)+pow(RPCEdeltaYMin,2));
      if(fabs(RPCEdeltaRMin) < fabs(deltaRPCR) && fabs(RPCEdeltaRMin) < 50. ){
	aMatch = true;
	//	deltaRPCPhi = RPCEdeltaPhiMin;
	deltaRPCR = RPCEdeltaRMin;
	rpcHitOut = rpcIt;
      }
    }
  } // loop over RPCE rechits
   
  return aMatch;
}
*/
//////////////  Get the matching with LCTs...
LocalPoint * TPTrackMuonSys::matchTTwithLCTs(Float_t xPos, Float_t yPos, Short_t ec, Short_t st, Short_t &rg, Short_t cham, 
					     edm::Handle<CSCCorrelatedLCTDigiCollection> mpclcts, Float_t &dRTrkLCT, Int_t &lctBX ) {
  LocalPoint *interSect=NULL;
  
  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator detMPCUnitIt = mpclcts->begin(); 
       detMPCUnitIt != mpclcts->end(); detMPCUnitIt++) {
    CSCDetId id = (*detMPCUnitIt).first;
    
    if(ec != id.endcap())continue;
    if(st != id.station())continue;
    if(cham != id.chamber())continue;

    Bool_t ed1 = (st == 1) && ((rg == 1 || rg == 4) && (id.ring() == 1 || id.ring() == 4));
    Bool_t ed2 = (st == 1) && ((rg == 2 && id.ring() == 2) || (rg == 3 && id.ring() == 3));
    Bool_t ed3 = (st != 1) && (rg == id.ring());
    if ( !(ed1 || ed2 || ed3) ) continue;
    const CSCCorrelatedLCTDigiCollection::Range& MPCrange = (*detMPCUnitIt).second;
    for (CSCCorrelatedLCTDigiCollection::const_iterator mpcIt = MPCrange.first; mpcIt != MPCrange.second; mpcIt++) {
      Bool_t lct_valid = (*mpcIt).isValid();
      if(!lct_valid)continue;
      //In CSC offline/hlt software in general, is COUNT FROM ONE, such as CSCGeometry.
      //However, the LCT software counts from zero. strip_id is from 0 to 159. wireGroup_id is from 0 to 31(ME13),47(ME11),63(ME1234/2),95(ME31,41),111(ME21). So here we need to plus one.
      Byte_t wireGroup_id = (*mpcIt).getKeyWG()+1;
      Byte_t strip_id=(*mpcIt).getStrip()/2+1;
      //if ( id.ring()==4 ) cerr<<"Alert id.ring()==4"<<endl;
      Bool_t me11=(st == 1) && (id.ring() == 1 || id.ring() == 4),
	me11a = me11 && strip_id>64;
      //me11b = me11 && strip_id<=64;
      //http://cmssdt.cern.ch/SDT/lxr/source/CondFormats/CSCObjects/src/CSCChannelTranslator.cc
      //Translate a raw strip channel in range 1-80, iraw,  into corresponding geometry-oriented channel in which increasing
      //channel number <-> strip number increasing with +ve local x.
      //However, the LCT firmwire has already done the number flipping but it numbers ME11B strips from  0 to 63, ME11A strips from 64-79. (confirmed by Slava)
      //That is the reason we commented out the following number flipping codes.
      if ( me11a ) {
	strip_id-=64;
	// The CSCCorrelatedLCTDigi DetId does NOT distinguish ME11A and B. All of the DetIDs are labelled as ME11B (all ME11, none ME14)
	// However, stripWireGroupIntersection must know that since ME11A has 48 strips and ME11B has 64.
	id=CSCDetId(ec, 1, 4, cham, 3);
	//if ( id.endcap()==1 ) strip_id=17-strip_id;
      }
      //if ( me11b && id.endcap()!=1 ) strip_id=65-strip_id;
      const CSCLayerGeometry *layerGeom = cscGeom->chamber(id)->layer (3)->geometry ();
      for(Byte_t ii = 0; ii < 3; ii++){
	// if ( strip_id>64 ) LogWarning("Strip_id") << "Got "<<strip_id<<", but there are "<< Nstrips <<" strips in total." <<m_hlt.process();
	LocalPoint interSect_ = layerGeom->stripWireGroupIntersection(strip_id, wireGroup_id);
	//	printf( "ME%d/%d: %.1f/%d, %d/%d: xLCT-xTT=%.2f-%.2f; yLCT-yTT=%.2f-%.2f \n",st,id.ring(),strip_id,Nstrips,wireGroup_id,layerGeom->numberOfWireGroups(),interSect_.x(),xPos,interSect_.y(),yPos);
	Float_t DeltaR_ = sqrt(pow((interSect_.x()-xPos),2) + pow((interSect_.y()-yPos),2));
	if( DeltaR_ < fabs(dRTrkLCT) ) {
	  delete interSect;
	  interSect=new LocalPoint(interSect_);
	  dRTrkLCT =  DeltaR_ ;
	  lctBX = (*mpcIt).getBX();
	  if (me11a) rg=4;
	  else rg=id.ring();
	  //cout << "1: BX = " << (*mpcIt).getBX() << " BX0 = " << (*mpcIt).getBX0() << std::endl;
	} // for the matching if statement...
	if (me11a) strip_id+=16;
	else break;
      }// end iteration over of ME11A triplet
    }// end iteration over lcts_mpc_data in one chamber
  }// end iteration over lcts_mpc_data
  return interSect;
}

inline Float_t TPTrackMuonSys::TrajectoryDistToSeg( TrajectoryStateOnSurface *TrajSuf, CSCSegmentCollection::const_iterator segIt) {  
  if (!TrajSuf->isValid()) return 9999.;
  const GeomDet* gdet=cscGeom->idToDet( (CSCDetId)(*segIt).cscDetId() );
  LocalPoint localTTPos = gdet->surface().toLocal(TrajSuf->freeState()->position());
  LocalPoint localSegPos = (*segIt).localPosition();
  Float_t CSCdeltaX = localSegPos.x() - localTTPos.x();
  Float_t CSCdeltaY = localSegPos.y() - localTTPos.y();
  return sqrt(pow(CSCdeltaX,2)+pow(CSCdeltaY,2));
}

TrajectoryStateOnSurface *TPTrackMuonSys::matchTTwithCSCSeg( reco::TrackRef trackRef, edm::Handle<CSCSegmentCollection> cscSegments, 
							     CSCSegmentCollection::const_iterator &cscSegOut, CSCDetId & idCSC ) {
  TrajectoryStateOnSurface *TrajSuf=NULL;
  Float_t deltaCSCR = 9999.;
  for(CSCSegmentCollection::const_iterator segIt=cscSegments->begin(); segIt != cscSegments->end(); segIt++) {
    CSCDetId id  = (CSCDetId)(*segIt).cscDetId();
    if(idCSC.endcap() != id.endcap())continue;
    if(idCSC.station() != id.station())continue;
    if(idCSC.chamber() != id.chamber())continue;
      
    Bool_t ed1 = (idCSC.station() == 1) && ((idCSC.ring() == 1 || idCSC.ring() == 4) && (id.ring() == 1 || id.ring() == 4));
    Bool_t ed2 = (idCSC.station() == 1) && ((idCSC.ring() == 2 && id.ring() == 2) || (idCSC.ring() == 3 && id.ring() == 3));
    Bool_t ed3 = (idCSC.station() != 1) && (idCSC.ring() == id.ring());
    Bool_t TMCSCMatch = (ed1 || ed2 || ed3);
    if(! TMCSCMatch)continue;
    
    const CSCChamber* cscchamber = cscGeom->chamber(id);
    if (!cscchamber) continue;
    TrajectoryStateOnSurface TrajSuf_ = surfExtrapTrkSam(trackRef, cscchamber->toGlobal( (*segIt).localPosition() ).z());

    Float_t dR_= fabs( TrajectoryDistToSeg( &TrajSuf_, segIt ) );
    if ( dR_ < deltaCSCR ){
      delete TrajSuf;
      TrajSuf=new TrajectoryStateOnSurface(TrajSuf_);
      deltaCSCR = dR_;
      cscSegOut = segIt;
    }
  }//loop over segments

  return TrajSuf;
}
//////////////  Get the matching with CSC-sgements...
Bool_t TPTrackMuonSys::matchTTwithCSCSeg(Bool_t trackDir, Int_t st_, reco::TrackRef trackRef, edm::Handle<CSCSegmentCollection> cscSegments, 
					 CSCSegmentCollection::const_iterator &cscSegOut ) {
  Bool_t aMatch = false;
  //  Float_t deltaCSCX = 9999.0, deltaCSCY = 9999.0;
  Float_t deltaCSCR = 50.0;
  for(CSCSegmentCollection::const_iterator segIt=cscSegments->begin(); segIt != cscSegments->end(); segIt++) {
    CSCDetId id  = (CSCDetId)(*segIt).cscDetId();
    if((trackDir && id.endcap() == 2) || (!trackDir && id.endcap() == 1)) continue;
    if(
       (st_ == 0 && id.station() == 1 && ((id.ring() == 1) || (id.ring() == 4 ))) ||  // ME1/1,4
       (st_ == 1 && id.station() == 1 && (id.ring() == 2)) ||  // ME1/2
       (st_ == 2 && id.station() == 1 && (id.ring() == 3)) ||  // ME1/3
       (st_ == 3 && id.station() == 2) ||                      // ME2/1,2
       (st_ == 4 && id.station() == 3) ||                      // ME3/1,2
       (st_ == 5 && id.station() == 4)                         // ME4/1
       ) {
      const CSCChamber* cscchamber = cscGeom->chamber(id);
      if (!cscchamber) continue;
      TrajectoryStateOnSurface TrajSuf = surfExtrapTrkSam(trackRef, cscchamber->toGlobal( (*segIt).localPosition() ).z());
      
      Float_t dR_= fabs( TrajectoryDistToSeg( &TrajSuf, segIt ) );
      
      if ( dR_ < deltaCSCR ){
	aMatch = true;
	deltaCSCR = dR_;
	cscSegOut = segIt;
      }    
    } 
  } // loop over segments
  return aMatch;
}

////
Int_t TPTrackMuonSys::getNLayerMatchedCSCSeg(CSCSegmentCollection::const_iterator &cscSegMatch,
					     edm::Handle<CSCRecHit2DCollection> recHits,
					     Float_t *delRecSegX,
					     Float_t *delRecSegY,
					     Int_t &nGhits
					     )
{
  Int_t nhits = 0;
  CSCDetId idCSC  = (CSCDetId)(*cscSegMatch).cscDetId();
  Int_t ec = idCSC.endcap(), st = idCSC.station(), rg = idCSC.ring(), cham = idCSC.chamber();    
  LocalPoint localPos = (*cscSegMatch).localPosition();

  Float_t cscx = localPos.x();
  Float_t cscy = localPos.y();

  for(Int_t i = 0; i < 6; i++){
    delRecSegX[i]  = 999.;
    delRecSegY[i]  = 999.;
  }
  Float_t minDeltaR = 999.;

  for (CSCRecHit2DCollection::const_iterator recIt = recHits->begin(); recIt != recHits->end(); recIt++) {
    CSCDetId id = (CSCDetId)(*recIt).cscDetId();
    if(ec   != id.endcap())continue;
    if(st   != id.station())continue;
    if(rg   != id.ring())continue;
    if(cham != id.chamber())continue;
    LocalPoint rhitlocal = (*recIt).localPosition();  
    Float_t recx   = rhitlocal.x();
    Float_t recy   = rhitlocal.y();
    Float_t myDeltaR = sqrt(pow((cscx - recx),2)+ pow((cscy - recy),2));

    Int_t layer = id.layer();

    if(myDeltaR < minDeltaR){
      minDeltaR = myDeltaR;
      delRecSegX[layer-1] = cscx - recx;
      delRecSegY[layer-1] = cscy - recy;
    }
    nhits = nhits + 1;
  }

  nGhits = 0;
  for(Int_t i = 0; i < 6; i++){
    Float_t myDeltaR = sqrt(pow(delRecSegX[i],2)+ pow(delRecSegY[i],2));
    if(myDeltaR < 15.)nGhits = nGhits + 1; /// look for this hit is within 15 cm.
  }


  return nhits; // this is probably not the right number of hits.
}




#ifdef GetCSCHitsBefore500
void TPTrackMuonSys::getCSCSegWkeyHalfStrip(const std::vector<CSCRecHit2D> &theseRecHits, Float_t &cStrp, Float_t &ckWG){
  
#ifdef m_debug
  cout << " Start of the function " << endl;
#endif

  Float_t m_cStrp = 999.0, m_ckWG = 999.;

  Float_t strpSum =0.0, wkeySum = 0.0;
  Int_t nhits =0 ;
  Bool_t layer3 = false;
  Bool_t me11a = false;

  for ( std::vector<CSCRecHit2D>::const_iterator iRH = theseRecHits.begin(); iRH != theseRecHits.end(); iRH++) {
    CSCDetId idRH = (CSCDetId)(*iRH).cscDetId();
    Int_t mx_layer   = idRH.layer();

    if (idRH.ring() == 4) me11a = true;

    /////////////////////////    
    //
    // Get the strip number... ------> m_cStrp 
    //
    /////////////////////////    
    LocalPoint rhitlocal = (*iRH).localPosition();  
    Float_t pWS  = (*iRH).positionWithinStrip();
    // Find the strip containing this hit
    CSCRecHit2D::ChannelContainer hitstrips = (*iRH).channels();

    Int_t nStrips   =  hitstrips.size();
    //  Int_t centerid  =  nStrips/2 + 1;    

    if(nStrips == 3) m_cStrp  = (int) (2.0*(hitstrips[1] + pWS - 0.5 ));
    else m_cStrp  = (int) ( 2.0*(hitstrips[0] - pWS - 0.5 ));

    Bool_t evenLayer =  (mx_layer % 2 == 0);
    if ( evenLayer )m_cStrp -= 1;
    if ( (idRH.station() == 1) && (idRH.layer() != 3) && evenLayer ) m_cStrp += 1;
    //m_cStrp    =  hitstrips[centerid - 1];
    /////////////////////////    
    //
    // Get the wireKey number...------> m_ckWG
    //
    /////////////////////////    

    CSCRecHit2D::ChannelContainer wiresg = (*iRH).wgroups();
    //  m_ckWG = wiresg[0]-1;
    m_ckWG = wiresg[wiresg.size()/2]-1;//corrected from to wiresg[0] to [wiresg.size()/2]

    if(mx_layer == 3){
      cStrp = m_cStrp;
      ckWG = m_ckWG;
      layer3 = true;
    }
    strpSum = strpSum + m_cStrp;
    wkeySum = wkeySum + m_ckWG;
    nhits = nhits + 1;
  }
  if(! layer3 && nhits != 0){
    cStrp = strpSum/nhits;
    ckWG  = wkeySum/nhits;
  }
  if (me11a)cStrp = ((int)cStrp - 1)%32 + 1;
}
#else
void TPTrackMuonSys::getCSCSegWkeyHalfStrip(const vector<CSCRecHit2D> &theseRecHits, Float_t &cStrp, Float_t &ckWG){
  Float_t strpSum = 0., wkeySum = 0.;
  Int_t nhits = 0;  Bool_t me11a=false;
  for ( vector<CSCRecHit2D>::const_iterator iRH = theseRecHits.begin(); iRH != theseRecHits.end(); iRH++,nhits++ ) {
    CSCDetId idRH = (CSCDetId)(*iRH).cscDetId();
    if (idRH.ring() == 4) me11a = true;
    Byte_t layer = idRH.layer();
    Float_t pWS  = (*iRH).positionWithinStrip();
    Int_t m_cStrp  = (Int_t) (2.0*(iRH->channels(1) + pWS - 0.5 ));
    if ( layer%2!=0 &&
	 !( (idRH.station() == 1) && (idRH.layer() != 3) ) ) m_cStrp--;

    if( layer == 3 ) {
      cStrp=m_cStrp;
      ckWG=iRH->hitWire()-1;
      nhits=1;
      break;
    }
    else {
      strpSum += m_cStrp;
      wkeySum += iRH->hitWire()-1;
    }
#ifdef m_debug
    clog << "ME" << idRH.station() << idRH.ring() <<":" << endl
	 << "cstrips ("<< iRH->nStrips() <<"): " //should be always 3
	 << iRH->channels(0) << "," << iRH->channels(1) << "," << iRH->channels(2) 
	 << "; cwire:" << iRH->hitWire() << endl;
#endif
  }
  cStrp = strpSum/nhits;
  ckWG  = wkeySum/nhits;
  if (me11a) cStrp = ( (Int_t) cStrp - 1)%32 + 1;
#ifdef m_debug
  clog << "cStrip: "<< cStrp << "; cWire: " << ckWG << endl;
#endif
}
#endif// end of GetCSCHitsBefore500
/*
Bool_t TPTrackMuonSys::matchCSCSegWithLCT(edm::Handle<CSCCorrelatedLCTDigiCollection> mpclcts, 
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
					  )
{
  /// Input 
  // idCSC, halfstrip and wireKey

  /// Output 
  /// Whether there is a match based on halfstrip and wireKey
  /// or based on phi.
  /// also: minPhi, minDAngle;

  for(Int_t i = 0; i < 2; i++){
    xMatch[i] = false;
    mDAngle[i]  = 999.;
    diffTrkEta[i] = 999.;
    diffTrkPhi[i] = 999.;
    delHStrp[i] = 99999.;
    delWkey[i] = 99999.;
  }

  Int_t ec = idCSC.endcap(), st = idCSC.station(), rg = idCSC.ring(), cham = idCSC.chamber();
  Bool_t me11a = false;
  if (rg == 4) me11a = true;

  //// For the minimization:

  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator detMPCUnitIt = mpclcts->begin(); detMPCUnitIt != mpclcts->end(); detMPCUnitIt++) {
    const CSCDetId& id = (*detMPCUnitIt).first;
	     
    if(ec != id.endcap())continue;
    if(st != id.station())continue;
    if(cham != id.chamber())continue;

    Bool_t ed1 = (st == 1) && ((rg == 1 || rg == 4) && (id.ring() == 1 || id.ring() == 4));
    Bool_t ed2 = (st == 1) && ((rg == 2 && id.ring() == 2) || (rg == 3 && id.ring() == 3));
    Bool_t ed3 = (st != 1) && (rg == id.ring());
    Bool_t TMCSCMatch = (ed1 || ed2 || ed3);
    if(! TMCSCMatch)continue;

    const CSCCorrelatedLCTDigiCollection::Range& MPCrange = (*detMPCUnitIt).second;
    for (CSCCorrelatedLCTDigiCollection::const_iterator mpcIt = MPCrange.first; mpcIt != MPCrange.second; mpcIt++) {
      Bool_t lct_valid = (*mpcIt).isValid();
      if(!lct_valid)continue;
     
      Int_t station = id.station()-1;
      Int_t sector  = id.triggerSector()-1;
      Int_t subSector = CSCTriggerNumbering::triggerSubSectorFromLabels(id);
      Int_t fpga    = ( subSector ? subSector-1 : station+1 );
      Int_t endcap = id.endcap()-1;
      lclphidat lclPhi = srLUTs_[endcap][sector][fpga]->localPhi((*mpcIt).getStrip(),(*mpcIt).getPattern(),(*mpcIt).getQuality(),(*mpcIt).getBend() );
      gblphidat gblPhi = srLUTs_[endcap][sector][fpga]->globalPhiME( lclPhi.phi_local, (*mpcIt).getKeyWG(), id.triggerCscId() );
      gbletadat gblEta = srLUTs_[endcap][sector][fpga]->globalEtaME(lclPhi.phi_bend_local, lclPhi.phi_local,(*mpcIt).getKeyWG(), id.triggerCscId() );
      
      UInt_t mpcphi = gblPhi.global_phi;
      UInt_t mpceta = gblEta.global_eta;
      double geta = theTriggerScales->getRegionalEtaScale(2)->getCenter(mpceta); //Type 2 is CSC
      double radphi = ((mpcphi)/4096.)*(62./180.)*M_PI;    // cscphi in rad
      
      double gphi   = radphi + sector*M_PI/3.+ (14*M_PI/180.);//0.24434609;// Global Phi with respect to 0 degrees.
      if(sector == 5){
	if(radphi <=45*(M_PI/180.)) gphi = radphi + 315*M_PI/180.;
	else gphi   = radphi - 45*M_PI/180.;
      }
      if(TrkPhi < 0)TrkPhi = TrkPhi + 2*M_PI;
      if(gphi < 0)gphi = gphi + 2*M_PI;
      Float_t DeltaAngle = sqrt(pow((TrkEta-geta),2) + pow((TrkPhi-gphi),2));

      Int_t m_strip  = (*mpcIt).getStrip(); // halfstrip that goes from 0 to 31 
      if (me11a)m_strip = (m_strip-1)%32+1;
      Float_t m_keyWG  = (*mpcIt).getKeyWG(); // 
      Float_t m_delHStrp = c1 - m_strip;
      Float_t m_delWkey  = w1 - m_keyWG;

      if (me11a){ // the ganging of ME11a causes wraparound effects at the boundaries for delta strip 
	if (m_delHStrp > 16) m_delHStrp -= 32;
	if (m_delHStrp < -16)m_delHStrp += 32;
      }

      Float_t phiDiff = TrkPhi-gphi;
      if( fabs(phiDiff) < fabs(diffTrkPhi[0]) ){
	xMatch[0]     = true;
	mDAngle[0]    = DeltaAngle;
	diffTrkEta[0] = TrkEta-geta; 
	diffTrkPhi[0] = phiDiff;
	delHStrp[0]   = m_delHStrp;
	delWkey[0]    = m_delWkey;
	mpcItOut      = mpcIt;
      } // for the matching if statement...
      
      if( (fabs(m_delHStrp) < fabs(delHStrp[1])) && (fabs(m_delWkey) < fabs(delWkey[1])) ){ // match strips at <=10 and wirekeys at <=5
	xMatch[1]     = true;
	mDAngle[1]    = DeltaAngle;
	diffTrkEta[1] = TrkEta-geta;
	diffTrkPhi[1] = TrkPhi-gphi;
	delHStrp[1]   = m_delHStrp;
	delWkey[1]    = m_delWkey;
	mpcHsWkOut    = mpcIt;
      } // for the matching if statement...

#ifdef m_debug
      std::cout << "MPC E: " << id.endcap() << " R:" << id.ring() << " S: " << id.station() << " C: " << id.chamber() 
		<< std::endl;
#endif
    }
  }// end iteration over lcts_mpc_data
  return (xMatch[0] || xMatch[1]);
} 
*/
void TPTrackMuonSys::chamberCandidates(Int_t station, Float_t feta, Float_t phi, std::vector <int> &coupleOfChambers){
  // yeah, hardcoded again...
  coupleOfChambers.clear();

  Int_t ring = ringCandidate(station, feta, phi);
  if(ring != -9999){
    Float_t phi_zero = 0.;// check! the phi at the "edge" of Ch 1
    Float_t phi_const = 2.*M_PI/36.;
    Int_t first_chamber = 1, last_chamber = 36;
    if(1 != station && 1==ring){ // 18 chambers in the ring
      phi_const*=2;
      last_chamber /= 2;
    }
    if(phi < 0.) phi += 2*M_PI;
    
    Float_t chamber_float = (phi - phi_zero)/phi_const;
    Int_t chamber_int = int(chamber_float);
    if (chamber_float - Float_t(chamber_int) -0.5 <0.){

      if(0!=chamber_int ) coupleOfChambers.push_back(chamber_int);
      else coupleOfChambers.push_back(last_chamber);
      coupleOfChambers.push_back(chamber_int+1); 

    }else{
      coupleOfChambers.push_back(chamber_int+1);
      if(last_chamber!=chamber_int+1) coupleOfChambers.push_back(chamber_int+2);
      else coupleOfChambers.push_back(first_chamber);
    }
  }
}


Int_t TPTrackMuonSys::ringCandidate(Int_t station, Float_t feta, Float_t phi){
  switch (station){
  case 1:
    if(fabs(feta)>=0.85 && fabs(feta)<1.18){//ME13
      return 3;
    }      
    if(fabs(feta)>=1.18 && fabs(feta)<=1.5){//ME12 if(fabs(feta)>1.18 && fabs(feta)<1.7){//ME12
      return 2;
    }
    if(fabs(feta)>1.5 && fabs(feta)<2.1){//ME11
      return 1;
    }
    if(fabs(feta)>=2.1 && fabs(feta)<2.45){//ME11
      return 4;
    }
    break;
  case 2:
    if(fabs(feta)>0.95 && fabs(feta)<1.6){//ME22
      return 2;      
    }  
    if(fabs(feta)>1.55 && fabs(feta)<2.45){//ME21
      return 1;      
    }
    break;
  case 3:
    if(fabs(feta)>1.08 && fabs(feta)<1.72){//ME32
      return 2;            
    }  
    if(fabs(feta)>1.69 && fabs(feta)<2.45){//ME31
      return 1; 
    }
    break;
  case 4:
    if(fabs(feta)>1.78 && fabs(feta) <2.45){//ME41
      return 1; 
    }
    if(fabs(feta)>1.15 && fabs(feta) <=1.78){//ME42
      return 2; 
    }
    break;
  default:
    LogError("")<<"Invalid station: "<<station<<endl;
    break;
  }
  return -9999;
}

Short_t TPTrackMuonSys::thisChamberCandidate(Short_t station, Short_t ring, Float_t phi){
  //search for chamber candidate based on CMS IN-2007/024
  //10 deg chambers are ME1,ME22,ME32,ME42 chambers; 20 deg chambers are ME21,31,41 chambers
  //Chambers one always starts from approx -5 deg.
  const Short_t nVal = (station>1 && ring == 1)?18:36;
  const Float_t ChamberSpan=2*M_PI/nVal;
  Float_t dphi = phi + M_PI/36;
  while (dphi >= 2*M_PI) dphi -= 2*M_PI;
  while (dphi < 0) dphi += 2*M_PI;
  Short_t ChCand=floor(dphi/ChamberSpan)+1;
  return ChCand>nVal?nVal:ChCand;
}

///////////////////////
// to get the track position info at a particular rho
TrajectoryStateOnSurface TPTrackMuonSys::cylExtrapTrkSam(reco::TrackRef track, double rho)
{
  Cylinder::PositionType pos(0, 0, 0);
  Cylinder::RotationType rot;
  Cylinder::CylinderPointer myCylinder = Cylinder::build(pos, rot, rho);
  try{
    FreeTrajectoryState recoStart = freeTrajStateMuon(track);
    TrajectoryStateOnSurface recoProp = propagatorAlong->propagate(recoStart, *myCylinder);
    if (!recoProp.isValid()) {
      recoProp = propagatorOpposite->propagate(recoStart, *myCylinder);
    }
    return recoProp;
  }
  catch(cms::Exception){
    edm::LogError("")<<"invalid track extrapolation to cylinder"<<endl;
    TrajectoryStateOnSurface recoProp;
    return recoProp;
  }
}

// to get track position at a particular (xy) plane given its z
TrajectoryStateOnSurface TPTrackMuonSys::surfExtrapTrkSam(reco::TrackRef track, double z)
{
  Plane::PositionType pos(0, 0, z);
  Plane::RotationType rot;
  Plane::PlanePointer myPlane = Plane::build(pos, rot);
  try{
    FreeTrajectoryState recoStart = freeTrajStateMuon(track);
    TrajectoryStateOnSurface recoProp = propagatorAlong->propagate(recoStart, *myPlane);
    if (!recoProp.isValid()) {
      recoProp = propagatorOpposite->propagate(recoStart, *myPlane);
    }
    return recoProp;
  }
  catch(cms::Exception){
    edm::LogError("")<<"invalid track extrapolation to plane"<<endl;
    TrajectoryStateOnSurface recoProp;
    return recoProp;
  }
}

FreeTrajectoryState TPTrackMuonSys::freeTrajStateMuon(reco::TrackRef track)
{
  GlobalPoint  innerPoint(track->innerPosition().x(), track->innerPosition().y(),  track->innerPosition().z());
  GlobalVector innerVec  (track->innerMomentum().x(), track->innerMomentum().y(),  track->innerMomentum().z());  
  
  GlobalTrajectoryParameters gtPars(innerPoint, innerVec, track->charge(), &*theBField);

  //FreeTrajectoryState recoStart(innerPoint, innerVec, track->charge(), &*theBField);
  //return recoStart;

  AlgebraicSymMatrix66 cov;
  cov *= 1e-20;

  CartesianTrajectoryError tCov(cov);
 
  return (cov.kRows == 6 ? FreeTrajectoryState(gtPars, tCov) : FreeTrajectoryState(gtPars)) ;
}


Bool_t TPTrackMuonSys::GetDecayChains(TrackingParticleRef tpr, HepMC::GenEvent *HepGenEvent, ULong64_t &truth_type, Int_t & truth_thesamewith, vector<vector< vector<Int_t> > > & ExistingSimChains) {
  Bool_t IsPileup=false;
  truth_thesamewith=-1;
  vector<TheTrackType> types;
  for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )  {
    //g4Track loop begin
    DChain.clear();SimChains.clear();
    const SimTrack *thisTrk=&(*g4Track_iter);
    SimTrackDaughtersTree( thisTrk,tpr );
    Bool_t ChainEnd; TrackingParticleRef tpr_tmp=tpr;
    do {
      ChainEnd=true;
      if ( !thisTrk->noVertex() ) {
	TrackingVertexRef tvr=tpr_tmp->parentVertex();
	for ( TrackingParticleRefVector::iterator parenttp=tvr->sourceTracks_begin();parenttp!=tvr->sourceTracks_end();parenttp++ ) {
	  for (vector<SimTrack>::const_iterator g4Trk_iter = (*parenttp)->g4Track_begin() ; g4Trk_iter != (*parenttp)->g4Track_end(); ++g4Trk_iter )
	    if ( SVC[thisTrk->vertIndex()].parentIndex()==Int_t(g4Trk_iter->trackId()) && g4Trk_iter->eventId().rawId()==thisTrk->eventId().rawId()) { 
	      thisTrk=&(*g4Trk_iter);tpr_tmp=*parenttp;
	      vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
	      if (SavedSimTrk_iter==SavedSimTrk.end()) {
		DChain.push_back( MCParticlesList.size() );
		SavedSimTrk.push_back(thisTrk);
		MCParticlesList.push_back( MCParticleInfo_Creator(thisTrk,tpr_tmp) );
	      }
	      else DChain.push_back( FindSimTrackInMCParticlesList(SavedSimTrk_iter-SavedSimTrk.begin()));
	      ChainEnd=false;
	      break;
	    }
	  if (!ChainEnd) break;
	}
      }
    }while (!ChainEnd);
    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
      SimChains_iter->insert(SimChains_iter->begin(),DChain.rbegin(),DChain.rend());
    DChain.clear();
    
    //HepMC Particles
    HepMCChains.clear();
    if (!thisTrk->noGenpart() ) {
      HepMC::GenParticle *genPar=HepGenEvent->barcode_to_particle(thisTrk->genpartIndex());
      if (genPar!=NULL) {
      	vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
	if (SavedSimTrk_iter!=SavedSimTrk.end()) MCParticlesList[FindSimTrackInMCParticlesList(SavedSimTrk_iter-SavedSimTrk.begin())].IsParticleFromGenerator=true;
	else LogError("CodeWrong")<<"Cannot find the simulated track in saved sim tracks.";
	HepMCParentTree( genPar );
      }
      else LogWarning("RefNull")<<"Either SimTrack::genpartIndex() or HepMC::GenParticle::barcode is wrong. Pilup track?";
    }

    //merge the HepMC and SimTrack Decay Chains
    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
      if ( !HepMCChains.empty() )
	for (vector< vector<Int_t> >::iterator HepMCChains_iter = HepMCChains.begin(); HepMCChains_iter !=  HepMCChains.end(); ++HepMCChains_iter ) {
	  vector<Int_t> thisChain(HepMCChains_iter->rbegin(),HepMCChains_iter->rend());
	  thisChain.insert(thisChain.end(),SimChains_iter->begin(),SimChains_iter->end());
	  IsPileup=SaveAndClassify(thisChain,types, truth_thesamewith, ExistingSimChains);
	}
      else IsPileup=SaveAndClassify( *SimChains_iter,types, truth_thesamewith, ExistingSimChains );
  }//g4Track loop end
  for (vector <TheTrackType>::iterator type_iter=types.begin();type_iter!=types.end();type_iter++)
    truth_type=truth_type*100+Long64_t(*type_iter);
  return IsPileup;
}

void TPTrackMuonSys::SimTrackDaughtersTree(const SimTrack * thisTrk, TrackingParticleRef tpr)
{
  // Find MC Truth Segment - the offical one use chi2 to match simtrk (MuonIdTruthInfo.cc) and it won't know the decay in flight segment truth
  // The particle type of the hit may differ from the particle type of the SimTrack with id trackId().
  // This happends if the hit was created by a secondary track(e.g. a delta ray) originating from the trackId() and not existing as aseparate SimTrack.
  // ( particle type match notice is from haiyun.teng@cern.ch )
  Bool_t ChainEnd=true;
  //To avoid duplicate particle saving
  vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
  if (SavedSimTrk_iter==SavedSimTrk.end()) {
    DChain.push_back( MCParticlesList.size() );
    SavedSimTrk.push_back(thisTrk);	
    MCParticlesList.push_back( MCParticleInfo_Creator(thisTrk,tpr) );
  }
  else DChain.push_back( FindSimTrackInMCParticlesList(SavedSimTrk_iter-SavedSimTrk.begin()) );

  for ( TrackingVertexRefVector::iterator tvr=tpr->decayVertices().begin();tvr!=tpr->decayVertices().end();tvr++ )
    for ( TrackingParticleRefVector::iterator daughtertp=(*tvr)->daughterTracks_begin();daughtertp!=(*tvr)->daughterTracks_end();daughtertp++ )
      for (vector<SimTrack>::const_iterator g4Trk_iter = (*daughtertp)->g4Track_begin() ; g4Trk_iter != (*daughtertp)->g4Track_end(); ++g4Trk_iter )
	if ( SVC[g4Trk_iter->vertIndex()].parentIndex()==Int_t(thisTrk->trackId()) && g4Trk_iter->eventId().rawId()==thisTrk->eventId().rawId()) { 
	  ChainEnd=false;
	  SimTrackDaughtersTree( &(*g4Trk_iter), *daughtertp );
	}
  if (ChainEnd) SimChains.push_back(DChain);
  DChain.pop_back();
}

void TPTrackMuonSys::HepMCParentTree(HepMC::GenParticle *genPar) {
  HepMC::GenVertex *thisVtx = genPar->production_vertex();
  Bool_t ChainEnd=true;
  if (thisVtx) {
    for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
      if ((*pgenD)->pdg_id()!=92)  {//Pythia special code for string, we only care about the particles after hadronization
	ChainEnd=false;
	vector<HepMC::GenParticle *>::iterator SavedHepPar_iter=find(SavedHepPar.begin(),SavedHepPar.end(),*pgenD);
	if (SavedHepPar_iter==SavedHepPar.end())
	  {
	    DChain.push_back(MCParticlesList.size());
	    SavedHepPar.push_back(*pgenD);
	    MCParticlesList.push_back( MCParticleInfo_Creator( (*pgenD) ) );
	  }
	else DChain.push_back(FindHepMCInMCParticlesList(SavedHepPar_iter-SavedHepPar.begin()));
	HepMCParentTree(*pgenD);
	DChain.pop_back();
      }
  }
  if (ChainEnd) HepMCChains.push_back(DChain);
}

Bool_t TPTrackMuonSys::SaveAndClassify(vector<Int_t> &Chain, vector<TheTrackType> &types, Int_t &truth_thesamewith, vector< vector< vector<Int_t> > > &ExistingSimChains) {
  //Find out if another track is just a decay product or mother particle of this track
  truth_thesamewith=-1;
  Bool_t truth_isPileup=false;
  for ( vector< vector< vector<Int_t> > >::const_iterator atrack=ExistingSimChains.begin(); atrack!=ExistingSimChains.end(); atrack++ ) {
    for ( vector< vector<Int_t> >::const_iterator adecaychain=atrack->begin(); adecaychain!=atrack->end(); adecaychain++ )
      if ( IstheSameDChain(Chain,*adecaychain) ) {
	truth_thesamewith=atrack-ExistingSimChains.begin();
	break;
      }
    if ( truth_thesamewith>-1 ) break;
  }

  //Save
  ExistingSimChains.back().push_back(Chain);
  
  //Classification
  TheTrackType type=Classify(Chain);
  if ( !MCParticlesList[Chain[0]].IsParticleFromGenerator ) truth_isPileup=true;
#ifdef jz_debug
  //  if ( Int_t(type)<2||Int_t(type)>100|| truth_isPileup ) {
  cerr<<"======================================================="<<endl;
  for (vector<Int_t>::reverse_iterator iter = Chain.rbegin(); iter != Chain.rend(); iter++)
    cerr<<"("<<MCParticlesList[*iter].pdgId<<","<<MCParticlesList[*iter].DoesParticleHaveMuonHit<<") <--|";
  cerr<<endl;
  cerr<<"MuTagpdgID:"<<MuTagtracktruth_id<<"; pdgID:"<<tracktruth_id;
  cerr<<"; type is "<<Int_t(type)<<"; pileup: "<<Int_t(truth_isPileup)<<endl
      <<"======================================================="<<endl;
  //  }
#endif
  if ( types.empty() ) types.push_back(type);
  else if (type!=Others) { 
    if ( find(types.begin(),types.end(),type)==types.end() ) {
      vector<TheTrackType>::iterator DeterminedAsOthers=find(types.begin(),types.end(),Others);
      if ( DeterminedAsOthers!=types.end() ) *DeterminedAsOthers=type;
      else {
	vector<TheTrackType>::iterator DeterminedAsNoMuSysHit=find(types.begin(),types.end(),Others);
	if ( DeterminedAsNoMuSysHit!=types.end() ) *DeterminedAsOthers=type;
	else types.push_back(type);
      }
    }
  }
  return truth_isPileup;
}

inline TPTrackMuonSys::ParticleType TPTrackMuonSys::ParticleCata(Int_t pid)
{
  pid=pid>0?pid:pid*-1;
  if (pid==13) return Muon;
  if (pid==443) return JPsi;
  if (pid==23) return Z;
  if (pid==24) return W;
  if ((pid>110&&pid<900)||(pid>10000))
    {
      pid=pid%1000;
      if (pid<400) return LightMeson;
      if (pid>400&&pid<440) return CharmedMeson;
      if (pid>440&&pid<500) return ccbarMeson;
      if (pid>500&&pid<550) return BottomMeson;
      if (pid>550&&pid<600) return bbarMeson;
    }
  if (pid>1111&&pid<3350&&(pid%100)/10>0) return LightBaryon;
  if (pid>4111&&pid<4450&&(pid%100)/10>0) return CharmedBaryon;
  if (pid>5111&&pid<5560&&(pid%100)/10>0) return BottomBaryon;
  if (pid>1000&&pid<6000&&(pid%100)/10==0) return DiQuarks;
  if (pid>10&&pid<19) return Lepton;
  return Other;
}

TPTrackMuonSys::TheTrackType TPTrackMuonSys::Classify(vector<Int_t> &Chain) {
  Int_t MuPos=-1;
  for (vector<Int_t>::reverse_iterator iter = Chain.rbegin(); iter != Chain.rend(); iter++) {
    ParticleType ParType=ParticleCata(MCParticlesList[*iter].pdgId);
    if (ParType==Muon&&MuPos<0) {//only care of the last muon
      MuPos=*iter;
      if ( !MCParticlesList[MuPos].DoesParticleHaveMuonHit ) return NoMuSysHit;
      continue;
    }
    if (MuPos>=0&&(ParType==W||ParType==Z)) {
      if ( MCParticlesList[MuPos].IsParticleFromGenerator ) {
	if (ParType==W) return PromptMuFromW;
	else return PromptMuFromZ;
      }
      else return NotPromptMufromWZ;//should never happen
    }
    if (MuPos>=0&&(ParType==JPsi)) return PromptMuFromJPsi;
    if (MuPos<0&&ParType>=LightMeson&&ParType<=DiQuarks) {
      if ( MCParticlesList[*iter].DoesParticleHaveMuonHit ) return PunchThrough;
      else return NoMuSysHit;
    }
    if (MuPos>=0&&ParType>=LightMeson&&ParType<=BottomBaryon) {
      if ( MCParticlesList[MuPos].IsParticleFromGenerator||MCParticlesList[MuPos].IsParticleBornInsideOfARegion ) {//inside pixel detector
	if (ParType==LightMeson) return PromptMuFromLightMeson;
	if (ParType>LightMeson&&ParType<LightBaryon) return PromptMuFromHeavyMeson;
	if (ParType==LightBaryon) return PromptMuFromLightBaryon;
	if (ParType>LightBaryon&&ParType<=BottomBaryon) return PromptMuFromHeavyBaryon;
      }
      if ( MCParticlesList[*iter].DoesParticleHaveMuonHit ) return PunchThroughAndDecayInFlight;
      else {
	if (ParType==LightMeson) return DecayInFlightFromLightMeson;
	if (ParType>LightMeson&&ParType<LightBaryon) return DecayInFlightFromHeavyMeson;
	if (ParType==LightBaryon) return DecayInFlightFromLightBaryon;
	if (ParType>LightBaryon&&ParType<=BottomBaryon) return DecayInFlightFromHeavyBaryon;
      }
    }
  }
  if (MuPos>=0) return PromptMuFromOthers;
  else return Others;
}

Bool_t TPTrackMuonSys::IstheSameDChain(const vector<Int_t> &Chain1,const vector<Int_t> &Chain2) {//whether two chains include each other
  Bool_t ChainIncluded=false;
  vector<Int_t>::const_iterator Chain1_Particle = Chain1.begin(),Chain2_Particle = Chain2.begin();
  for (; Chain2_Particle != Chain2.end()&&Chain1_Particle != Chain1.end(); Chain2_Particle++) {
    if (ChainIncluded&&*Chain1_Particle!=*Chain2_Particle) ChainIncluded=false;
    if (Chain1.front() == *Chain2_Particle) {
      ChainIncluded=true;
      Chain1_Particle = Chain1.begin();
    }
    if (ChainIncluded) Chain1_Particle++;
  }
  if (!ChainIncluded) {
    Chain2_Particle = Chain2.begin(); Chain1_Particle = Chain1.begin();
    for (; Chain2_Particle != Chain2.end()&&Chain1_Particle != Chain1.end(); Chain1_Particle++) {
      if (ChainIncluded&&*Chain1_Particle!=*Chain2_Particle) ChainIncluded=false;
      if (Chain2.front() == *Chain1_Particle) {
	Chain2_Particle = Chain2.begin();
	ChainIncluded=true;
      }
      if (ChainIncluded) Chain2_Particle++;
    }
  }
  return ChainIncluded;
}

TPTrackMuonSys::MCParticleInfo TPTrackMuonSys::MCParticleInfo_Creator( const SimTrack * thisTrk, TrackingParticleRef tpr ) {
  MCParticleInfo TBA;
  TBA.IsThisFromSimTrk=true;
  TBA.IsParticleFromGenerator=false;TBA.DoesParticleHaveMuonHit=false;
  TBA.IsParticleBornInsideOfARegion=false;
  TBA.IsPileup=false;
  TBA.pdgId=thisTrk->type();
  for ( vector<PSimHit>::const_iterator g4Hit_iter=tpr->pSimHit_begin();g4Hit_iter!=tpr->pSimHit_end();g4Hit_iter++ )
    if ( g4Hit_iter->trackId()==thisTrk->trackId() && g4Hit_iter->eventId().rawId()==thisTrk->eventId().rawId() && g4Hit_iter->particleType() == thisTrk->type() ) {
      DetId DetectorId( g4Hit_iter->detUnitId() );
      if ( DetectorId.det() == DetId::Muon ) {
	TBA.DoesParticleHaveMuonHit=true;
	break;
      }
    }
  if (!thisTrk->noVertex()) {
    SimVertex thisVtx=SVC[thisTrk->vertIndex()];
    //if ( (thisVtx.position().Perp2()<161604.&&deltaZ<568.)||(deltaRPhi2>82024.96&&deltaRPhi2<161604.&&deltaZ<666.)) //inside HCAL
    if ( thisVtx.position().Perp2()<100&&
	 fabs(thisVtx.position().z())<30 ) //inside pixel detector
      TBA.IsParticleBornInsideOfARegion=true;
  }
  else TBA.IsPileup=true;
  return TBA;
}

//deadzone center is according to http://cmssdt.cern.ch/SDT/lxr/source/RecoLocalMuon/CSCEfficiency/src/CSCEfficiency.cc#605
//wire spacing is according to CSCTDR
Float_t TPTrackMuonSys::YDistToHVDeadZone(Float_t yLocal, Int_t StationAndRing){
  //the ME11 wires are not parallel to x, but no gap
  //chamber edges are not included.
  const Float_t deadZoneCenterME1_2[2] = {-32.88305,32.867423};
  const Float_t deadZoneCenterME1_3[2] = {-22.7401,27.86665};
  const Float_t deadZoneCenterME2_1[2] = {-27.47,33.67};
  const Float_t deadZoneCenterME3_1[2] = {-36.21,23.68};
  const Float_t deadZoneCenterME4_1[2] = {-26.14,23.85};
  const Float_t deadZoneCenterME234_2[4] = {-81.8744,-21.18165,39.51105,100.2939};
  const Float_t *deadZoneCenter;
  Float_t deadZoneHeightHalf=0.32*7/2;// wire spacing * (wires missing + 1)/2
  Float_t minY=999999.;
  Byte_t nGaps=2;
  switch (abs(StationAndRing)) {
  case 11:
  case 14:
    return 162;//the height of ME11
    break;
  case 12:
    deadZoneCenter=deadZoneCenterME1_2;
    break;
  case 13:
    deadZoneCenter=deadZoneCenterME1_3;
    break;
  case 21:
    deadZoneCenter=deadZoneCenterME2_1;
    break;
  case 31:
    deadZoneCenter=deadZoneCenterME3_1;
    break;
  case 41:
    deadZoneCenter=deadZoneCenterME4_1;
    break;
  default:
    deadZoneCenter=deadZoneCenterME234_2;
    nGaps=4;
  }
  for ( Byte_t iGap=0;iGap<nGaps;iGap++ ) {
    Float_t newMinY=yLocal<deadZoneCenter[iGap]?deadZoneCenter[iGap]-deadZoneHeightHalf-yLocal:yLocal-(deadZoneCenter[iGap]+deadZoneHeightHalf);
    if ( newMinY<minY ) minY=newMinY;
  }
  return minY;
} 

//define this as a plug-in
DEFINE_FWK_MODULE(TPTrackMuonSys);


