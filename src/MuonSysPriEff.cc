// -*- C++ -*-
//
// Package:    MuonSysPriEff
// 
/**\class MuonSysPriEff MuonSysPriEff.cc UserCode/MuonSysPriEff/src/MuonSysPriEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhang Jinzhong
//         Created:  Mon Jun  7 22:19:50 CEST 2010
//

#include "../interface/MuonSysPriEff.h"
#define Muon_Mass 0.105658367

#ifdef TrackingParticles
#define trackIdLink(TrkID)  vector<SimTrack>::iterator Trk_iter = STC.begin(); \
  for (; Trk_iter != STC.end(); ++Trk_iter )	\
    if ((int) Trk_iter->trackId() == (int) TrkID) break;		\
  if (Trk_iter==STC.end()) LogWarning("RefNULL")<<"parentIndex/trackId Error"; \
  thisTrk=&*Trk_iter;

#define GenSimMomentum(ppt,peta,pphi,ppdg) Gen_pdgId->push_back(ppdg);Gen_pt->push_back(ppt);Gen_eta->push_back(peta);Gen_phi->push_back(pphi)

#define GenSimVertex(vtx,vty,vtz,vtt) Gen_vx->push_back(vtx);Gen_vy->push_back(vty);Gen_vz->push_back(vtz);Gen_vt->push_back(vtt)

#define RecordSimTrack(thisTrk) SavedSimTrk.push_back(thisTrk);			\
  IsParInHep->push_back(false);						\
  GenSimMomentum(thisTrk->momentum().pt(),thisTrk->momentum().eta(),thisTrk->momentum().phi(),thisTrk->type());	\
  if (!thisTrk->noVertex()) {						\
    SimVertex thisVtx=SVC[thisTrk->vertIndex()];			\
    GenSimVertex(thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9); \
  }									\
  else {GenSimVertex(0,0,0,0);}

#define RecordHepMC(GenParticle) SavedHepPar.push_back(GenParticle);	\
  IsParInHep->push_back(true);						\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}
#else
#define RecordHepMC(GenParticle) SavedHepPar.push_back(GenParticle);	\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}
#endif

// ------------ method called on each new Event  ------------

void MuonSysPriEff::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  NumberOfEvents++;
//####in the first event, check the special event content which we may want####
  //IsMC?HasDigi? based on event content
  //CSC LCT Digis
  edm::Handle<CSCCorrelatedLCTDigiCollection> CSCLCTs;
  if (FirstEvent) {
    //Is it Data or MC?
    Handle<edm::HepMCProduct> HepMCH;
    event.getByLabel(HepMCTag, HepMCH);
    if (HepMCH.failedToGet()){
      LogInfo("DataFormat")<<"It is real data. HepMC collection is disabled";
      if ( !event.isRealData() ) LogWarning("DataLost")<<"Corrupted MC samples! HepMC collection is missing";
      IsMC=false;
    }
    //Does it has Muon Digis information?(reconstructed from RAW?)
    try{
      event.getByLabel(CSCDigisTag,CSCLCTs);
    }
    catch (cms::Exception){
      LogWarning("DataFormat")<< "Can't get CSC-LCTdigi by label. ";
      HasDigi=false;
    }
    //Does it has TrackingParticle collection? (RAWDEBUG,RECODEBUG or tracking particle building process has been taken)
    if (IsMC) {
      try{
	Handle<TrackingParticleCollection> TPTest;
	event.getByType(TPTest);
      }catch (cms::Exception){
	LogWarning("DataFormat")<< "No tracking particle collection is found.";
	HasTrackingParticle=false;
      }
    }
    FirstEvent=false;
  }
  //########## preselect muons according to optimized pt() #####################
  Handle<reco::MuonCollection> Muon;
  event.getByLabel("muons", Muon);
  typedef pair<reco::MuonCollection::const_iterator,double> MuonPreSelection;
  vector<MuonPreSelection> MuonQueue;
  reco::MuonCollection const & muons = *Muon;
  Handle <reco::TrackToTrackMap> tevMapH1,tevMapH2,tevMapH3;
  event.getByLabel("tevMuons", "default", tevMapH1);
  const reco::TrackToTrackMap tevMap1 = *(tevMapH1.product());
  event.getByLabel("tevMuons", "firstHit", tevMapH2);
  const reco::TrackToTrackMap tevMap2 = *(tevMapH2.product());
  event.getByLabel("tevMuons", "picky", tevMapH3);
  const reco::TrackToTrackMap tevMap3 = *(tevMapH3.product());
  Bool_t GoodEvent=true;
  if (Muon.failedToGet()) GoodEvent=false;
  else {
    for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter) {
      Float_t pt,eta;
      if ( iter->globalTrack().isNonnull() ) {
	TrackRef tmptrack=muon::tevOptimized(*iter,tevMap1,tevMap2,tevMap3);
	pt=tmptrack->pt();
	eta=tmptrack->eta();
      }
      else {
	pt=iter->pt();
	eta=iter->eta();
      }
      if (pt>MuonPtCut&&fabs(eta)<MuonEtaMax&&fabs(eta)>MuonEtaMin) MuonQueue.push_back(make_pair(iter,pt));
    }
    if (MuonQueue.size()<MinNumberOfMuons) GoodEvent=false;
  }
  if (ThrowBadEvents&&!GoodEvent) return;
//############### standard filters ############################
  isHLTTriggerred=HLTFilter->TriggerFilter(event,iSetup);
  isGoodVertex=VertexFilter->filter(event,iSetup);
  isNoScrapping=ScrapingFilter->filter(event);
  if (ThrowBadEvents&&( (IsMC&&!isHLTTriggerred) ||!isGoodVertex||!isNoScrapping)) return;//We keep all triggerred data with good vertex and no scrapping
//########## Generation Information #####################
 
  HepMC::GenEvent *HepGenEvent=NULL;
  if (IsMC) {
    ClearVecs_HepMC();
    Handle<edm::HepMCProduct> HepMCH;
    event.getByLabel(HepMCTag, HepMCH);
    HepGenEvent=new HepMC::GenEvent(*(HepMCH->GetEvent()));
    Handle<GenEventInfoProduct> hEvtInfo;
    event.getByLabel("generator", hEvtInfo);
    GenEventWeight = hEvtInfo->weight();
  }
  else GenEventWeight = 1.;

//############### sort Muon according to Tuned pt ############
  ClearVecs_RECO();
  if (MuonQueue.size()>1)
    for (vector<MuonPreSelection>::iterator iter_queue=MuonQueue.begin();iter_queue!=MuonQueue.end()-1; ++iter_queue )
      for (vector<MuonPreSelection>::iterator iter_queue2=iter_queue+1;iter_queue2!=MuonQueue.end(); ++iter_queue2 )
	if (iter_queue->second<iter_queue2->second) swap(*iter_queue,*iter_queue2);

//############### event information ##########################
  Info.RUN   = event.id ().run ();
  Info.EVENT = event.id ().event();
  Info.LS    = event.luminosityBlock ();
  Info.ORBIT = event.orbitNumber ();
  Info.BX = event.bunchCrossing ();
  isRealData = event.isRealData();

  //Beamspot
  Handle<reco::BeamSpot> beamSpotHandle;
  if (!event.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) LogWarning("DataLost")<<"BeamSpot Information is not found.";

  //Primary Vertex
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(PrimaryVerticesTag.c_str(),recVtxs);
  if (recVtxs.isValid())
    for(reco::VertexCollection::const_iterator v=recVtxs->begin(); v!=recVtxs->end(); ++v) {
      vx->push_back(v->x());	   vxError->push_back(v->xError());
      vy->push_back(v->y());	   vyError->push_back(v->yError());
      vz->push_back(v->z());	   vzError->push_back(v->zError());
    }
  else LogWarning("InvalidTag")<<PrimaryVerticesTag.c_str()<<" information is not valid.";
  
  // this is needed by the IPTools methods from the tracking group
  ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);
  //Pileup
  if (IsMC) {
    edm::InputTag PileupSrc_("addPileupInfo");
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    if ( !event.getByLabel(PileupSrc_, PupInfo) ) LogWarning("DataLost")<<"Pileup Information is not found.";
    else {
      numberOfPUVertices=0;
      for(std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
	numberOfPUVertices+=PVI->getPU_NumInteractions();
    }
  }

  //HLT
  Handle<TriggerResults> hltTriggerResults;
  event.getByLabel(TriggerResultsTag,hltTriggerResults);
  if (hltTriggerResults.isValid())
    for (UInt_t i=0; i < HLTSize; i++)
      HLTacceptance->push_back(hltTriggerResults->accept(i));
  else LogWarning("InvalidTag")<<"TriggerResultsTag: "<<TriggerResultsTag.label()<<"is not valid";

  //HLT Objects
  Handle<trigger::TriggerEvent> trgEvent;
  event.getByLabel(TriggerEventTag,trgEvent);
  if (trgEvent.isValid()) {
    const trigger::TriggerObjectCollection& TOC = trgEvent->getObjects();
    int Total_Filters = trgEvent->sizeFilters();
    for( unsigned int i =0; i < num_HLTsSaveObjs;i++) 
      if (HLTFilterNamesAcceptenceIndex[i]>=0) {
	trigger::size_type pos_Filter = trgEvent->filterIndex(HLT_ModuleNames[i]);
	if (pos_Filter < Total_Filters) {
	  const trigger::Keys& KEYS(trgEvent->filterKeys(pos_Filter));
	  unsigned int num_thisHLTObjs=KEYS.size();
	  for(unsigned int j = 0; j < num_thisHLTObjs; j++) {
	    const trigger::TriggerObject& TO = TOC[KEYS[j]];
	    //printf("pt%f\teta%f\tphi%f\n",TO.pt(),TO.eta(),TO.phi());
	    HLTObj_pt[i]->push_back(TO.pt());
	    HLTObj_eta[i]->push_back(TO.eta());
	    HLTObj_phi[i]->push_back(TO.phi());
	  }
	}
	else if ((*HLTacceptance)[HLTFilterNamesAcceptenceIndex[i]])
	  LogWarning("DataLost")<<(*HLTNamesSet)[HLTFilterNamesAcceptenceIndex[i]]<<" passed but "<<HLT_ModuleNames[i].label()<<" is not found.";
      }
  }
  else LogWarning("InvalidTag")<<"TriggerEventTag "<<TriggerEventTag.label()<<" is not valid.";

//################# Tracks Information #####################
#ifdef TrackingParticles
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(tracksTag,trackCollectionHV);
  //Simulated Vertices: the vertexId() is just the index
  Handle<SimVertexContainer> SVCollectionH;
  //TrackingParticles
  Handle<TrackingParticleCollection> TPCollectionH ;
  RecoToSimCollection RecoToSimByHits,RecoToSimByChi2;
  if (IsMC&&HasTrackingParticle) {
    //Simulated Vertices
    event.getByLabel("g4SimHits", SVCollectionH);
    SVC = *SVCollectionH.product();
    //Simulated Tracks: the trackId() is NOT the index
    Handle<SimTrackContainer> STCollectionH;
    event.getByLabel("g4SimHits", STCollectionH);
    SimTrackContainer STC = *STCollectionH.product();
    //find daughters of each SimTrack
    Daughters.clear();
    for (vector<SimTrack>::iterator FindTrack = STC.begin(); FindTrack != STC.end(); ++FindTrack )
      if (!FindTrack->noVertex()) {
	SimVertex thisVtx=SVC[FindTrack->vertIndex()];
	if (!thisVtx.noParent()) {
	  SimTrack *thisTrk;
	  trackIdLink(thisVtx.parentIndex())
	    Daughters[thisTrk].push_back(&*FindTrack);
	}
      } 
    //TrackingParticles
    event.getByType(TPCollectionH);
    //SimToReco Tracks Association
    ESHandle<TrackAssociatorBase> AssociatorByHits;
    iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", AssociatorByHits);
    RecoToSimByHits = AssociatorByHits->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
    //Match by chi2
    ESHandle<TrackAssociatorBase> AssociatorByChi2;
    iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2", AssociatorByChi2);
    RecoToSimByChi2 = AssociatorByChi2->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
  }
#endif
  //CSC Geometry
  ESHandle<CSCGeometry> cscGeometry;
  if (HasDigi) {
    iSetup.get<MuonGeometryRecord>().get(cscGeometry);
    event.getByLabel(CSCDigisTag,CSCLCTs);
  }
  //DEDX
  Handle< ValueMap<reco::DeDxData> > dEdxTrackHandle;
  event.getByLabel(dEdxTag, dEdxTrackHandle);
  const ValueMap<reco::DeDxData> dEdxTrack = *(dEdxTrackHandle.product());
//################# Main Muon LOOP #########################
  for (vector<MuonPreSelection>::iterator iter_queue=MuonQueue.begin();iter_queue!=MuonQueue.end(); ++iter_queue ) {//muon loop begin
    reco::MuonCollection::const_iterator iter=iter_queue->first;
    //muon basic information (pt,eta,phi,charge)
    //Global Muon uses Cocktail track selection, Track Muon uses tracker track
    TrackRef theMuonTrack;
    TrackRef InnerTrack = iter->innerTrack();
    if ( iter->isGlobalMuon() ) {
      theMuonTrack = muon::tevOptimized(iter->globalTrack(), InnerTrack, tevMap1, tevMap2, tevMap3);//change the global track to the optimized track (TuneP algo)
      isTrackerMu->push_back(iter->isTrackerMuon());
      isGlobalMu->push_back(true);
      numberOfValidMuonHits->push_back( theMuonTrack->hitPattern().numberOfValidMuonHits() );
      GlobalTrack_chi2->push_back(theMuonTrack->chi2());
      GlobalTrack_ndof->push_back(theMuonTrack->ndof());
    }
    else {
      if ( iter->isTrackerMuon() ) {
	theMuonTrack=InnerTrack;
	isTrackerMu->push_back(true);
	isGlobalMu->push_back(false);
	numberOfValidMuonHits->push_back(0);
	GlobalTrack_chi2->push_back(0);
	GlobalTrack_ndof->push_back(0);
      }
      else continue;
    }
    pt->push_back(theMuonTrack->pt());
    ptError->push_back(theMuonTrack->ptError());
    eta->push_back(theMuonTrack->eta());
    etaError->push_back(theMuonTrack->etaError());
    phi->push_back(theMuonTrack->phi());
    phiError->push_back(theMuonTrack->phiError());
    if (iter->charge()==-1) chargeMinus->push_back(true);
    else chargeMinus->push_back(false);
    //Muon Selectors
    for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
      SelectorPassed[whichcut]->push_back(muon::isGoodMuon(*iter,OfficialMuonSelectors[whichcut]));
    //isHLTObj
    for (unsigned int i=0; i<num_HLTsSaveObjs; i++)
      isHLTObj[i]->push_back(false);
    if ( trgEvent.isValid() )
      for (unsigned int i=0; i<num_HLTsSaveObjs; i++) {
	unsigned int num_thisHLTObjs=HLTObj_pt[i]->size();
	for (unsigned int j=0; j<num_thisHLTObjs; j++) {
	  float dR=(eta->back()-(*HLTObj_eta[i])[j])*(eta->back()-(*HLTObj_eta[i])[j])+(phi->back()-(*HLTObj_phi[i])[j])*(phi->back()-(*HLTObj_phi[i])[j]);
	  if (dR<MaxHLTObjDeviation2) {
	    isHLTObj[i]->back()=true;
	    break;
	  }
	}
      }
    //Muon Vertex
    Vertex_x->push_back(iter->vertex().X());
    Vertex_y->push_back(iter->vertex().Y());
    Vertex_z->push_back(iter->vertex().Z());
    //numberofmatches
    numberOfMatchedSegments->push_back( iter->numberOfMatches() );
    numberOfMatchedStations->push_back( iter->numberOfMatchedStations() );
    //Isolation
    isoR03sumPt->push_back(iter->isolationR03().sumPt);
    isoR05sumPt->push_back(iter->isolationR05().sumPt);
    //Kinks
    TrkKink->push_back(iter->combinedQuality().trkKink);
    GlbKink->push_back(iter->combinedQuality().glbKink);
    //Combined track quality
    TrkRelChi2->push_back(iter->combinedQuality().trkRelChi2);
   //match to RecoTrk
    dEdx->push_back(dEdxTrack[InnerTrack].dEdx());
    dEdxError->push_back(dEdxTrack[InnerTrack].dEdxError());
    dEdx_numberOfSaturatedMeasurements->push_back(dEdxTrack[InnerTrack].numberOfSaturatedMeasurements());
    dEdx_numberOfMeasurements->push_back(dEdxTrack[InnerTrack].numberOfMeasurements());
    InnerTrack_nValidTrackerHits->push_back(InnerTrack->hitPattern().numberOfValidTrackerHits());
    InnerTrack_nValidPixelHits->push_back(InnerTrack->hitPattern().numberOfValidPixelHits());
    InnerTrack_nLostPixelHits->push_back(InnerTrack->hitPattern().numberOfLostPixelHits());
    InnerTrack_nLostTrackerHits->push_back(InnerTrack->hitPattern().numberOfLostTrackerHits());
    InnerTrack_chi2->push_back(InnerTrack->chi2());
    InnerTrack_ndof->push_back(InnerTrack->ndof());
    //DXY using beamspot
    DXYwtBS->push_back( InnerTrack->dxy(beamSpotHandle->position()) );
    DZwtBS->push_back( InnerTrack->dz(beamSpotHandle->position()) );
    //DXY using first PV
    std::pair<bool,Measurement1D> result = IPTools::absoluteTransverseImpactParameter(trackBuilder->build(InnerTrack),recVtxs->at(0));
    DXYwtPV->push_back( result.second.value() );
    DXYErrwtPV->push_back( result.second.error() );

//#################Muon Chambers#######################
    vector<const MuonChamberMatch *> Chambers;
    for( vector<MuonChamberMatch>::const_iterator chamberMatch = iter->matches().begin();chamberMatch != iter->matches().end(); chamberMatch++ )
      if ( chamberMatch->detector()!=MuonSubdetId::RPC ) Chambers.push_back( &(*chamberMatch) );
    //align save position from DT to CSC, from - to +, from inner to outter,from smaller chamber # to big
    if (Chambers.size()>1)
      for( vector<const MuonChamberMatch *>::iterator chamberMatch1 = Chambers.begin();chamberMatch1 != Chambers.end()-1; chamberMatch1++ ) {
	Bool_t IsCSC1=(*chamberMatch1)->detector()==MuonSubdetId::CSC;
	Int_t wheelstation1=0,stationring1=0,sectorchamber1=0;
	if ( IsCSC1 ) {
	  const CSCDetId ChamberID( (*chamberMatch1)->id.rawId() );
	  wheelstation1=ChamberID.endcap()==1?-ChamberID.station():ChamberID.station();
	  stationring1=ChamberID.ring()%4;
	  sectorchamber1=ChamberID.chamber();
	}
	else {
	  const DTChamberId ChamberID( (*chamberMatch1)->id.rawId() );
	  wheelstation1=ChamberID.station();
	  stationring1=ChamberID.wheel();
	  sectorchamber1=ChamberID.sector();
	}
	for( vector<const MuonChamberMatch *>::iterator chamberMatch2 = chamberMatch1+1;chamberMatch2 != Chambers.end(); chamberMatch2++ ) {
	  Bool_t IsCSC2=(*chamberMatch2)->detector()==MuonSubdetId::CSC;
	  if ( IsCSC2&&! IsCSC1 ) continue;
	  Int_t wheelstation2=0,stationring2=0,sectorchamber2=0;
	  if ( IsCSC2 ) {
	    const CSCDetId ChamberID( (*chamberMatch2)->id.rawId() );
	    wheelstation2=ChamberID.endcap()==1?-ChamberID.station():ChamberID.station();
	    stationring2=ChamberID.ring()%4;
	    sectorchamber2=ChamberID.chamber();
	  }
	  else {
	    const DTChamberId ChamberID( (*chamberMatch2)->id.rawId() );
	    wheelstation2=ChamberID.station();
	    stationring2=ChamberID.wheel();
	    sectorchamber2=ChamberID.sector();
	  }
	  if ( IsCSC2==IsCSC1 ) {
	    if ( wheelstation2<wheelstation1 ) continue;
	    if ( wheelstation2==wheelstation1 ) {
	      if ( stationring2>stationring1 ) continue;
	      if ( ( stationring2==stationring1 )&&( sectorchamber2>sectorchamber1 ) ) continue;
	    }
	  }
	  swap(*chamberMatch1,*chamberMatch2);
	}
      }
    //save the chambers and the segments in the chambers information 
    for( vector<const MuonChamberMatch *>::iterator chamberMatch = Chambers.begin();chamberMatch != Chambers.end(); chamberMatch++ ) {//chamberloop begin
      std::vector<MuonSegmentMatch>::const_iterator segmentMatch = (*chamberMatch)->segmentMatches.begin();
      for( ;segmentMatch != (*chamberMatch)->segmentMatches.end(); segmentMatch++ ) 
	if (segmentMatch->isMask(MuonSegmentMatch::BestInStationByDR)) break;
      if ( segmentMatch == (*chamberMatch)->segmentMatches.end() &&  (*chamberMatch)->dist()==999999 ) continue;
      Bool_t IsCSC=(*chamberMatch)->detector()==MuonSubdetId::CSC;
      //Chamber Information
      Byte_t NumberOfLCTs=0; Float_t XthisLCT=999999.,YthisLCT=999999.;
      IsCSCChamber->push_back(IsCSC);
      if ( IsCSC ) {
	const CSCDetId ChamberID( (*chamberMatch)->id.rawId() );
	StationRing->push_back( ChamberID.endcap()==1?ChamberID.station()*10+ChamberID.ring():-ChamberID.station()*10-ChamberID.ring() );
	SectorChamber->push_back( ChamberID.chamber() );
	dYToDeadZone->push_back( YDistToHVDeadZone( (*chamberMatch)->y,StationRing->back() ) );
	if (HasDigi) {
	  Float_t dR2=999999.;
	  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator CSCDigi_iter = CSCLCTs->begin(); CSCDigi_iter != CSCLCTs->end(); CSCDigi_iter++) {
	    const CSCDetId& LCTDetID = (*CSCDigi_iter).first;
	    if ( ChamberID.endcap() != LCTDetID.endcap() || ChamberID.station() != LCTDetID.station() || ChamberID.chamber() != LCTDetID.chamber() || ChamberID.ring() != LCTDetID.ring() ) continue;
	    const CSCLayer *csclayer=cscGeometry->chamber(LCTDetID)->layer(3);
	    const CSCLayerGeometry *layerGeometry=csclayer->geometry();
	    const CSCCorrelatedLCTDigiCollection::Range& LCTRange = (*CSCDigi_iter).second;
	    for (CSCCorrelatedLCTDigiCollection::const_iterator LCT_iter = LCTRange.first; LCT_iter != LCTRange.second; LCT_iter++) {
	      NumberOfLCTs++;
	      Int_t WGCenter = layerGeometry->middleWireOfGroup( LCT_iter->getKeyWG() );
	      Int_t halfstrip_id = LCT_iter->getStrip();
#define GETLCTPOS LocalPoint LCTPos = layerGeometry->intersectionOfStripAndWire( halfstrip_id/2. , WGCenter ); \
	      Float_t new_dR2=pow(LCTPos.x()-(*chamberMatch)->x,2)+pow(LCTPos.y()-(*chamberMatch)->y,2); \
	      if ( new_dR2<dR2 ) {dR2=new_dR2;XthisLCT=LCTPos.x();YthisLCT=LCTPos.y();}
	      if( ( abs(StationRing->back())==14 || abs(StationRing->back())==11 ) && (*chamberMatch)->x!=999999 ){
		halfstrip_id=halfstrip_id%32;
		for(Byte_t ME14ReadOut = 0; ME14ReadOut < 3; ME14ReadOut++){
		  GETLCTPOS
		  //printf("ME%d:HS%3.1f LCT(%3.2f,%3.2f)--(%3.2fcm)--Track(%3.2f,%3.2f)\n",thisRing,halfstrip_id/2.,LCTPos.x(),LCTPos.y(),sqrt(new_dR2),(*chamberMatch)->x,(*chamberMatch)->y);
		  halfstrip_id += 32;
		} 
	      }
	      else {GETLCTPOS}
	    }//end of LCT_iter
	  }//end of CSCDigi_iter
	}//end of HasDigi
      }
      else {
	const DTChamberId ChamberID( (*chamberMatch)->id.rawId() );
	StationRing->push_back( ChamberID.wheel()>0?ChamberID.wheel()*10+ChamberID.station():ChamberID.wheel()*10-ChamberID.station() );
	SectorChamber->push_back( ChamberID.sector() );
	dYToDeadZone->push_back(999999.);
      }
      MuonIndex->push_back( eta->size()-1 );
      TrackDistToChamberEdge->push_back( (*chamberMatch)->dist() );
      TrackDistToChamberEdgeErr->push_back( (*chamberMatch)->distErr() );
      
      XTrack->push_back( (*chamberMatch)->x );
      YTrack->push_back( (*chamberMatch)->y );
      XErrTrack->push_back( (*chamberMatch)->xErr );
      YErrTrack->push_back( (*chamberMatch)->yErr );
      //LCT and segment information
      NumberOfLCTsInChamber->push_back( NumberOfLCTs );
      XLCT->push_back(XthisLCT);
      YLCT->push_back(YthisLCT);
      XSegment->push_back(999999);
      YSegment->push_back(999999);
      XErrSegment->push_back(999999);
      YErrSegment->push_back(999999);
      NumberOfHitsInSegment->push_back(0);
      IsSegmentOwnedExclusively->push_back( false );
      IsSegmentBelongsToTrackByDR->push_back( false );
      IsSegmentBelongsToTrackByCleaning->push_back( false );
      if ( segmentMatch != (*chamberMatch)->segmentMatches.end() ) {
	if (IsCSC) 
	  if (segmentMatch->cscSegmentRef.isNonnull()) NumberOfHitsInSegment->back() = segmentMatch->cscSegmentRef->nRecHits();
	  else LogWarning("RefNull")<<"ME "<<StationRing->back()<<", Segment Ref not found";
	else
	  if (segmentMatch->dtSegmentRef.isNonnull()) NumberOfHitsInSegment->back() = segmentMatch->dtSegmentRef->recHits().size();
	  else LogWarning("RefNull")<<"MB "<<StationRing->back()<<", Segment Ref not found";	   
	XSegment->back() = segmentMatch->x;
	YSegment->back() = segmentMatch->y;
	XErrSegment->back() = segmentMatch->xErr;
	YErrSegment->back() = segmentMatch->yErr;
	UInt_t CurrentChamber=XTrack->size()-1;
	IsSegmentOwnedExclusively->back()=true;
	for (UInt_t chamber=0; chamber<CurrentChamber&&MuonIndex->at(chamber)!=MuonIndex->back(); chamber++)
	  if ( IsTheSameSegment(chamber,CurrentChamber) ) {
	    IsSegmentOwnedExclusively->back()=false;
	    IsSegmentOwnedExclusively->at(chamber)=false;
	  }
	IsSegmentBelongsToTrackByDR->back() = segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR);
	IsSegmentBelongsToTrackByCleaning->back() = segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning);
      }
    }//chamber loop end
//################# Tracker Track and Segments MC Truth Match #######################
#ifdef TrackingParticles
    if (IsMC&&HasTrackingParticle) {
      SavedTP.clear();
      IsPUTrack->push_back(false);
      TrkParticles_pt->push_back(0);
      TrkParticles_eta->push_back(0);
      TrkParticles_phi->push_back(0);
      TrkParticles_pdgId->push_back(0);
      TrkParticles_charge->push_back(-100.0);
      SharedHitsRatio->push_back(-100.0);
      MCMatchChi2->push_back(-100.0);
      TTTruthDChains->push_back(-2);
      TTTruthMuType->push_back(0);
      SegTruthDChains->push_back(-2);
      SegTruthMuType->push_back(0);
      NumMisMatch++;
      ChamberSimHits.clear();
      RefToBase<Track> trk(InnerTrack);
      if (InnerTrack_nValidTrackerHits->back()>=minTrackHits&&RecoToSimByHits.find(trk) != RecoToSimByHits.end()&&RecoToSimByChi2.find(trk) != RecoToSimByChi2.end()) {
	pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trk].front();
	vector<pair<TrackingParticleRef, double> > TPCByChi2=RecoToSimByChi2[trk];
	vector<pair<TrackingParticleRef, double> >::iterator TPCByChi2_iter=TPCByChi2.begin();
	for (;TPCByChi2_iter!=TPCByChi2.end();TPCByChi2_iter++)
	  if (BestMatch.first==TPCByChi2_iter->first) break;
	if (TPCByChi2_iter!=TPCByChi2.end()) {
	  TrackingParticleRef tpr = BestMatch.first;
	  //Simulated Tracks
	  TrkParticles_pt->back()=tpr->pt();
	  TrkParticles_eta->back()=tpr->eta();
	  TrkParticles_phi->back()=tpr->phi();
	  TrkParticles_charge->back()=tpr->charge();
	  TrkParticles_pdgId->back()=tpr->pdgId();
	  SharedHitsRatio->back()=BestMatch.second;
	  MCMatchChi2->back()=-1.*TPCByChi2_iter->second;
	  NumMisMatch--;
	  //Get the decay chain of this track
	  vector <TheMuonType> type;
	  IsPUTrack->back()=true;
	  GetDecayChains(tpr,TTTruthDChains,type,HepGenEvent);
	  Long64_t TypeRecord=0;
	  for (vector <TheMuonType>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	    if (*type_iter) TypeRecord=TypeRecord*100+Long64_t(*type_iter);
	  TTTruthMuType->back()=TypeRecord;
	}
      }
      //CSC Segment Truth
      vector <TheMuonType> type;SavedTP.clear();
      UInt_t NumOfSegments=XSegment->size();
      for (UInt_t Seg=0; Seg<NumOfSegments; Seg++)
	if ( XSegment->at(Seg)<99999. && MuonIndex->at(Seg)==eta->size()-1 && IsSegmentOwnedExclusively->at(Seg) && IsCSCChamber->at(Seg) ) {
	  Float_t DR2=999999.; Int_t tpIndex=-2;
	  UInt_t TPCollectionSize = TPCollectionH.product()->size();
	  for( UInt_t tp_iter = 0; tp_iter < TPCollectionSize; tp_iter++ ) {
	    TrackingParticleRef tpr(TPCollectionH,tp_iter);
	    if ( find(SavedTP.begin(),SavedTP.end(),tpr)!=SavedTP.end() ) continue;
	    for ( vector<PSimHit>::const_iterator g4Hit_iter=tpr->pSimHit_begin();g4Hit_iter!=tpr->pSimHit_end();g4Hit_iter++ ) {
	      DetId SimHitDetectorId( g4Hit_iter->detUnitId() );
	      if ( SimHitDetectorId.det() != DetId::Muon ) continue;
	      if ( SimHitDetectorId.subdetId() == MuonSubdetId::CSC ) {
		const CSCDetId SimHitChamberID( g4Hit_iter->detUnitId() );
		Int_t SimHitRing = SimHitChamberID.station()*10+SimHitChamberID.ring();
		if ( SimHitChamberID.endcap()==2 ) SimHitRing=-SimHitRing;
		if ( StationRing->at(Seg) != SimHitRing || SectorChamber->at(Seg) != SimHitChamberID.chamber() ) continue;
		Local3DPoint entryPoint = g4Hit_iter->entryPoint(), exitPoint = g4Hit_iter->exitPoint();
		LocalVector direction = exitPoint - entryPoint;
		if ( fabs(direction.z()) <= 0.001) continue;
		LocalPoint projection = entryPoint - direction*(entryPoint.z()/direction.z());
		if ( fabs(projection.z()) > 0.001 ) continue;
		Float_t newDR2=pow(projection.x()-XSegment->at(Seg),2)+pow(projection.y()-YSegment->at(Seg),2);
		if (newDR2<DR2) {
		  DR2=newDR2;
		  tpIndex=Int_t(tp_iter);
		}
	      }
	    }
	  }
	  if ( tpIndex>=0 ) GetDecayChains(TrackingParticleRef(TPCollectionH,tpIndex),SegTruthDChains,type,HepGenEvent);
	}
      ULong64_t TypeRecord=0;
      for (vector <TheMuonType>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	if (*type_iter) TypeRecord=TypeRecord*100+Long64_t(*type_iter);
      SegTruthMuType->back()=TypeRecord;
    }//end of IsMC
#endif
  }//muon loop end
#ifdef DiMuonAnalysis
  TLorentzVector P1,P2,Pc;
  unsigned int num_muons=pt->size();
  for (unsigned int Mu1=0; Mu1<num_muons; Mu1++)
    for (unsigned int Mu2=Mu1+1; Mu2<num_muons; Mu2++) {
      P1.SetPtEtaPhiM((*pt)[Mu1],(*eta)[Mu1],(*phi)[Mu1],Muon_Mass);
      P2.SetPtEtaPhiM((*pt)[Mu2],(*eta)[Mu2],(*phi)[Mu2],Muon_Mass);
      Float_t  InvarMass=(P1+P2).M();
      DiMuonInvariantMass->push_back(InvarMass);
      Pc=P1+P2;
      Float_t Gamma=Pc.Gamma();
      TLorentzRotation l;
      l.Boost(Pc.Px()/Gamma/InvarMass,Pc.Py()/Gamma/InvarMass,Pc.Pz()/Gamma/InvarMass);
      l.Invert();
      P1.Transform(l);
      CosThetaStar->push_back(P1.CosTheta());
      AngleBetweenDiMuon->push_back( acos( -(P1.Px()*P2.Px()+P1.Py()*P2.Py()+P1.Pz()*P2.Pz())/P1.P()/P2.P() ) ); 
    }
#endif
  Muons_Tree->Fill();
  delete HepGenEvent;
  return;
}

void MuonSysPriEff::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup) {
  Bool_t changed=true;
  if (NumberOfEvents!=0) RunInfo_Tree->Fill();
  NumberOfEvents=0;
  if (HLTConfig.init(iRun,iSetup,HLTProc,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
      HLTTableName->assign(HLTConfig.tableName());
      *HLTNamesSet=HLTConfig.triggerNames();
      LogInfo("DataFormat")<<"HLT table is changed to";
      HLTConfig.dump("TableName");
      for (UInt_t j=0;j<num_HLTsSaveObjs;j++) {
	UInt_t i,length=HLTObj_HLTNames[j].rfind("_v");
	if (length<6) length=HLTObj_HLTNames[j].size();
	HLTSize=HLTNamesSet->size();
	for (i=0;i<HLTSize;i++)
	  if ( HLTObj_HLTNames[j].compare(0,length,HLTNamesSet->at(i),0,length)==0 ) {HLTFilterNamesAcceptenceIndex[j]=i;break;}//ignore the versions of the HLT name
	if (i==HLTSize) {
	  LogWarning("DataLost")<<"the HLT object request:"<<HLTObj_HLTNames[j].substr(0,length)<<" does not exist in the current HLT scope";
	  HLTFilterNamesAcceptenceIndex[j]=-1;
	}
	else //for (vector<string>::const_iterator iter=HLTConfig.moduleLabels(HLTFilterNamesAcceptenceIndex[j]).begin();iter!=HLTConfig.moduleLabels(HLTFilterNamesAcceptenceIndex[j]).end();iter++)
	  HLT_ModuleNames[j]=InputTag((HLTConfig.moduleLabels(HLTFilterNamesAcceptenceIndex[j]).end()-2)->c_str(),"",HLTProc.c_str());
      }
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    LogWarning("DataLost") << " HLT config extraction failure with process name " <<HLTProc;
    // In this case, all access methods will return empty values!
  }
}

//deadzone center is according to http://cmslxr.fnal.gov/lxr/source/RecoLocalMuon/CSCEfficiency/src/CSCEfficiency.cc#605
//wire spacing is according to CSCTDR
const Float_t deadZoneCenterME1_2[2] = {-32.88305,32.867423};
const Float_t deadZoneCenterME1_3[2] = {-22.7401,27.86665};
const Float_t deadZoneCenterME2_1[2] = {-27.47,33.67};
const Float_t deadZoneCenterME3_1[2] = {-36.21,23.68};
const Float_t deadZoneCenterME4_1[2] = {-26.14,23.85};
const Float_t deadZoneCenterME234_2[4] = {-81.8744,-21.18165,39.51105,100.2939};
Float_t MuonSysPriEff::YDistToHVDeadZone(Float_t yLocal, Int_t StationAndRing){
  const Float_t *deadZoneCenter;
  Float_t deadZoneHeightHalf=3.16*7./2.,minY=999999.;
  Byte_t nGaps=2;
  switch (abs(StationAndRing)) {
  case 11:
  case 14:
    return minY;
  case 12:
    deadZoneCenter=deadZoneCenterME1_2;
    break;
  case 13:
    deadZoneCenter=deadZoneCenterME1_3;
    break;
  case 21:
    deadZoneCenter=deadZoneCenterME2_1;
    deadZoneHeightHalf=3.12*7./2.;
    break;
  case 31:
    deadZoneCenter=deadZoneCenterME3_1;
    deadZoneHeightHalf=3.12*7./2.;
    break;
  case 41:
    deadZoneCenter=deadZoneCenterME4_1;
    deadZoneHeightHalf=3.12*7./2.;
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

#define MakeVecBranch(Name,Var,Type) Var=new vector<Type>();Muons_Tree->Branch(Name,&Var)
MuonSysPriEff::MuonSysPriEff(const edm::ParameterSet& pset) {

//###############Get the input parameters###################
//Standard filters
  HLTFilter=new TriggerResultsFilter(pset.getParameter<edm::ParameterSet>("HLTFilterPSet"));
  VertexFilter=new GoodVertexFilter(pset.getParameter<edm::ParameterSet>("VertexFilterPSet"));
  ScrapingFilter=new FilterOutScraping(pset.getParameter<edm::ParameterSet>("noscrapingPSet"));
//Output filename
  FileName = pset.getParameter<string>("FileName");
//Muon Cut
  MuonPtCut = pset.getParameter<double>("MuonPtCut");
  MinNumberOfMuons = pset.getUntrackedParameter<unsigned int>("MinNumberOfMuons",1);
  MuonEtaMin = pset.getParameter<double>("MuonEtaMin");
  MuonEtaMax = pset.getParameter<double>("MuonEtaMax");
  ThrowBadEvents = pset.getParameter<bool>("ThrowBadEvents");
 //  HasTrackingParticle = pset.getUntrackedParameter<bool>("HasTrackingParticle",false);
//standard muon types e.g. TMLastStation......
  vector<string> Cuts;
  Cuts = pset.getUntrackedParameter< vector<string> >("StandardMuonCuts",Cuts);//default is empty
  string MuonArbitrationTypeStr;
  MuonArbitrationTypeStr = pset.getUntrackedParameter<string>("MuonArbitrationType","SegmentArbitration");
  MuonArbitrationType=MuonArbitrationTypeFromString(MuonArbitrationTypeStr.c_str());
//PV
  PrimaryVerticesTag = pset.getUntrackedParameter<string>("PrimaryVertices","offlinePrimaryVertices");//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
//tracks
  tracksTag = pset.getUntrackedParameter<InputTag>("tracksTag",InputTag("generalTracks"));
  TPInputTag = pset.getUntrackedParameter<InputTag>("TPInputTag",InputTag("mergedtruth","MergedTrackTruth") );
  CSCDigisTag = pset.getParameter<InputTag>("CSCDigisTag");
  dEdxTag = pset.getUntrackedParameter<string>("dEdxTag","dedxHarmonic2");//Other options already available in RECO files are dedxMedian and dedxTruncated40. 
//Triggers
  TriggerResultsTag = pset.getParameter<InputTag>("TriggerResultsTag");
  HLTProc=TriggerResultsTag.process();
  TriggerEventTag = pset.getUntrackedParameter<InputTag>("triggerEventTag",InputTag("hltTriggerSummaryAOD","",HLTProc.c_str()));
  HLTObj_HLTNames = pset.getParameter< vector<string> >("HLTObj_HLTNames");
  num_HLTsSaveObjs = HLTObj_HLTNames.size();
  HLT_ModuleNames.assign(num_HLTsSaveObjs,InputTag(""));
  if (num_HLTsSaveObjs>maxFilterObjects) throw cms::Exception("You need to increase the constant maxFilterObjects.\n");;
//Generation Information
  HepMCTag = pset.getUntrackedParameter<string>("HepMCTag","generator");

//#################### Make a root file and plant TTrees #########################
  file=new TFile(FileName.c_str(),"RECREATE");
  Muons_Tree = new TTree ("MuTrkCand","TrksCandidates") ;
  RunInfo_Tree = new TTree ("RunInfo","RunInfo") ;
  Muons_Tree->SetCircular(500000);
  RunInfo_Tree->SetCircular(500000);//the max events in a single root file
//#################### Build Branches #########################
  //General event information
  Muons_Tree->Branch("Event_Info",&Info.RUN,"RUN/l:EVENT:LumiBlock:ORBIT:BunchCrossing");
  Muons_Tree->Branch("isRealData",&isRealData,"isRealData/O");
  Muons_Tree->Branch("isHLTTriggerred",&isHLTTriggerred,"isHLTTriggerred/O");
  Muons_Tree->Branch("isGoodVertex",&isGoodVertex,"isGoodVertex/O");
  Muons_Tree->Branch("isNoScrapping",&isNoScrapping,"isNoScrapping/O");
  Muons_Tree->Branch("GenEventWeight",&GenEventWeight,"GenEventWeight/D");
  numberOfPUVertices=0;
  Muons_Tree->Branch("numberOfPUVertices",&numberOfPUVertices,"numberOfPUVertices/I");
  //Muon Basic Information
  MakeVecBranch("pt",pt,Float_t); MakeVecBranch("ptError",ptError,Float_t); 
  MakeVecBranch("eta",eta,Float_t); MakeVecBranch("etaError",etaError,Float_t);  
  MakeVecBranch("phi",phi,Float_t); MakeVecBranch("phiError",phiError,Float_t);
  MakeVecBranch("chargeMinus",chargeMinus,Bool_t); MakeVecBranch("isGlobalMu",isGlobalMu,Bool_t);  MakeVecBranch("isTrackerMu",isTrackerMu,Bool_t);
  //Muon Vertex
  MakeVecBranch("Vertex_x",Vertex_x,Float_t);  MakeVecBranch("Vertex_y",Vertex_y,Float_t);  MakeVecBranch("Vertex_z",Vertex_z,Float_t);
  MakeVecBranch("DXYwtBS",DXYwtBS,Float_t);
  MakeVecBranch("DZwtBS",DZwtBS,Float_t);
  MakeVecBranch("DXYwtPV",DXYwtPV,Float_t); MakeVecBranch("DXYErrwtPV",DXYErrwtPV,Float_t);
  //Muon Isolation
  MakeVecBranch("isoR03sumPt",isoR03sumPt,Float_t);
  MakeVecBranch("isoR05sumPt",isoR05sumPt,Float_t);
  MakeVecBranch("TrkKink",TrkKink,Float_t); MakeVecBranch("GlbKink",GlbKink,Float_t); MakeVecBranch("TrkRelChi2",TrkRelChi2,Float_t);
  //trackquality
  MakeVecBranch("numberOfMatchedSegments",numberOfMatchedSegments,Int_t);
  MakeVecBranch("numberOfMatchedStations",numberOfMatchedStations,Int_t);
  MakeVecBranch("dEdx",dEdx,Float_t);  MakeVecBranch("dEdxError",dEdxError,Float_t);
  MakeVecBranch("dEdx_numberOfSaturatedMeasurements",dEdx_numberOfSaturatedMeasurements,Int_t);
  MakeVecBranch("dEdx_numberOfMeasurements",dEdx_numberOfMeasurements,Int_t);
  MakeVecBranch("Track_nValidTrackerHits",InnerTrack_nValidTrackerHits,UInt_t);
  MakeVecBranch("Track_nValidPixelHits",InnerTrack_nValidPixelHits,UInt_t);
  MakeVecBranch("Track_nLostTrackerHits",InnerTrack_nLostTrackerHits,UInt_t);
  MakeVecBranch("Track_nLostPixelHits",InnerTrack_nLostPixelHits,UInt_t);
  MakeVecBranch("numberOfValidMuonHits",numberOfValidMuonHits,UInt_t);
  MakeVecBranch("InnerTrack_chi2",InnerTrack_chi2,Float_t);  MakeVecBranch("InnerTrack_ndof",InnerTrack_ndof,UInt_t);
  MakeVecBranch("GlobalTrack_chi2",GlobalTrack_chi2,Float_t);  MakeVecBranch("GlobalTrack_ndof",GlobalTrack_ndof,UInt_t);
#ifdef DiMuonAnalysis
  DiMuonInvarMassCut = pset.getUntrackedParameter<double>("DiMuonInvarMassCut",0.);
  //Information of two muons
  MakeVecBranch("DiMuonInvariantMass",DiMuonInvariantMass,Float_t);
  MakeVecBranch("CosThetaStar",CosThetaStar,Float_t);
  MakeVecBranch("AngleBetweenDiMuon",AngleBetweenDiMuon,Float_t);
#endif
  //Muon System - Segment/Chambers
  MakeVecBranch("MuonIndex",MuonIndex,Byte_t);
  MakeVecBranch("IsCSCChamber",IsCSCChamber,Bool_t);
  MakeVecBranch("SectorChamber",SectorChamber,Byte_t);
  MakeVecBranch("StationRing",StationRing,Char_t);
  MakeVecBranch("TrackDistToChamberEdge",TrackDistToChamberEdge,Float_t);
  MakeVecBranch("TrackDistToChamberEdgeErr",TrackDistToChamberEdgeErr,Float_t);
  MakeVecBranch("XTrack",XTrack,Float_t);
  MakeVecBranch("YTrack",YTrack,Float_t);
  MakeVecBranch("XErrTrack",XErrTrack,Float_t);
  MakeVecBranch("YErrTrack",YErrTrack,Float_t);
  MakeVecBranch("dYToDeadZone",dYToDeadZone,Float_t);
  //CSC LCT Information
  MakeVecBranch("NumberOfLCTsInChamber",NumberOfLCTsInChamber,Byte_t);
  MakeVecBranch("XLCT",XLCT,Float_t);
  MakeVecBranch("YLCT",YLCT,Float_t);
  //Muon Segment Information
  MakeVecBranch("XSegment",XSegment,Float_t);
  MakeVecBranch("YSegment",YSegment,Float_t);
  MakeVecBranch("XErrSegment",XErrSegment,Float_t);
  MakeVecBranch("YErrSegment",YErrSegment,Float_t);
  MakeVecBranch("IsSegmentOwnedExclusively",IsSegmentOwnedExclusively,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByDR",IsSegmentBelongsToTrackByDR,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByCleaning",IsSegmentBelongsToTrackByCleaning,Bool_t);
  MakeVecBranch("NumberOfHitsInSegment",NumberOfHitsInSegment,Byte_t);
  //HLT
  MakeVecBranch("HLTacceptance",HLTacceptance,Bool_t);
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    MakeVecBranch((HLTObj_HLTNames[i]+"_pt").c_str(),HLTObj_pt[i],Float_t);
    MakeVecBranch((HLTObj_HLTNames[i]+"_eta").c_str(),HLTObj_eta[i],Float_t);
    MakeVecBranch((HLTObj_HLTNames[i]+"_phi").c_str(),HLTObj_phi[i],Float_t);
    MakeVecBranch(("is"+HLTObj_HLTNames[i]+"Obj").c_str(),isHLTObj[i],Bool_t);
  }
  //PV
  MakeVecBranch("PVx",vx,Float_t);  MakeVecBranch("PVy",vy,Float_t);  MakeVecBranch("PVz",vz,Float_t);
  MakeVecBranch("PVxError",vxError,Float_t);  MakeVecBranch("PVyError",vyError,Float_t);  MakeVecBranch("PVzError",vzError,Float_t);
  //Cuts
  num_Cuts=Cuts.size();
  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++) 
    if (muon::selectionTypeFromString(Cuts[whichcut])!=(muon::SelectionType)-1) {
      MakeVecBranch(Cuts[whichcut].c_str(),SelectorPassed[whichcut],Bool_t); 
      OfficialMuonSelectors.push_back(muon::selectionTypeFromString(Cuts[whichcut]));
    }
  num_Cuts=OfficialMuonSelectors.size();
  
  //Generation and Simulation Particles
  MakeVecBranch("Gen_pt",Gen_pt,Float_t);  MakeVecBranch("Gen_eta",Gen_eta,Float_t);  MakeVecBranch("Gen_phi",Gen_phi,Float_t);   MakeVecBranch("Gen_pdgId",Gen_pdgId,Int_t);
  MakeVecBranch("Gen_vx",Gen_vx,Float_t);  MakeVecBranch("Gen_vy",Gen_vy,Float_t);  MakeVecBranch("Gen_vz",Gen_vz,Float_t);  MakeVecBranch("Gen_vt",Gen_vt,Float_t);
  MakeVecBranch("IsParInHep",IsParInHep,Bool_t);

#ifdef TrackingParticles  
  //Simulated Tracks
  minTrackHits = pset.getUntrackedParameter<uint>("minTrackHits",3);
  MakeVecBranch("IsParHasMuonHits",IsParHasMuonHits,Bool_t);
  MakeVecBranch("IsPUTrack",IsPUTrack,Bool_t);
  MakeVecBranch("TrkParticles_pt",TrkParticles_pt,Float_t);  MakeVecBranch("TrkParticles_eta",TrkParticles_eta,Float_t);  MakeVecBranch("TrkParticles_phi",TrkParticles_phi,Float_t);
  MakeVecBranch("TrkParticles_pdgId",TrkParticles_pdgId,Int_t); MakeVecBranch("TrkParticles_charge",TrkParticles_charge,Int_t);
  MakeVecBranch("SharedHitsRatio",SharedHitsRatio,Double_t); MakeVecBranch("MCMatchChi2",MCMatchChi2,Double_t); 
  MakeVecBranch("theSameWithMuon",theSameWithMuon,Int_t);
  MakeVecBranch("TTTruthDChains",TTTruthDChains,Int_t);
  MakeVecBranch("TTTruthMuType",TTTruthMuType,Long64_t);
  MakeVecBranch("SegTruthDChains",SegTruthDChains,Int_t);
  MakeVecBranch("SegTruthMuType",SegTruthMuType,Long64_t);
  NumMisMatch=0;
#endif

  //Summary
  HLTTableName = new string();
  RunInfo_Tree->Branch("RunNumber",&Info.RUN,"RunNum/l");
  RunInfo_Tree->Branch("TableName",&HLTTableName);
  HLTNamesSet = new vector<string>();
  RunInfo_Tree->Branch("HLTNamesSet",&HLTNamesSet);
  RunInfo_Tree->Branch("NumberOfEvents",&NumberOfEvents,"NumberOfEvents/l");
  FirstEvent=true;IsMC=true;HasDigi=true;HasTrackingParticle=true;
}

Muon::ArbitrationType MuonSysPriEff::MuonArbitrationTypeFromString( const std::string &label ) {
  struct MuonArbitrationTypeStringToEnum { const char *label; reco::Muon::ArbitrationType value; };
  static MuonArbitrationTypeStringToEnum MuonArbitrationTypeStringToEnumMap[] = {
    {"NoArbitration",Muon::NoArbitration},
    {"SegmentArbitration",Muon::SegmentArbitration},
    {"SegmentAndTrackArbitration",Muon::SegmentAndTrackArbitration},
    {"SegmentAndTrackArbitrationCleaned",Muon::SegmentAndTrackArbitrationCleaned},
    { 0, (Muon::ArbitrationType)-1 }
  };
  Muon::ArbitrationType value = (Muon::ArbitrationType)-1;
  bool found = false;
  for(int i = 0; MuonArbitrationTypeStringToEnumMap[i].label && (! found); ++i)
    if (! strcmp(label.c_str(), MuonArbitrationTypeStringToEnumMap[i].label)) {
      found = true;
      value = MuonArbitrationTypeStringToEnumMap[i].value;
    }
  
  // in case of unrecognized selection type
  if (! found) throw cms::Exception("MuonSelectorError") << label << " is not a recognized SelectionType";
  return value;
}

void MuonSysPriEff::ClearVecs_HepMC() {
  Gen_pt->clear();  Gen_eta->clear();  Gen_phi->clear(); Gen_pdgId->clear();
  Gen_vx->clear();  Gen_vy->clear();  Gen_vz->clear();  Gen_vt->clear();
  SavedHepPar.clear();
}

void MuonSysPriEff::ClearVecs_RECO() {
  //basic information
  pt->clear(); eta->clear();  phi->clear();
  ptError->clear(); etaError->clear();  phiError->clear();
  chargeMinus->clear();  isGlobalMu->clear();  isTrackerMu->clear();
  //muon vertex
  Vertex_x->clear();  Vertex_y->clear();  Vertex_z->clear();
  DXYwtBS->clear();
  DZwtBS->clear();
  DXYwtPV->clear(); DXYErrwtPV->clear();
  //isolation
  isoR03sumPt->clear();
  isoR05sumPt->clear();
  //track quality
  TrkKink->clear();  GlbKink->clear();  TrkRelChi2->clear();
  dEdx->clear();  dEdxError->clear();
  dEdx_numberOfSaturatedMeasurements->clear();
  dEdx_numberOfMeasurements->clear();
  InnerTrack_nValidTrackerHits->clear();
  InnerTrack_nValidPixelHits->clear();  
  InnerTrack_nLostTrackerHits->clear();
  InnerTrack_nLostPixelHits->clear(); 
  InnerTrack_chi2->clear(); InnerTrack_ndof->clear();
  numberOfValidMuonHits->clear();
  numberOfMatchedSegments->clear();
  numberOfMatchedStations->clear();
  GlobalTrack_chi2->clear(); GlobalTrack_ndof->clear();
#ifdef DiMuonAnalysis
  //Information of two muons
  DiMuonInvariantMass->clear();  CosThetaStar->clear(); AngleBetweenDiMuon->clear();
#endif
  //HLT
  HLTacceptance->clear();
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    HLTObj_pt[i]->clear();
    HLTObj_eta[i]->clear();
    HLTObj_phi[i]->clear();
    isHLTObj[i]->clear();
  }

  //Chamber
  TrackDistToChamberEdge->clear();  TrackDistToChamberEdgeErr->clear();
  XTrack->clear();YTrack->clear();  XErrTrack->clear();YErrTrack->clear();
  IsCSCChamber->clear(); StationRing->clear();SectorChamber->clear();MuonIndex->clear();
  dYToDeadZone->clear();
  //LCT
  NumberOfLCTsInChamber->clear();
  XLCT->clear();  YLCT->clear();

  //Segment
  XSegment->clear();  YSegment->clear();
  XErrSegment->clear();  YErrSegment->clear();
  IsSegmentOwnedExclusively->clear();
  IsSegmentBelongsToTrackByDR->clear();  IsSegmentBelongsToTrackByCleaning->clear();
  NumberOfHitsInSegment->clear();

  vx->clear(); vy->clear(); vz->clear();
  vxError->clear(); vyError->clear(); vzError->clear();

  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    SelectorPassed[whichcut]->clear(); 

#ifdef TrackingParticles
  TrkParticles_pt->clear();  TrkParticles_eta->clear(); TrkParticles_phi->clear();
  TrkParticles_pdgId->clear(); TrkParticles_charge->clear(); 
  SharedHitsRatio->clear(); MCMatchChi2->clear(); SavedSimTrk.clear(); 
  theSameWithMuon->clear();
  IsParInHep->clear(); IsParHasMuonHits->clear();
  IsPUTrack->clear();
  TTTruthDChains->clear(); TTTruthMuType->clear();
  SegTruthDChains->clear(); SegTruthMuType->clear();
#endif
}

// ------- method called once each job just before starting event loop  ------------
void MuonSysPriEff::beginJob(){}

// -------- method called once each job just after ending the event loop  ------------
void MuonSysPriEff::endJob() {
  RunInfo_Tree->Fill();
  file->Write();
  file->Close();
}

MuonSysPriEff::~MuonSysPriEff() {
  //standard filters
  HLTFilter->~TriggerResultsFilter();
  VertexFilter->~GoodVertexFilter();
  ScrapingFilter->~FilterOutScraping();
  //Generation Information
  delete Gen_pt; delete Gen_eta; delete Gen_phi; delete Gen_pdgId;
  delete Gen_vx; delete Gen_vy;  delete Gen_vz;  delete Gen_vt;
 //basic information
  delete pt; delete eta; delete  phi;
  delete ptError; delete etaError; delete phiError;
  delete chargeMinus;
  delete  isGlobalMu; delete  isTrackerMu;
  //muon vertex
  delete Vertex_x; delete  Vertex_y; delete  Vertex_z;
  delete DXYwtBS;
  delete DZwtBS;
  delete  DXYwtPV; delete DXYErrwtPV;
  //isolation
  delete isoR03sumPt;
  delete isoR05sumPt;
  //track quality
  delete TrkKink; delete GlbKink; delete TrkRelChi2;
  delete dEdx; delete  dEdxError;
  delete dEdx_numberOfSaturatedMeasurements;
  delete dEdx_numberOfMeasurements;
  delete InnerTrack_nValidTrackerHits;
  delete InnerTrack_nValidPixelHits;
  delete InnerTrack_nLostTrackerHits;
  delete InnerTrack_nLostPixelHits;
  delete InnerTrack_chi2; delete InnerTrack_ndof;
  delete numberOfValidMuonHits;
  delete numberOfMatchedSegments;
  delete numberOfMatchedStations;
  delete GlobalTrack_chi2; delete GlobalTrack_ndof;
#ifdef DiMuonAnalysis
  //Information of two muons
  delete DiMuonInvariantMass;
  delete  CosThetaStar; delete AngleBetweenDiMuon;
#endif
  //HLT
  delete HLTTableName;
  delete HLTacceptance;
  for( unsigned int i =0; i < num_HLTsSaveObjs;i++) {
    delete HLTObj_pt[i];
    delete HLTObj_eta[i];
    delete HLTObj_phi[i];
    delete isHLTObj[i];
  }

  //Chamber
  delete TrackDistToChamberEdge; delete TrackDistToChamberEdgeErr;
  delete XTrack; delete YTrack; delete XErrTrack; delete YErrTrack;
  delete IsCSCChamber; delete StationRing; delete SectorChamber;
  delete MuonIndex;
  delete dYToDeadZone;
  //LCT
  delete NumberOfLCTsInChamber;
  delete XLCT;  delete YLCT;
  //Segment
  delete XSegment; delete  YSegment;
  delete XErrSegment; delete  YErrSegment;
  delete IsSegmentOwnedExclusively; delete IsSegmentBelongsToTrackByDR;
  delete IsSegmentBelongsToTrackByCleaning; delete NumberOfHitsInSegment;
 
  delete vx; delete vy; delete vz;
  delete vxError; delete vyError; delete vzError;

  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    delete SelectorPassed[whichcut];
#ifdef TrackingParticles
  delete TrkParticles_pt; delete TrkParticles_eta; delete TrkParticles_phi;
  delete TrkParticles_pdgId; delete  TrkParticles_charge;
  delete SharedHitsRatio; delete MCMatchChi2;
  delete theSameWithMuon;
  delete IsParInHep; delete IsParHasMuonHits;
  delete IsPUTrack;
  delete TTTruthDChains;
  delete TTTruthMuType;
  delete SegTruthDChains;
  delete SegTruthMuType;
#endif
}

#ifdef TrackingParticles
//particle classifiction based on pdg_id
MuonSysPriEff::ParticleType MuonSysPriEff::ParticleCata(int pid)
{
  pid=pid>0?pid:pid*-1;
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
//Classification
MuonSysPriEff::TheMuonType MuonSysPriEff::classification(vector<Int_t> &Chain)
{
  int MuPos=-1;
  for (vector<Int_t>::reverse_iterator iter = Chain.rbegin(); iter != Chain.rend(); iter++) {
    int pid=(*Gen_pdgId)[*iter];
    ParticleType ParType=ParticleCata(pid);
    if (abs(pid)==13) MuPos=*iter;
    if (MuPos>=0&&(ParType==W||ParType==Z)) {
      if ((*IsParInHep)[MuPos]) return PromptMuFromWZ;
      else return NotPromptMufromWZ;
    }
    if (MuPos<0&&ParType>=LightMeson&&ParType<=DiQuarks) {
      if ( (*IsParHasMuonHits)[*iter] ) return PunchThrough;
      else return NoMuSysHit;
    }
    if (MuPos>=0&&ParType>=LightMeson&&ParType<=BottomBaryon) {
      float deltaRPhi2=((*Gen_vx)[MuPos])*((*Gen_vx)[MuPos])+((*Gen_vy)[MuPos])*((*Gen_vy)[MuPos]);
      float deltaZ=abs((*Gen_vz)[MuPos]);
      if ( (*IsParInHep)[MuPos]||(deltaRPhi2<100&&deltaZ<30) ) {//inside pixel detector
	if (ParType==LightMeson) return PromptMuFromLightMeson;
	if (ParType>LightMeson&&ParType<LightBaryon) return PromptMuFromHeavyMeson;
	if (ParType==LightBaryon) return PromptMuFromLightBaryon;
	if (ParType>LightBaryon&&ParType<=BottomBaryon) return PromptMuFromHeavyBaryon;
      }
      //if ( (RPhi2<161604.&&Z<568.)||(RPhi2>82024.96&&RPhi2<161604.&&Z<666.)) //inside HCAL
      if ( (*IsParHasMuonHits)[*iter] ) return PunchThroughAndDecayInFlight;
      else {
	if (ParType==LightMeson) return DecayInFlightFromLightMeson;
	if (ParType>LightMeson&&ParType<LightBaryon) return DecayInFlightFromHeavyMeson;
	if (ParType==LightBaryon) return DecayInFlightFromLightBaryon;
	if (ParType>LightBaryon&&ParType<=BottomBaryon) return DecayInFlightFromHeavyBaryon;
      }
    }
  }
  if (MuPos>=0) return PromptMuFromOthers;
  return Others;
}

//changed from MuonIdTruthInfo.cc: that one saves the layer but not chamber infomation, the return value is whether or not it is hit in muon system
Bool_t MuonSysPriEff::SimHitsToSegment(const PSimHit& hit) {
   // find the hit position projection at the reference surface of the layer:
   // first get entry and exit point of the hit in the global coordinates, then
   // get local coordinates of these points wrt the chamber and then find the
   // projected X-Y coordinates
  DetId DetectorId( hit.detUnitId() );
  if ( DetectorId.det() != DetId::Muon ) return false;
  else return true;
  /*Not using
  if ( DetectorId.subdetId() != MuonSubdetId::CSC ) return true;
  const CSCDetId CSCChamberID( hit.detUnitId() );
  Byte_t thisLayer=CSCChamberID.layer();
  Int_t ThisChamber=CSCChamberID.station()*10+CSCChamberID.ring();
  //reject non-existing CSC Chambers (incorrect geometry) In ME42, only ME+42/9-13 exist
  if ( ThisChamber==-42 ) return true;
  if ( ThisChamber==42 && ( CSCChamberID.chamber()<9 || CSCChamberID.chamber()>13 ) ) return true;
  if ( ThisChamber==-14 ) ThisChamber=-11;
  if ( ThisChamber==14 ) ThisChamber=11;//treat all ME14 as ME11 (they are the same ring)
  ThisChamber=ThisChamber*100+CSCChamberID.chamber();
  if ( CSCChamberID.endcap()==2 ) ThisChamber=-ThisChamber;

  Local3DPoint entryPoint = hit.entryPoint(), exitPoint = hit.exitPoint();
  LocalVector direction = exitPoint - entryPoint;
  if ( fabs(direction.z()) <= 0.001) return true;//won't save as sim seg. It is not a good sim hit, but it is a muon hit
  LocalPoint projection = entryPoint - direction*(entryPoint.z()/direction.z());
  if ( fabs(projection.z()) > 0.001 ) {
    LogWarning("IncorrectCalculation")<<"z coordinate of the hit projection must be zero and it's not!";
    return true;
  }
  
  vector<CSCChamberSimHitsInfo>::iterator thisChamberSimHits=ChamberSimHits.end(); Float_t DR2=64.;//max 8 cm to consider that they are in the same segment
  for ( vector<CSCChamberSimHitsInfo>::iterator ChamberSimHits_iter=ChamberSimHits.begin();ChamberSimHits_iter!=ChamberSimHits.end();ChamberSimHits_iter++)
    if ( ChamberSimHits_iter->ChamberID==ThisChamber ) 
      if ( ! (1<<thisLayer & ChamberSimHits_iter->HitsMask) ) {
	Float_t new_DR2=pow(projection.x()-(ChamberSimHits_iter->FirstHit_X+ChamberSimHits_iter->LastHit_X)/2.,2)+pow(projection.y()-(ChamberSimHits_iter->FirstHit_Y+ChamberSimHits_iter->LastHit_Y)/2.,2);
	if ( new_DR2<DR2 ) {
	  DR2=new_DR2;
	  thisChamberSimHits=ChamberSimHits_iter;
	}
      }

  if ( thisChamberSimHits == ChamberSimHits.end() ) {
    ChamberSimHits.push_back( CSCChamberSimHitsInfo( ThisChamber,projection.x(),projection.y(),projection.x(),projection.y(),1<<thisLayer ) );
    return true;
  }
  if ( ! (thisChamberSimHits->HitsMask>>thisLayer) ) {
    thisChamberSimHits->LastHit_X=projection.x();
    thisChamberSimHits->LastHit_Y=projection.y();
  }
  else if ( ! (thisChamberSimHits->HitsMask<<(8-thisLayer) ) ) {
    thisChamberSimHits->FirstHit_X=projection.x();
    thisChamberSimHits->FirstHit_Y=projection.y();
  }
  thisChamberSimHits->HitsMask |= 1<<thisLayer;
  return true;
  */
}

void MuonSysPriEff::GetDecayChains(TrackingParticleRef tpr,vector<int> *DChains, vector <TheMuonType> &type, HepMC::GenEvent *HepGenEvent) {
  for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )  {//g4Track loop begin
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
	      Bool_t HasMuonHits=false;
	      for ( vector<PSimHit>::const_iterator g4Hit_iter=tpr_tmp->pSimHit_begin();g4Hit_iter!=tpr_tmp->pSimHit_end();g4Hit_iter++ )
		if ( g4Hit_iter->trackId()==thisTrk->trackId() && g4Hit_iter->eventId().rawId()==thisTrk->eventId().rawId() && g4Hit_iter->particleType() == thisTrk->type() )  HasMuonHits=SimHitsToSegment( *g4Hit_iter );
	      vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
	      if (SavedSimTrk_iter==SavedSimTrk.end())
		{
		  DChain.push_back(IsParInHep->size());
		  IsParHasMuonHits->push_back(HasMuonHits);
		  RecordSimTrack(thisTrk)
		}
	      else DChain.push_back(FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin()));
	      SavedTP.push_back(tpr_tmp);
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
    if (!thisTrk->noGenpart()) {
      vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
      if (SavedSimTrk_iter!=SavedSimTrk.end()) (*IsParInHep)[FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin())]=true;
      else LogError("CodeWrong")<<"Cannot find the simulated track in saved sim tracks.";
      HepMC::GenEvent::particle_iterator genPar = HepGenEvent->particles_begin();
      for (int count=1; count<thisTrk->genpartIndex()&&genPar != HepGenEvent->particles_end(); count++ )
	genPar++;
      if (genPar != HepGenEvent->particles_end()) HepMCParentTree(*genPar);
      else LogWarning("RefNull")<<"genpartIndex() Error or HepMC is empty";
    }
    
#define SaveAndClassifyDChain(thisChain) int Muref=-1;			\
    for (vector<int>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end()&&Muref<(int) eta->size()-1; DChain_iter++ ) { \
      if (*DChain_iter==-1) {						\
	DChain_iter++;							\
	vector<int>::iterator DChain_begin=DChain_iter;			\
	for (; DChain_iter != DChains->end()&&*DChain_iter!=-2&&*DChain_iter!=-1; DChain_iter++ ) ; \
	if (DChain_begin!=DChain_iter)					\
	  if (IstheSameDChain(thisChain,* new vector<int> (DChain_begin,DChain_iter))) { \
	    theSameWithMuon->back()=Muref;				\
	    break;							\
	  }								\
	DChain_iter--;							\
      }									\
      if (*DChain_iter==-2) Muref++;					\
    }									\
    DChains->push_back(-1);						\
    DChains->insert(DChains->end(),thisChain.begin(),thisChain.end());	\
    TheMuonType newtype=classification(thisChain);			\
    if ( find(type.begin(),type.end(),newtype)==type.end() ) {		\
      vector <TheMuonType>::iterator tmp_iter=type.end();		\
      switch (newtype) {						\
          case Others:							\
	    if ( type.empty() ) type.push_back(newtype);		\
	    break;							\
          case NoMuSysHit:						\
	    tmp_iter=find(type.begin(),type.end(),Others);		\
            if ( tmp_iter==type.end() ) break;				\
          default:							\
	    if ( tmp_iter==type.end() ) tmp_iter=find(type.begin(),type.end(),NoMuSysHit); \
	    if ( tmp_iter==type.end() ) type.push_back(newtype);	\
	    else *tmp_iter=newtype;					\
      }									\
    }									\
    if ( (*IsParInHep)[thisChain.front()] ) IsPUTrack->back()=false

    //merge the HepMC and SimTrack Decay Chains
    for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
      if ( !HepMCChains.empty() )
	for (vector< vector<Int_t> >::iterator HepMCChains_iter = HepMCChains.begin(); HepMCChains_iter !=  HepMCChains.end(); ++HepMCChains_iter ) {
	  vector<Int_t> thisChain(HepMCChains_iter->rbegin(),HepMCChains_iter->rend());
	  thisChain.insert(thisChain.end(),SimChains_iter->begin(),SimChains_iter->end());
	  SaveAndClassifyDChain(thisChain);
	}
      else {SaveAndClassifyDChain((*SimChains_iter));}
  }//g4Track loop end
}

bool MuonSysPriEff::IstheSameDChain(vector<int> &Chain1,vector<int> &Chain2) {//whether two chains include each other
  bool ChainIncluded=false;
  vector<int>::iterator Chain1_Particle = Chain1.begin(),Chain2_Particle = Chain2.begin();
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

void MuonSysPriEff::SimTrackDaughtersTree(const SimTrack * thisTrk, TrackingParticleRef tpr)
{
  // Find MC Truth Segment - the offical one use chi2 to match simtrk (MuonIdTruthInfo.cc) and it won't know the decay in flight segment truth
  // The particle type of the hit may differ from the particle type of the SimTrack with id trackId().
  // This happends if the hit was created by a secondary track(e.g. a delta ray) originating from the trackId() and not existing as aseparate SimTrack.
  // ( particle type match notice is from haiyun.teng@cern.ch )
  Bool_t ChainEnd=true, HasMuonHits=false;
  for ( vector<PSimHit>::const_iterator g4Hit_iter=tpr->pSimHit_begin();g4Hit_iter!=tpr->pSimHit_end();g4Hit_iter++ )
    if ( g4Hit_iter->trackId()==thisTrk->trackId() && g4Hit_iter->eventId().rawId()==thisTrk->eventId().rawId() && g4Hit_iter->particleType() == thisTrk->type() )
      HasMuonHits=SimHitsToSegment( *g4Hit_iter);
  
  //To avoid duplicate particle saving
  vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
  if (SavedSimTrk_iter==SavedSimTrk.end()) {
      DChain.push_back(IsParInHep->size());
      IsParHasMuonHits->push_back(HasMuonHits);
      RecordSimTrack(thisTrk)
  }
  else DChain.push_back(FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin()));
  SavedTP.push_back(tpr);

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

void MuonSysPriEff::HepMCParentTree(HepMC::GenParticle *genPar) {
  HepMC::GenVertex *thisVtx = genPar->production_vertex();
  Bool_t ChainEnd=true;
  if (thisVtx) {
      for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	if ((*pgenD)->pdg_id()!=92)  {//Pythia special code for string, we only care about the particles after hadronization
	  ChainEnd=false;
	  vector<HepMC::GenParticle *>::iterator SavedHepPar_iter=find(SavedHepPar.begin(),SavedHepPar.end(),*pgenD);
	  if (SavedHepPar_iter==SavedHepPar.end())
	    {
	      DChain.push_back(IsParInHep->size());
	      IsParHasMuonHits->push_back(false);
	      RecordHepMC((*pgenD))
	    }
	  else DChain.push_back(FindHepMCRecordingPosition(SavedHepPar_iter-SavedHepPar.begin()));
	  HepMCParentTree(*pgenD);
	  DChain.pop_back();
	}
  }
  if (ChainEnd) HepMCChains.push_back(DChain);
}
#endif
//define this as a plug-in
DEFINE_FWK_MODULE(MuonSysPriEff);
