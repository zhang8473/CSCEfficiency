// -*- C++ -*-
//
// Package:    CSCPriEff
// 
/**\class CSCPriEff CSCPriEff.cc UserCode/CSCPriEff/src/CSCPriEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhang Jinzhong
//         Created:  Mon Jun  7 22:19:50 CEST 2010
// $Id: CSCPriEff.cc,v 1.7 2010/10/26 01:20:10 zhangjin Exp $
//
//
//#define LocalRun
//#define FilterOnly
#include "../interface/CSCPriEff.h"
using namespace edm;
using namespace std;
using namespace reco;

#define ReportError(code,MSG) Error.ErrorCode=UInt_t(code);	\
  Error.Run=Info.RUN;			      \
  Error.Event=Info.EVENT;		      \
  ErrorMsg_Tree->Fill();		      \
  cerr<<MSG<<endl;			      

#define GenSimMomentum(ppt,peta,pphi,ppdg) Gen_pdgId->push_back(ppdg);Gen_pt->push_back(ppt);Gen_eta->push_back(peta);Gen_phi->push_back(pphi)

#define GenSimVertex(vtx,vty,vtz,vtt) Gen_vx->push_back(vtx);Gen_vy->push_back(vty);Gen_vz->push_back(vtz);Gen_vt->push_back(vtt)

#define RecordSimTrack(thisTrk) SavedSimTrk.push_back(thisTrk);			\
  IsParInHep->push_back(false);						\
  IsParSavedAsHep.push_back(false);						\
  GenSimMomentum(thisTrk->momentum().pt(),thisTrk->momentum().eta(),thisTrk->momentum().phi(),thisTrk->type());	\
  if (!thisTrk->noVertex()) {						\
    SimVertex thisVtx=SVC[thisTrk->vertIndex()];			\
    GenSimVertex(thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9); \
  }									\
  else {GenSimVertex(0,0,0,0);}

#define RecordHepMC(GenParticle) SavedHepPar.push_back(GenParticle);	\
  IsParInHep->push_back(true);						\
  IsParSavedAsHep.push_back(true);						\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}

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
  if ( type.empty() ) type.push_back(newtype);				\
  else if ( type.size()==1 && type.back()==Others ) type.back()=newtype;	\
  else if ( type.size()==1 && type.back()==NoMuSysHit && newtype!=Others ) type.back()=newtype; \
  else if (find(type.begin(),type.end(),newtype)==type.end() && newtype!=Others && newtype!=NoMuSysHit ) type.push_back(newtype)

#define NtupleMuonSelection (iter->isTrackerMuon()||iter->isGlobalMuon())&&(abs(iter->eta())>0.8)

// ------------ method called on each new Event  ------------

Bool_t CSCPriEff::filter(edm::Event& event, const edm::EventSetup& iSetup) {
  Summarization.Total_Events++;
  
  //Determine if keeps the event
  pt->clear();
  //Set Muons Handle
  Handle<reco::MuonCollection> Muon;
  event.getByLabel("muons", Muon);
  if (!Muon.isValid()) return false;
  reco::MuonCollection const & muons = *Muon;
  
  for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
    if (NtupleMuonSelection) pt->push_back(iter->pt());
  if (pt->size()==0) return false;
#ifndef FilterOnly
  //general event information
  Info.RUN   = event.id ().run ();
  Info.EVENT = event.id ().event();
  Info.LS    = event.luminosityBlock ();
  Info.ORBIT = event.orbitNumber ();
  Info.BX = event.bunchCrossing ();

  ClearVecs();
  //dEdx
  Handle< ValueMap<reco::DeDxData> > dEdxTrackHandle;
  event.getByLabel(dEdxTag, dEdxTrackHandle);
  const ValueMap<reco::DeDxData> dEdxTrack = *(dEdxTrackHandle.product());

  //HepMC Particles
  Handle<edm::HepMCProduct> HepMCH;
  event.getByLabel(HepMCTag, HepMCH);
  HepMC::GenEvent HepGenEvent(*(HepMCH->GetEvent()));
  for (unsigned int i=0;i<4;i++)
    HepMCFourVec[i]=999999.;
  HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent.particles_begin();
  for (;GenParticle_iter != HepGenEvent.particles_end();GenParticle_iter++) { 
    HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
    if (thisVtx) {
      HepMCFourVec[0]=thisVtx->position().x()/10.;
      HepMCFourVec[1]=thisVtx->position().y()/10.;
      HepMCFourVec[2]=thisVtx->position().z()/10.;
      HepMCFourVec[3]=thisVtx->position().t()/299.792458;
      break;
    }
  }
  if (GenParticle_iter == HepGenEvent.particles_end()) {ReportError(7,"No HepMC(Core) Vertex Information")}

  //Reco Tracks
  Handle<reco::TrackCollection> trackCollectionH;
  event.getByLabel(tracksTag,trackCollectionH);
  reco::TrackCollection tC = *trackCollectionH.product();
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(tracksTag,trackCollectionHV);
  //CSC LCT Digis
  edm::Handle<CSCCorrelatedLCTDigiCollection> CSCLCTs;
  try{
    event.getByLabel(CSCDigisTag,CSCLCTs);
  }
  catch (cms::Exception){
    LogError("")<< "Error! Can't get LCTdigi by label. ";
  }

  //Simulated Vertices: the vertexId() is just the position
  Handle<SimVertexContainer> SVCollectionH;
  event.getByLabel("g4SimHits", SVCollectionH);
  SVC = *SVCollectionH.product();
    
  //Simulated Tracks: the trackId() is not the position
  //Handle<SimTrackContainer> STCollectionH;
  //event.getByLabel("g4SimHits", STCollectionH);
  //SimTrackContainer STC = *STCollectionH.product();

  //Tracking Particles (collection of SimTracks and Hits)
  Handle<TrackingParticleCollection> TPCollectionH ;
  //event.getByType(TPCollectionH);
  event.getByLabel(TPInputTag,TPCollectionH);

  //SimToReco Tracks Association
  ESHandle<TrackAssociatorBase> AssociatorByHits;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits", AssociatorByHits);
  //SimToRecoCollection SimToRecoByHits = AssociatorByHits->associateSimToReco(trackCollectionHV,TPCollectionH,&event);
  RecoToSimCollection RecoToSimByHits = AssociatorByHits->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
 
  //Match by chi2
  ESHandle<TrackAssociatorBase> AssociatorByChi2;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByChi2", AssociatorByChi2);
  //SimToRecoCollection SimToRecoByChi2 = AssociatorByChi2->associateSimToReco(trackCollectionHV,TPCollectionH,&event);
  RecoToSimCollection RecoToSimByChi2 = AssociatorByChi2->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
  //Primary Vertex
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(PrimaryVerticesTag.c_str(),recVtxs);
  if (recVtxs.isValid())
    for(reco::VertexCollection::const_iterator v=recVtxs->begin(); v!=recVtxs->end(); ++v) {
      vx->push_back(v->x());	   vxError->push_back(v->xError());
      vy->push_back(v->y());	   vyError->push_back(v->yError());
      vz->push_back(v->z());	   vzError->push_back(v->zError());
    }
  else if (FirstEntry) {ReportError(1,PrimaryVerticesTag.c_str()<<" information is not valid.")}
  for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
    //muon loop begin
    if (NtupleMuonSelection) {
      //muon basic information (pt,eta,phi,charge)
      eta->push_back(iter->eta());
      phi->push_back(iter->phi());
      if (iter->charge()==-1) chargeMinus->push_back(true);
      else chargeMinus->push_back(false);
      //Muon Selectors
      if (iter->isTrackerMuon()) Summarization.Total_TrackerMuons++;
      isTrackerMu->push_back(iter->isTrackerMuon());
      if (iter->isGlobalMuon()) Summarization.Total_GlobalMuon++;
      isGlobalMu->push_back(iter->isGlobalMuon());
      if (!iter->isTrackerMuon()&&iter->isGlobalMuon()) Summarization.Total_GlobalnotTrackerMuon++;
      for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
	SelectorPassed[whichcut]->push_back(muon::isGoodMuon(*iter,OfficialMuonSelectors[whichcut]));
      MySelector->push_back(muon::isGoodMuon(*iter,muon::TMLastStation,1,3,3,3,3,maxChamberDist,maxChamberDistPull,MuonArbitrationType,false,true));
      //Station and Segment Matches
      StationMask->push_back(iter->stationMask(MuonArbitrationType));
      RequiredStationMask->push_back(muon::RequiredStationMask(*iter,maxChamberDist,maxChamberDistPull,MuonArbitrationType));
      //chambers being considered
      vector<const MuonChamberMatch *> Chambers;
      for( vector<MuonChamberMatch>::const_iterator chamberMatch = iter->matches().begin();chamberMatch != iter->matches().end(); chamberMatch++ ) {
	if (chamberMatch->detector()!=MuonSubdetId::CSC) continue;
	Chambers.push_back( &(*chamberMatch) );
      }
      //align save position from inner to outter,from smaller chamber # to big
      if (Chambers.size()>1)
	for( vector<const MuonChamberMatch *>::iterator chamberMatch1 = Chambers.begin();chamberMatch1 != Chambers.end()-1; chamberMatch1++ ) {
	  const CSCDetId ChamberID1( (*chamberMatch1)->id.rawId() );
	  for( vector<const MuonChamberMatch *>::iterator chamberMatch2 = chamberMatch1+1;chamberMatch2 != Chambers.end(); chamberMatch2++ ) {
	    const CSCDetId ChamberID2( (*chamberMatch2)->id.rawId() );
	    if ( ChamberID1.endcap()>ChamberID2.endcap() ) continue;
	    if ( ChamberID1.station()<ChamberID2.station() ) continue;
	    if ( ChamberID1.station()==ChamberID2.station() ) {
	      if ( ChamberID1.ring()%4<ChamberID2.ring()%4 ) continue;
	      if ( ChamberID1.ring()==ChamberID2.ring() && ChamberID1.chamber()<ChamberID2.chamber() ) continue;
	    }
	    swap(*chamberMatch1,*chamberMatch2);
	  }
	}
      for( vector<const MuonChamberMatch *>::iterator chamberMatch = Chambers.begin();chamberMatch != Chambers.end(); chamberMatch++ ) {
	const CSCDetId ChamberID( (*chamberMatch)->id.rawId() );
	Byte_t NumberOfLCTs=0;
	for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator CSCDigi_iter = CSCLCTs->begin(); CSCDigi_iter != CSCLCTs->end(); CSCDigi_iter++) {
	  const CSCDetId& LCTDetID = (*CSCDigi_iter).first;
	  if ( ChamberID.endcap() != LCTDetID.endcap() || ChamberID.station() != LCTDetID.station() || ChamberID.chamber() != LCTDetID.chamber() || ChamberID.ring() != LCTDetID.ring() ) continue;
	  NumberOfLCTs++;
	}
	std::vector<MuonSegmentMatch>::const_iterator segmentMatch = (*chamberMatch)->segmentMatches.begin();
	for( ;segmentMatch != (*chamberMatch)->segmentMatches.end(); segmentMatch++ ) 
	  if (segmentMatch->isMask(MuonSegmentMatch::BestInChamberByDR)) break;
	if ( segmentMatch == (*chamberMatch)->segmentMatches.end() && !NumberOfLCTs && (*chamberMatch)->dist()==999999 ) continue;//drop the chamber with nothing in

	//Chamber Information
	Ring->push_back( ChamberID.station()*10+ChamberID.ring() );
	if ( ChamberID.endcap()==2 ) Ring->back()=-Ring->back();
	Chamber->push_back( ChamberID.chamber() );
	MuonIndex->push_back( eta->size()-1 );
	TrackDistToChamberEdge->push_back( (*chamberMatch)->dist() );
	TrackDistToChamberEdgeErr->push_back( (*chamberMatch)->distErr() );

	XTrack->push_back( (*chamberMatch)->x );
	YTrack->push_back( (*chamberMatch)->y );
	XErrTrack->push_back( (*chamberMatch)->xErr );
	YErrTrack->push_back( (*chamberMatch)->yErr );
	//LCT and segment information
	NumberOfLCTsInChamber->push_back( NumberOfLCTs );
	XSegment->push_back(999999);
	YSegment->push_back(999999);
	XErrSegment->push_back(999999);
	YErrSegment->push_back(999999);
	DRTrackToSegment->push_back(999999);
	DRErrTrackToSegment->push_back(999999);
	NumberOfHitsInSegment->push_back(0);
	IsSegmentOwnedExclusively->push_back( false );
	IsSegmentBestInStationByDR->push_back( false );
	IsSegmentBelongsToTrackByDR->push_back( false );
	IsSegmentBelongsToTrackByCleaning->push_back( false );
	if ( segmentMatch != (*chamberMatch)->segmentMatches.end() ) {
	  if (segmentMatch->cscSegmentRef.isNonnull()) NumberOfHitsInSegment->back() = segmentMatch->cscSegmentRef->specificRecHits().size();
	  else {ReportError(3,"CSC Station"<<(*chamberMatch)->station()<<", Segment Ref not found")}
	  XSegment->back() = segmentMatch->x;
	  YSegment->back() = segmentMatch->y;
	  XErrSegment->back() = segmentMatch->xErr;
	  YErrSegment->back() = segmentMatch->yErr;
	  Float_t DXTrackToSegment=abs( XSegment->back()-XTrack->back() ),DYTrackToSegment=abs( YSegment->back()-YTrack->back() ),DXErrTrackToSegment=sqrt( XErrSegment->back()*XErrSegment->back()+XErrTrack->back()*XErrTrack->back() ),DYErrTrackToSegment=sqrt( YErrSegment->back()*YErrSegment->back()+YErrTrack->back()*YErrTrack->back() );
	  DRTrackToSegment->back() = sqrt(DXTrackToSegment*DXTrackToSegment+DYTrackToSegment*DYTrackToSegment);
	  DRErrTrackToSegment->back() = sqrt(DXTrackToSegment*DXTrackToSegment*DXErrTrackToSegment*DXErrTrackToSegment+DYTrackToSegment*DYTrackToSegment*DYErrTrackToSegment*DYErrTrackToSegment)/DRTrackToSegment->back();
	  UInt_t CurrentChamber=XTrack->size()-1;
	  IsSegmentOwnedExclusively->back()=true;
	  for (UInt_t chamber=0; chamber<CurrentChamber&&MuonIndex->at(chamber)!=MuonIndex->back(); chamber++)
	    if ( IsTheSameSegment(chamber,CurrentChamber) ) {
	      IsSegmentOwnedExclusively->back()=false;
	      IsSegmentOwnedExclusively->at(chamber)=false;
	    }
	  IsSegmentBestInStationByDR->back() = segmentMatch->isMask(MuonSegmentMatch::BestInStationByDR);
	  IsSegmentBelongsToTrackByDR->back() = segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR);
	  IsSegmentBelongsToTrackByCleaning->back() = segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning);
	}
      }//chamber loop end
      //miscellany
      Vertex_x->push_back(iter->vertex().X());
      Vertex_y->push_back(iter->vertex().Y());
      Vertex_z->push_back(iter->vertex().Z());
      isoR03sumPt->push_back(iter->isolationR03().sumPt);
      isoR03emEt->push_back(iter->isolationR03().emEt);
      isoR03hadEt->push_back(iter->isolationR03().hadEt);
      isoR03hoEt->push_back(iter->isolationR03().hoEt);
      isoR03nJets->push_back(iter->isolationR03().nJets);
      isoR03nTracks->push_back(iter->isolationR03().nTracks);
      isoR05sumPt->push_back(iter->isolationR05().sumPt);
      isoR05emEt->push_back(iter->isolationR05().emEt);
      isoR05hadEt->push_back(iter->isolationR05().hadEt);
      isoR05hoEt->push_back(iter->isolationR05().hoEt);
      isoR05nJets->push_back(iter->isolationR05().nJets);
      isoR05nTracks->push_back(iter->isolationR05().nTracks);
      isoemVetoEt->push_back(iter->isolationR03().emVetoEt);
      isohadVetoEt->push_back(iter->isolationR03().hadVetoEt);
      isohoVetoEt->push_back(iter->isolationR03().hoVetoEt);
      TrkKink->push_back(iter->combinedQuality().trkKink);
      GlbKink->push_back(iter->combinedQuality().glbKink);
      TrkRelChi2->push_back(iter->combinedQuality().trkRelChi2);
      CaloE_emMax->push_back(iter->calEnergy().emMax);
      CaloE_emS9->push_back(iter->calEnergy().emS9);
      CaloE_emS25->push_back(iter->calEnergy().emS25);
      CaloE_hadMax->push_back(iter->calEnergy().hadMax);
      CaloE_hadS9->push_back(iter->calEnergy().hadS9);
      Calo_emPos_R->push_back(iter->calEnergy().ecal_position.R());
      Calo_emPos_eta->push_back(iter->calEnergy().ecal_position.eta());
      Calo_emPos_phi->push_back(iter->calEnergy().ecal_position.phi());
      Calo_hadPos_R->push_back(iter->calEnergy().hcal_position.R());
      Calo_hadPos_eta->push_back(iter->calEnergy().hcal_position.eta());
      Calo_hadPos_phi->push_back(iter->calEnergy().hcal_position.phi());
      //match to RecoTrk
      TrackRef innertrack = iter->innerTrack();
      dEdx->push_back(dEdxTrack[innertrack].dEdx());
      dEdxError->push_back(dEdxTrack[innertrack].dEdxError());
      dEdx_numberOfSaturatedMeasurements->push_back(dEdxTrack[innertrack].numberOfSaturatedMeasurements());
      dEdx_numberOfMeasurements->push_back(dEdxTrack[innertrack].numberOfMeasurements());
      theSameWithMuon->push_back(-1);
      InnerTrack_nValidTrackerHits->push_back(innertrack->hitPattern().numberOfValidTrackerHits());
      InnerTrack_nValidPixelHits->push_back(innertrack->hitPattern().numberOfValidPixelHits());
      InnerTrack_nLostPixelHits->push_back(innertrack->hitPattern().numberOfLostPixelHits());
      InnerTrack_nLostTrackerHits->push_back(innertrack->hitPattern().numberOfLostTrackerHits());
      InnerTrack_chi2->push_back(innertrack->chi2());
      InnerTrack_ndof->push_back(innertrack->ndof());
      //Match to Tracking Particle Collection
      TrkParticles_pt->push_back(0);
      TrkParticles_eta->push_back(0);
      TrkParticles_phi->push_back(0);
      TrkParticles_pdgId->push_back(0);
      TrkParticles_charge->push_back(-100.0);
      SharedHitsRatio->push_back(-100.0);
      MCMatchChi2->push_back(-100.0);
      MuonType->push_back(0);
      NumMisMatch++;
      unsigned int innertrack_pos;  DChains->push_back(-2);
      vector <TheMuonType> type; FirstLayer.clear();LastLayer.clear();
      for(innertrack_pos=0; innertrack_pos<tC.size(); innertrack_pos++)  //recotrk loop begin
	if (innertrack == TrackRef(trackCollectionH,innertrack_pos) ) {
	  RefToBase<Track> trk(trackCollectionHV, innertrack_pos);
	  if (InnerTrack_nValidTrackerHits->back()<minTrackHits||RecoToSimByHits.find(trk) == RecoToSimByHits.end()||RecoToSimByChi2.find(trk) == RecoToSimByChi2.end()) break;
	  pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trk].front();
	  vector<pair<TrackingParticleRef, double> > TPCByChi2=RecoToSimByChi2[trk];
	  vector<pair<TrackingParticleRef, double> >::iterator TPCByChi2_iter=TPCByChi2.begin();
	  for (;TPCByChi2_iter!=TPCByChi2.end();TPCByChi2_iter++)
	    if (BestMatch.first==TPCByChi2_iter->first) break;
	  if (TPCByChi2_iter==TPCByChi2.end()) break;
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
		      vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
		      if (SavedSimTrk_iter==SavedSimTrk.end())
			{
			  DChain.push_back(IsParInHep->size());
			  IsParHasMuonHits->push_back(false);
	      		  RecordSimTrack(thisTrk)
			}
		      else DChain.push_back(FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin()));
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
	      else {ReportError(5,"Code is Wrong")}
	      HepMC::GenEvent::particle_iterator genPar = HepGenEvent.particles_begin();
		for (int count=1; count<thisTrk->genpartIndex()&&genPar != HepGenEvent.particles_end(); count++ )
		  genPar++;
		if (genPar != HepGenEvent.particles_end()) HepMCParentTree(*genPar);
		else {ReportError(6,"genpartIndex() Error or HepMC is empty")}
	    }
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
	  Long64_t TypeRecord=0;
	  for (vector <TheMuonType>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	    if (*type_iter) TypeRecord=TypeRecord*100+Long64_t(*type_iter);
	  MuonType->back()=TypeRecord;
	  break;
	}//recotrk loop end
      if ( eta->back()>0 )
	for ( map<Int_t,CSCSimHitInfo>::iterator SimSeg=FirstLayer.begin(); SimSeg != FirstLayer.end(); SimSeg++) {
	  RingSimSegment->push_back( Char_t(SimSeg->first/100) );
	  ChamberSimSegment->push_back( abs(SimSeg->first)%100 );
	  MuonIndexSimSegment->push_back( eta->size()-1 );
	  XSimSegment->push_back( (LastLayer[SimSeg->first].X+SimSeg->second.X)/2. );
	  YSimSegment->push_back( (LastLayer[SimSeg->first].Y+SimSeg->second.Y)/2. );
	}
      else
	for ( map<Int_t,CSCSimHitInfo>::reverse_iterator SimSeg=FirstLayer.rbegin(); SimSeg != FirstLayer.rend(); SimSeg++) {
	  RingSimSegment->push_back( Char_t(SimSeg->first/100) );
	  ChamberSimSegment->push_back( abs(SimSeg->first)%100 );
	  MuonIndexSimSegment->push_back( eta->size()-1 );
	  XSimSegment->push_back( (LastLayer[SimSeg->first].X+SimSeg->second.X)/2. );
	  YSimSegment->push_back( (LastLayer[SimSeg->first].Y+SimSeg->second.Y)/2. );
	}
    }//muon loop end
#ifdef LocalRun
  Bool_t Print=false;
  /*for (vector<Int_t>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end(); DChain_iter++ )
    if (*DChain_iter==-2) 
      for (vector<Int_t>::iterator DChain_iter2 = DChain_iter+1; DChain_iter2 != DChains->end()&&*DChain_iter2!=-2; DChain_iter2++ )
	  for (vector<Int_t>::iterator DChain_iter3 = DChains->begin(); DChain_iter3 != DChain_iter; DChain_iter3++ )
	    if (*DChain_iter2==*DChain_iter3&&*DChain_iter2>0) Print=true;
  int muref=-1;
  for (vector<Int_t>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end(); DChain_iter++ ) {
    if (*DChain_iter==-2) muref++;
    if ((*MuonType)[muref]==21) {Print=true;break;}
    }*/
  if (true) {
    cout<<"Run "<<Info.RUN<<"; Event "<<Info.EVENT<<endl;
    for (unsigned int parref = 0; parref<Gen_pdgId->size(); parref++)
      printf("%d: %d (%3.2fcm,%3.2fcm,%3.2fcm,%3.2fns) ---(%3.3f,%3.3f,%3.3f,%3.3f)GeV/c---> %s;\n",parref,Gen_pdgId->at(parref),Gen_vx->at(parref),Gen_vy->at(parref),Gen_vz->at(parref),Gen_vt->at(parref),Gen_pt->at(parref)*cosh(Gen_eta->at(parref)),Gen_pt->at(parref)*cos(Gen_phi->at(parref)),Gen_pt->at(parref)*sin(Gen_phi->at(parref)),Gen_pt->at(parref)*sinh(Gen_eta->at(parref)),IsParInHep->at(parref)?"HEPMC":"SimTrk");
    unsigned int Muref=0;
    for (vector<Int_t>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end(); DChain_iter++ )
      {
	if (*DChain_iter==-1) printf("\n");
	else if (*DChain_iter==-2) {
	  printf("\n Mu (%d): (%3.2fcm,%3.2fcm,%3.2fcm) --- (%3.2fGeV/c,%3.2f,%3.2f) ---> SameWith: Mu %d, Type",Muref,Vertex_x->at(Muref),Vertex_y->at(Muref),Vertex_z->at(Muref),pt->at(Muref)*cosh(eta->at(Muref)),eta->at(Muref),phi->at(Muref),theSameWithMuon->at(Muref));
	  cout<<MuonType->at(Muref)<<endl;
	  Muref++;
	}
	else printf("%d --> ", *DChain_iter);
      }
    printf("\n----------------------------------------------------\n");
  }
#endif
  FirstEntry=false;
  Tracks_Tree->Fill();
#endif
  return true;
}

//whether two chains include each other
Bool_t CSCPriEff::IstheSameDChain(vector<int> &ThisChain,vector<int> &AnotherChain)
{
  Bool_t ChainIncluded=false;
  vector<int>::iterator ThisChain_Particle = ThisChain.begin(),AnotherChain_Particle = AnotherChain.begin();
  for (; AnotherChain_Particle != AnotherChain.end()&&ThisChain_Particle != ThisChain.end(); AnotherChain_Particle++) {
    if (ChainIncluded&&*ThisChain_Particle!=*AnotherChain_Particle) ChainIncluded=false;
    if (ThisChain.front() == *AnotherChain_Particle) {
      ChainIncluded=true;
      ThisChain_Particle = ThisChain.begin();
    }
    if (ChainIncluded) ThisChain_Particle++;
  }
  if (!ChainIncluded) {
    AnotherChain_Particle = AnotherChain.begin(); ThisChain_Particle = ThisChain.begin();
    for (; AnotherChain_Particle != AnotherChain.end()&&ThisChain_Particle != ThisChain.end(); ThisChain_Particle++) {
      if (ChainIncluded&&*ThisChain_Particle!=*AnotherChain_Particle) ChainIncluded=false;
      if (AnotherChain.front() == *ThisChain_Particle) {
	AnotherChain_Particle = AnotherChain.begin();
	ChainIncluded=true;
      }
      if (ChainIncluded) AnotherChain_Particle++;
    }
  }
  return ChainIncluded;
}

//particle classifiction based on pdg_id
CSCPriEff::ParticleType CSCPriEff::ParticleCata(int pid)
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
CSCPriEff::TheMuonType CSCPriEff::classification(vector<Int_t> &Chain)
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
      float deltaRPhi2=((*Gen_vx)[MuPos]-HepMCFourVec[0])*((*Gen_vx)[MuPos]-HepMCFourVec[0])+((*Gen_vy)[MuPos]-HepMCFourVec[1])*((*Gen_vy)[MuPos]-HepMCFourVec[1]);
      float deltaZ=abs((*Gen_vz)[MuPos]-HepMCFourVec[2]);
      //float RPhi2=((*Gen_vx)[MuPos])*((*Gen_vx)[MuPos])+((*Gen_vy)[MuPos])*((*Gen_vy)[MuPos]);
      //float Z=abs((*Gen_vz)[MuPos]);
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

void CSCPriEff::SimTrackDaughtersTree(const SimTrack * thisTrk, TrackingParticleRef tpr)
{
  // Find MC Truth Segment - the offical one use chi2 to match simtrk (MuonIdTruthInfo.cc) and it won't know the decay in flight segment truth
  // The particle type of the hit may differ from the particle type of the SimTrack with id trackId().
  // This happends if the hit was created by a secondary track(e.g. a delta ray) originating from the trackId() and not existing as aseparate SimTrack.
  // ( particle type match notice is from haiyun.teng@cern.ch )
  Bool_t ChainEnd=true, HasMuonHits=false;
  for ( vector<PSimHit>::const_iterator g4Hit_iter=tpr->pSimHit_begin();g4Hit_iter!=tpr->pSimHit_end();g4Hit_iter++ )
    if ( g4Hit_iter->trackId()==thisTrk->trackId() && g4Hit_iter->eventId().rawId()==thisTrk->eventId().rawId() && g4Hit_iter->particleType() == thisTrk->type() )      HasMuonHits=SimHitToSegment( *g4Hit_iter);
  
  //To avoid duplicate particle saving
  vector<const SimTrack *>::iterator SavedSimTrk_iter=find(SavedSimTrk.begin(),SavedSimTrk.end(),thisTrk);
  if (SavedSimTrk_iter==SavedSimTrk.end()) {
      DChain.push_back(IsParInHep->size());
      IsParHasMuonHits->push_back(HasMuonHits);
      RecordSimTrack(thisTrk)
  }
  else DChain.push_back(FindSimTrackRecordingPosition(SavedSimTrk_iter-SavedSimTrk.begin()));
 
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
//changed from MuonIdTruthInfo.cc: that one saves the layer but not chamber infomation, the return value is whether or not it is hit in muon system
Bool_t CSCPriEff::SimHitToSegment(const PSimHit& hit) {
   // find the hit position projection at the reference surface of the layer:
   // first get entry and exit point of the hit in the global coordinates, then
   // get local coordinates of these points wrt the chamber and then find the
   // projected X-Y coordinates
  DetId DetectorId( hit.detUnitId() );
  if ( DetectorId.det() != DetId::Muon ) return false;
  if ( DetectorId.subdetId() != MuonSubdetId::CSC ) return true;
  const CSCDetId CSCChamberID( hit.detUnitId() );
  Byte_t thisLayer=CSCChamberID.layer();
  Int_t ThisChamber=CSCChamberID.station()*10+CSCChamberID.ring();
  //reject non-existing CSC Chambers (incorrect geometry) In ME42, only ME+42/9-13 exist
  if ( ThisChamber==-42 ) return true;
  if ( ThisChamber==42 && ( CSCChamberID.chamber()<9 || CSCChamberID.chamber()>13 ) ) return true;

  Bool_t SaveToFirstLayer=false,SaveToLastLayer=false;
  ThisChamber=ThisChamber*100+CSCChamberID.chamber();
  if ( CSCChamberID.endcap()==2 ) ThisChamber=-ThisChamber;
  map<Int_t,CSCSimHitInfo>::iterator thisChamberFirstLayer=FirstLayer.find(ThisChamber),thisChamberLastLayer=LastLayer.find(ThisChamber);
  if ( thisChamberFirstLayer!=FirstLayer.end() ) {
    if ( thisLayer<thisChamberFirstLayer->second.Layer) SaveToFirstLayer=true;
  }
  else SaveToFirstLayer=true;
  if ( thisChamberLastLayer!=LastLayer.end() ) {
    if ( thisLayer>thisChamberLastLayer->second.Layer) SaveToLastLayer=true;
  }
  else SaveToLastLayer=true;
  if ( !SaveToFirstLayer && !SaveToLastLayer ) return true;
  //  cerr<<"Layer: "<<CSCChamberID.layer()<<endl;
  Local3DPoint entryPoint = hit.entryPoint(), exitPoint = hit.exitPoint();
  LocalVector direction = exitPoint - entryPoint;
  if ( fabs(direction.z()) > 0.001) {
    LocalPoint projection = entryPoint - direction*(entryPoint.z()/direction.z());
    if ( fabs(projection.z()) > 0.001 ) {ReportError(2,"z coordinate of the hit projection must be zero and it's not!")}
    CSCSimHitInfo thisCSCSimHit={projection.x(),projection.y(),thisLayer};
    if (SaveToFirstLayer) FirstLayer[ThisChamber]=thisCSCSimHit;
    if (SaveToLastLayer) LastLayer[ThisChamber]=thisCSCSimHit;
  }
  return true;
}

void CSCPriEff::HepMCParentTree(HepMC::GenParticle *genPar) {
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

void CSCPriEff::ClearVecs() {
//clear all the other vectors except pt
//  MyMuons->clear();
  eta->clear();  phi->clear();  chargeMinus->clear();
  isGlobalMu->clear();  isTrackerMu->clear();
  Vertex_x->clear();  Vertex_y->clear();  Vertex_z->clear();
  isoR03sumPt->clear();  isoR03emEt->clear();  isoR03hadEt->clear();
  isoR03hoEt->clear();  isoR03nJets->clear();  isoR03nTracks->clear();
  isoR05sumPt->clear();  isoR05emEt->clear();  isoR05hadEt->clear();
  isoR05hoEt->clear();  isoR05nJets->clear();  isoR05nTracks->clear();
  isoemVetoEt->clear();  isohadVetoEt->clear();  isohoVetoEt->clear();
  TrkKink->clear();  GlbKink->clear();  TrkRelChi2->clear();
  CaloE_emMax->clear();  CaloE_emS9->clear();  CaloE_emS25->clear();
  CaloE_hadMax->clear();  CaloE_hadS9->clear();
  Calo_emPos_R->clear();  Calo_emPos_eta->clear();  Calo_emPos_phi->clear();
  Calo_hadPos_R->clear();  Calo_hadPos_eta->clear();  Calo_hadPos_phi->clear();

  dEdx->clear();  dEdxError->clear();
  dEdx_numberOfSaturatedMeasurements->clear();  dEdx_numberOfMeasurements->clear();

  InnerTrack_nValidTrackerHits->clear(); InnerTrack_nValidPixelHits->clear();  
  InnerTrack_nLostTrackerHits->clear(); InnerTrack_nLostPixelHits->clear(); 
  InnerTrack_chi2->clear(); InnerTrack_ndof->clear();
  //Chamber
  TrackDistToChamberEdge->clear();  TrackDistToChamberEdgeErr->clear();
  XTrack->clear();YTrack->clear();  XErrTrack->clear();YErrTrack->clear();
  Ring->clear();Chamber->clear();MuonIndex->clear();
  NumberOfLCTsInChamber->clear();
  //Segment
  XSegment->clear();  YSegment->clear();
  XErrSegment->clear();  YErrSegment->clear();
  DRTrackToSegment->clear();  DRErrTrackToSegment->clear();
  IsSegmentOwnedExclusively->clear();  IsSegmentBestInStationByDR->clear();
  IsSegmentBelongsToTrackByDR->clear();  IsSegmentBelongsToTrackByCleaning->clear();
  NumberOfHitsInSegment->clear();
  StationMask->clear();RequiredStationMask->clear();
  //SimSegment
  RingSimSegment->clear(); ChamberSimSegment->clear();  MuonIndexSimSegment->clear();
  XSimSegment->clear(); YSimSegment->clear();
  //Primary Vertex
  vx->clear();   vxError->clear();
  vy->clear();   vyError->clear();
  vz->clear();   vzError->clear();
  //Hep/Sim Particles
  Gen_pt->clear();  Gen_eta->clear();  Gen_phi->clear(); Gen_pdgId->clear();
  Gen_vx->clear();  Gen_vy->clear();  Gen_vz->clear();  Gen_vt->clear();
  
  MySelector->clear();
  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    SelectorPassed[whichcut]->clear(); 

  TrkParticles_pt->clear();  TrkParticles_eta->clear(); TrkParticles_phi->clear();
  TrkParticles_pdgId->clear(); TrkParticles_charge->clear(); 
  SharedHitsRatio->clear(); MCMatchChi2->clear();

  SavedSimTrk.clear(); SavedHepPar.clear(); IsParSavedAsHep.clear();
  IsParInHep->clear(); IsParHasMuonHits->clear();
  DChains->clear(); theSameWithMuon->clear(); MuonType->clear();
}

#define MakeVecBranch(Name,Var,Type) Var=new vector<Type>();Tracks_Tree->Branch(Name,&Var)

CSCPriEff::CSCPriEff(const edm::ParameterSet& pset) {
//---- Get the input parameters
  FileName = pset.getUntrackedParameter<string>("FileName","MuonTracks.root");
  maxChamberDist = pset.getUntrackedParameter<double>("maxChamberDist",-3.);
  maxChamberDistPull = pset.getUntrackedParameter<double>("maxChamberDistPull",-3.);
  vector<string> Cuts;
  Cuts = pset.getUntrackedParameter< vector<string> >("StandardMuonCuts",Cuts);//default is empty
  string MuonArbitrationTypeStr;
  MuonArbitrationTypeStr = pset.getUntrackedParameter<string>("MuonArbitrationType","SegmentArbitration");
  MuonArbitrationType=MuonArbitrationTypeFromString(MuonArbitrationTypeStr.c_str());
  minTrackHits = pset.getUntrackedParameter<uint>("minTrackHits",3);
  PrimaryVerticesTag = pset.getUntrackedParameter<string>("PrimaryVertices","offlinePrimaryVerticesWithBS");//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
  const InputTag tracksTag_default("generalTracks");
  tracksTag = pset.getUntrackedParameter<InputTag>("tracksTag",tracksTag_default);
  const InputTag TPInputTag_default("mergedtruth","MergedTrackTruth");
  TPInputTag = pset.getUntrackedParameter<InputTag>("TPInputTag",TPInputTag_default);
  CSCDigisTag = pset.getParameter<InputTag>("CSCDigisTag");
  dEdxTag = pset.getUntrackedParameter<string>("dEdxTag","dedxHarmonic2");//Other options already available in RECO files are dedxMedian and dedxTruncated40. 
  HepMCTag = pset.getUntrackedParameter<string>("HepMCTag","generator");
//Build a root file and a TTree
  file=new TFile(FileName.c_str(),"RECREATE");
  Tracks_Tree = new TTree ("MuTrkCand","TrksCandidates") ;
  Summarization_Tree = new TTree ("Summary","Summary") ;
  ErrorMsg_Tree = new TTree ("Errors","Errors") ;
  Tracks_Tree->SetCircular(500000);//the max events in a single root file
  ErrorMsg_Tree->SetCircular(500000);

//Build Branches
  //General event information
  Tracks_Tree->Branch("Event_Info",&Info.RUN,"RUN/l:EVENT:LumiBlock:ORBIT:BunchCrossing");
  //Tracker Muon Information
  MakeVecBranch("pt",pt,Float_t);  MakeVecBranch("eta",eta,Float_t);  MakeVecBranch("phi",phi,Float_t);
  MakeVecBranch("chargeMinus",chargeMinus,Bool_t); MakeVecBranch("isGlobalMu",isGlobalMu,Bool_t);  MakeVecBranch("isTrackerMu",isTrackerMu,Bool_t);
  MakeVecBranch("Vertex_x",Vertex_x,Float_t);  MakeVecBranch("Vertex_y",Vertex_y,Float_t);  MakeVecBranch("Vertex_z",Vertex_z,Float_t);
  MakeVecBranch("isoR03sumPt",isoR03sumPt,Float_t);  MakeVecBranch("isoR03emEt",isoR03emEt,Float_t);  MakeVecBranch("isoR03hadEt",isoR03hadEt,Float_t);
  MakeVecBranch("isoR03hoEt",isoR03hoEt,Float_t);  MakeVecBranch("isoR03nJets",isoR03nJets,Float_t);  MakeVecBranch("isoR03nTracks",isoR03nTracks,Float_t);
  MakeVecBranch("isoR05sumPt",isoR05sumPt,Float_t);  MakeVecBranch("isoR05emEt",isoR05emEt,Float_t);  MakeVecBranch("isoR05hadEt",isoR05hadEt,Float_t);
  MakeVecBranch("isoR05hoEt",isoR05hoEt,Float_t);  MakeVecBranch("isoR05nJets",isoR05nJets,Float_t);  MakeVecBranch("isoR05nTracks",isoR05nTracks,Float_t);
  MakeVecBranch("isoemVetoEt",isoemVetoEt,Float_t);  MakeVecBranch("isohadVetoEt",isohadVetoEt,Float_t);  MakeVecBranch("isohoVetoEt",isohoVetoEt,Float_t);
  MakeVecBranch("TrkKink",TrkKink,Float_t); MakeVecBranch("GlbKink",GlbKink,Float_t); MakeVecBranch("TrkRelChi2",TrkRelChi2,Float_t);
  MakeVecBranch("CaloE_emMax",CaloE_emMax,Float_t);  MakeVecBranch("CaloE_emS9",CaloE_emS9,Float_t);  MakeVecBranch("CaloE_emS25",CaloE_emS25,Float_t);
  MakeVecBranch("CaloE_hadMax",CaloE_hadMax,Float_t);  MakeVecBranch("CaloE_hadS9",CaloE_hadS9,Float_t);
  MakeVecBranch("Calo_emPos_R",Calo_emPos_R,Float_t);  MakeVecBranch("Calo_emPos_eta",Calo_emPos_eta,Float_t);  MakeVecBranch("Calo_emPos_phi",Calo_emPos_phi,Float_t);
  MakeVecBranch("Calo_hadPos_R",Calo_hadPos_R,Float_t);  MakeVecBranch("Calo_hadPos_eta",Calo_hadPos_eta,Float_t);  MakeVecBranch("Calo_hadPos_phi",Calo_hadPos_phi,Float_t);
  //InnerTrack
  MakeVecBranch("dEdx",dEdx,Float_t);
  MakeVecBranch("dEdxError",dEdxError,Float_t);
  MakeVecBranch("dEdx_numberOfSaturatedMeasurements",dEdx_numberOfSaturatedMeasurements,Int_t);
  MakeVecBranch("dEdx_numberOfMeasurements",dEdx_numberOfMeasurements,Int_t);
  MakeVecBranch("InnerTrack_nValidTrackerHits",InnerTrack_nValidTrackerHits,UInt_t);
  MakeVecBranch("InnerTrack_nValidPixelHits",InnerTrack_nValidPixelHits,UInt_t);
  MakeVecBranch("InnerTrack_nLostTrackerHits",InnerTrack_nLostTrackerHits,UInt_t);
  MakeVecBranch("InnerTrack_nLostPixelHits",InnerTrack_nLostPixelHits,UInt_t);
  MakeVecBranch("InnerTrack_chi2",InnerTrack_chi2,Float_t);
  MakeVecBranch("InnerTrack_ndof",InnerTrack_ndof,UInt_t);
  //Muon Chamber Information
  MakeVecBranch("TrackDistToChamberEdge",TrackDistToChamberEdge,Float_t);
  MakeVecBranch("TrackDistToChamberEdgeErr",TrackDistToChamberEdgeErr,Float_t);
  MakeVecBranch("XTrack",XTrack,Float_t);
  MakeVecBranch("YTrack",YTrack,Float_t);
  MakeVecBranch("XErrTrack",XErrTrack,Float_t);
  MakeVecBranch("YErrTrack",YErrTrack,Float_t);
  MakeVecBranch("Ring",Ring,Char_t);
  MakeVecBranch("Chamber",Chamber,Byte_t);
  MakeVecBranch("MuonIndex",MuonIndex,Byte_t);
  MakeVecBranch("NumberOfLCTsInChamber",NumberOfLCTsInChamber,Byte_t);
  //Muon Segment Information
  MakeVecBranch("XSegment",XSegment,Float_t);
  MakeVecBranch("YSegment",YSegment,Float_t);
  MakeVecBranch("XErrSegment",XErrSegment,Float_t);
  MakeVecBranch("YErrSegment",YErrSegment,Float_t);
  MakeVecBranch("DRTrackToSegment",DRTrackToSegment,Float_t);
  MakeVecBranch("DRErrTrackToSegment",DRErrTrackToSegment,Float_t);
  MakeVecBranch("IsSegmentOwnedExclusively",IsSegmentOwnedExclusively,Bool_t);
  MakeVecBranch("IsSegmentBestInStationByDR",IsSegmentBestInStationByDR,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByDR",IsSegmentBelongsToTrackByDR,Bool_t);
  MakeVecBranch("IsSegmentBelongsToTrackByCleaning",IsSegmentBelongsToTrackByCleaning,Bool_t);
  MakeVecBranch("NumberOfHitsInSegment",NumberOfHitsInSegment,Byte_t);
  MakeVecBranch("StationMask",StationMask,UInt_t);MakeVecBranch("RequiredStationMask",RequiredStationMask,UInt_t); 
  //Muon SimSegment Information
  MakeVecBranch("RingSimSegment",RingSimSegment,Char_t);
  MakeVecBranch("ChamberSimSegment",ChamberSimSegment,Byte_t);
  MakeVecBranch("MuonIndexSimSegment",MuonIndexSimSegment,Byte_t);
  MakeVecBranch("XSimSegment",XSimSegment,Float_t);
  MakeVecBranch("YSimSegment",YSimSegment,Float_t);
  //PV
  MakeVecBranch("PVx",vx,Float_t);  MakeVecBranch("PVy",vy,Float_t);  MakeVecBranch("PVz",vz,Float_t);
  MakeVecBranch("PVxError",vxError,Float_t);  MakeVecBranch("PVyError",vyError,Float_t);  MakeVecBranch("PVzError",vzError,Float_t);
  //Generation Level Muons
  MakeVecBranch("Gen_pt",Gen_pt,Float_t);  MakeVecBranch("Gen_eta",Gen_eta,Float_t);  MakeVecBranch("Gen_phi",Gen_phi,Float_t);   MakeVecBranch("Gen_pdgId",Gen_pdgId,Int_t);
  MakeVecBranch("Gen_vx",Gen_vx,Float_t);  MakeVecBranch("Gen_vy",Gen_vy,Float_t);  MakeVecBranch("Gen_vz",Gen_vz,Float_t);  MakeVecBranch("Gen_vt",Gen_vt,Float_t);
  MakeVecBranch("IsParInHep",IsParInHep,Bool_t);
  MakeVecBranch("IsParHasMuonHits",IsParHasMuonHits,Bool_t);

  Tracks_Tree->Branch("HepMCVertex",HepMCFourVec,"HepMCVertex[4]/F");

  //Simulated Tracks
  MakeVecBranch("TrkParticles_pt",TrkParticles_pt,Float_t);  MakeVecBranch("TrkParticles_eta",TrkParticles_eta,Float_t);  MakeVecBranch("TrkParticles_phi",TrkParticles_phi,Float_t); MakeVecBranch("TrkParticles_charge",TrkParticles_charge,Int_t);
  MakeVecBranch("TrkParticles_pdgId",TrkParticles_pdgId,Int_t);  MakeVecBranch("SharedHitsRatio",SharedHitsRatio,Double_t);
  MakeVecBranch("MCMatchChi2",MCMatchChi2,Double_t);
  MakeVecBranch("DChains",DChains,Int_t); MakeVecBranch("theSameWithMuon",theSameWithMuon,Int_t);
  MakeVecBranch("MuonType",MuonType,Long64_t);
  //Cuts
  num_Cuts=Cuts.size();
  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++) 
    if (muon::selectionTypeFromString(Cuts[whichcut])!=(muon::SelectionType)-1) {
      MakeVecBranch(Cuts[whichcut].c_str(),SelectorPassed[whichcut],Bool_t); 
      OfficialMuonSelectors.push_back(muon::selectionTypeFromString(Cuts[whichcut]));
    }
  num_Cuts=OfficialMuonSelectors.size();
  MakeVecBranch("MySelector",MySelector,Bool_t);
  //Summary
  ErrorMsg_Tree->Branch("ErrorMsg",&Error.ErrorCode,"ErrorCode/l:RunNum:EventNum");
  Summarization_Tree->Branch("Summarization",&Summarization.Total_Events,"Total_Events/l:Total_TrackerMuons:Total_GlobalMuon:Total_GlobalnotTrackerMuon");
  Summarization.Total_Events=0;Summarization.Total_TrackerMuons=0;Summarization.Total_GlobalMuon=0;Summarization.Total_GlobalnotTrackerMuon=0;
  FirstEntry=true;
  NumMisMatch=0;
}

Muon::ArbitrationType CSCPriEff::MuonArbitrationTypeFromString( const std::string &label ) {
  struct MuonArbitrationTypeStringToEnum { const char *label; reco::Muon::ArbitrationType value; };
  static MuonArbitrationTypeStringToEnum MuonArbitrationTypeStringToEnumMap[] = {
    {"NoArbitration",Muon::NoArbitration},
    {"SegmentArbitration",Muon::SegmentArbitration},
    {"SegmentAndTrackArbitration",Muon::SegmentAndTrackArbitration},
    {"SegmentAndTrackArbitrationCleaned",Muon::SegmentAndTrackArbitrationCleaned},
    { 0, (Muon::ArbitrationType)-1 }
  };
  Muon::ArbitrationType value = (Muon::ArbitrationType)-1;
  Bool_t found = false;
  for(int i = 0; MuonArbitrationTypeStringToEnumMap[i].label && (! found); ++i)
    if (! strcmp(label.c_str(), MuonArbitrationTypeStringToEnumMap[i].label)) {
      found = true;
      value = MuonArbitrationTypeStringToEnumMap[i].value;
    }
  
  // in case of unrecognized selection type
  if (! found) throw cms::Exception("MuonSelectorError") << label << " is not a recognized SelectionType";
  return value;
}

CSCPriEff::~CSCPriEff()
{
  printf("MisRate:%f%%\n",NumMisMatch/(float) (Summarization.Total_TrackerMuons+Summarization.Total_GlobalnotTrackerMuon)*100);
  //Summarizations
  Summarization_Tree->Fill();
  file->Write();
  file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCPriEff);
