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
// $Id: CSCPriEff.cc,v 1.2 2010/07/12 14:40:37 zhangjin Exp $
//
//
//#define LocalRun
//#define FilterOnly
#include "../interface/CSCPriEff.h"

using namespace edm;
using namespace std;
using namespace reco;

#define ReportError(code,MSG) Error.ErrorCode=code;    \
  Error.Run=Info.RUN;			      \
  Error.Event=Info.EVENT;		      \
  ErrorMsg_Tree->Fill();		      \
  cerr<<MSG<<endl;			      

#define trackIdLink(TrkID)  vector<SimTrack>::iterator Trk_iter = STC.begin(); \
  for (; Trk_iter != STC.end(); ++Trk_iter )	\
    if ((int) Trk_iter->trackId() == (int) TrkID) break;		\
  if (Trk_iter==STC.end()) {ReportError(2,"parentIndex/trackId Error")}	\
  thisTrk=&*Trk_iter;

#define GenSimMomentum(ppt,peta,pphi,ppdg) Gen_pdgId->push_back(ppdg);Gen_pt->push_back(ppt);Gen_eta->push_back(peta);Gen_phi->push_back(pphi)

#define GenSimVertex(vtx,vty,vtz,vtt) Gen_vx->push_back(vtx);Gen_vy->push_back(vty);Gen_vz->push_back(vtz);Gen_vt->push_back(vtt)

#define RecordSimTrack(thisTrk) ParToSim.push_back(thisTrk);			\
  IsParInHep->push_back(false);						\
  GenSimMomentum(thisTrk->momentum().pt(),thisTrk->momentum().eta(),thisTrk->momentum().phi(),thisTrk->type());	\
  if (!thisTrk->noVertex()) {						\
    SimVertex thisVtx=SVC[thisTrk->vertIndex()];			\
    GenSimVertex(thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9); \
  }									\
  else {GenSimVertex(0,0,0,0);}

#define RecordHepMC(GenParticle) ParToHep.push_back(GenParticle);	\
  IsParInHep->push_back(true);						\
  GenSimMomentum(sqrt(GenParticle->momentum().px()*GenParticle->momentum().px()+GenParticle->momentum().py()*GenParticle->momentum().py()),GenParticle->momentum().eta(),GenParticle->momentum().phi(),GenParticle->pdg_id()); \
  HepMC::GenVertex *thisVtx=GenParticle->production_vertex();		\
  if (thisVtx) {GenSimVertex(thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);} \
  else {GenSimVertex(0,0,0,0);}

// ------------ method called on each new Event  ------------
bool
CSCPriEff::filter(edm::Event& event, const edm::EventSetup& iSetup) {
  Summarization.Total_Events++;
  
  //Determine if keeps the event
  pt->clear();
  //Set Muons Handle
  Handle<reco::MuonCollection> Muon;
  event.getByLabel("muons", Muon);
  if (!Muon.isValid()) return false;
  reco::MuonCollection const & muons = *Muon;
  
  for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
    if (iter->isTrackerMuon()||iter->isGlobalMuon()) pt->push_back(iter->pt());
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
  HepMC::GenEvent * HepGenEvent = new HepMC::GenEvent(*(HepMCH->GetEvent()));
  HepMC::GenEvent::particle_iterator GenParticle_iter = HepGenEvent->particles_begin();
  for (;GenParticle_iter != HepGenEvent->particles_end();GenParticle_iter++) { 
    HepMC::GenVertex *thisVtx=(*GenParticle_iter)->production_vertex();
    if (thisVtx) {
      HepMCFourVec[0]=thisVtx->position().x()/10.;
      HepMCFourVec[1]=thisVtx->position().y()/10.;
      HepMCFourVec[2]=thisVtx->position().z()/10.;
      HepMCFourVec[3]=thisVtx->position().t()/299.792458;
      break;
    }
  }
  if (GenParticle_iter == HepGenEvent->particles_end()) {ReportError(7,"No HepMC(Core) Vertex Information")}

  //Reco Tracks
  Handle<reco::TrackCollection> trackCollectionH;
  event.getByLabel(tracksTag,trackCollectionH);
  reco::TrackCollection tC = *trackCollectionH.product();
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(tracksTag,trackCollectionHV);

  //Simulated Vertices: the vertexId() is just the position
  Handle<SimVertexContainer> SVCollectionH;
  event.getByLabel("g4SimHits", SVCollectionH);
  SVC = *SVCollectionH.product();

  //Simulated Tracks: the trackId() is not the position
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

  //Tracking Particles (collection of SimTracks and Hits)
  Handle<TrackingParticleCollection> TPCollectionH ;
  event.getByLabel("mergedtruthNoSimHits",TPCollectionH);

  //RecoToSim Muons Association
  ESHandle<TrackAssociatorBase> AssociatorByHits;
  iSetup.get<TrackAssociatorRecord>().get("muonAssociatorByHits_NoSimHits", AssociatorByHits);
  //const MuonAssociatorByHits * AssociatorByHits = dynamic_cast<const MuonAssociatorByHits *>(associatorBase.product());
  RecoToSimCollection RecoToSimByHits = AssociatorByHits->associateRecoToSim(trackCollectionHV,TPCollectionH,&event,&iSetup);
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
    if (iter->isTrackerMuon()||iter->isGlobalMuon()) {
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
	SelectorPassed[whichcut]->push_back(muon::isGoodMuon(*iter,muon::selectionTypeFromString(Cuts[whichcut])));
      MySelector->push_back(muon::isGoodMuon(*iter,muon::TMLastStation,1,3,3,3,3,maxChamberDist,maxChamberDistPull,Muon::SegmentAndTrackArbitration,false,true));
      //Station and Segment Matches
      StationMask->push_back(iter->stationMask(Muon::SegmentAndTrackArbitration));
      RequiredStationMask->push_back(muon::RequiredStationMask(*iter,maxChamberDist,maxChamberDistPull,Muon::SegmentAndTrackArbitration));
      for(int stationIdx = 0; stationIdx <4; ++stationIdx) {
	TrackDist[stationIdx]->push_back(iter->trackDist(stationIdx+1,2,Muon::SegmentAndTrackArbitration));//999999 means that there is no track
	TrackDistErr[stationIdx]->push_back(iter->trackDistErr(stationIdx+1,2,Muon::SegmentAndTrackArbitration));//999999 means that there is no track
	DXTrackToSegment[stationIdx]->push_back(999999);
	DYTrackToSegment[stationIdx]->push_back(999999);
	DXErrTrackToSegment[stationIdx]->push_back(999999);
	DYErrTrackToSegment[stationIdx]->push_back(999999);
	DRTrackToSegment[stationIdx]->push_back(999999);
	DRErrTrackToSegment[stationIdx]->push_back(999999);
	IsSegmentBelongsToTrackByDR[stationIdx]->push_back(false);
	IsSegmentBelongsToTrackByCleaning[stationIdx]->push_back(false);
      }
      for( std::vector<MuonChamberMatch>::const_iterator chamberMatch = iter->matches().begin();chamberMatch != iter->matches().end(); chamberMatch++ ) {
	if (chamberMatch->segmentMatches.empty()) continue;
	if (chamberMatch->detector()!=MuonSubdetId::CSC) continue;
	Byte_t stationIdx=chamberMatch->station()-1;
	for( std::vector<MuonSegmentMatch>::const_iterator segmentMatch = chamberMatch->segmentMatches.begin();segmentMatch != chamberMatch->segmentMatches.end(); segmentMatch++ ) {
	  if (!segmentMatch->isMask(MuonSegmentMatch::BestInChamberByDR)) continue;
	  IsSegmentBelongsToTrackByDR[stationIdx]->back()=segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByDR);
	  IsSegmentBelongsToTrackByCleaning[stationIdx]->back()=segmentMatch->isMask(MuonSegmentMatch::BelongsToTrackByCleaning);
	  DXTrackToSegment[stationIdx]->back()=abs(segmentMatch->x-chamberMatch->x);
	  DYTrackToSegment[stationIdx]->back()=abs(segmentMatch->y-chamberMatch->y);
	  DXErrTrackToSegment[stationIdx]->back()=sqrt(segmentMatch->xErr*segmentMatch->xErr+chamberMatch->xErr*chamberMatch->xErr);
	  DYErrTrackToSegment[stationIdx]->back()=sqrt(segmentMatch->yErr*segmentMatch->yErr+chamberMatch->yErr*chamberMatch->yErr);
	  DRTrackToSegment[stationIdx]->back()=sqrt(DXTrackToSegment[stationIdx]->back()*DXTrackToSegment[stationIdx]->back()+DYTrackToSegment[stationIdx]->back()*DYTrackToSegment[stationIdx]->back());
	  DRErrTrackToSegment[stationIdx]->back()=sqrt(DXTrackToSegment[stationIdx]->back()*DXTrackToSegment[stationIdx]->back()*DXErrTrackToSegment[stationIdx]->back()*DXErrTrackToSegment[stationIdx]->back()+DYTrackToSegment[stationIdx]->back()*DYTrackToSegment[stationIdx]->back()*DYErrTrackToSegment[stationIdx]->back()*DYErrTrackToSegment[stationIdx]->back())/DRTrackToSegment[stationIdx]->back();
	  break;
	}
      }
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
      TrkRecoChi2->push_back(iter->combinedQuality().trkRelChi2);
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
      reco::TrackRef innertrack = iter->innerTrack();
      dEdx->push_back(dEdxTrack[innertrack].dEdx());
      dEdxError->push_back(dEdxTrack[innertrack].dEdxError());
      dEdx_numberOfSaturatedMeasurements->push_back(dEdxTrack[innertrack].numberOfSaturatedMeasurements());
      dEdx_numberOfMeasurements->push_back(dEdxTrack[innertrack].numberOfMeasurements());
      theSameWithMuon->push_back(-1);
      
      reco::TrackRef glbtrack = iter->globalTrack();
      if (glbtrack.isNull()) cerr<<"InValidTrack"<<endl;
      unsigned int glbtrack_pos;  DChains->push_back(-2);
      vector <Byte_t> type;
      for(glbtrack_pos=0; glbtrack_pos<tC.size(); glbtrack_pos++)  //recotrk loop begin
	if (glbtrack == reco::TrackRef(trackCollectionH,glbtrack_pos) ) {
	  RefToBase<Track> trk(trackCollectionHV, glbtrack_pos);
	  nValidTrackerHits->push_back(trk->hitPattern().numberOfValidTrackerHits());
	  nValidPixelHits->push_back(trk->hitPattern().numberOfValidPixelHits());
	  nLostPixelHits->push_back(trk->hitPattern().numberOfLostPixelHits());
	  nLostTrackerHits->push_back(trk->hitPattern().numberOfLostTrackerHits());
	  nBadMuonHits->push_back(trk->hitPattern().numberOfBadMuonHits());
	  
	  if(RecoToSimByHits.find(trk) != RecoToSimByHits.end()) {//SimTrk begin
	    pair<TrackingParticleRef, double>  BestMatch=RecoToSimByHits[trk].front();
	    TrackingParticleRef tpr = BestMatch.first;
	    //Simulated Tracks
	    TrkParticles_pt->push_back(tpr->pt());
	    TrkParticles_eta->push_back(tpr->eta());
	    TrkParticles_phi->push_back(tpr->phi());
	    TrkParticles_charge->push_back(tpr->charge());
	    TrkParticles_pdgId->push_back(tpr->pdgId());
	    MatchQuality->push_back(BestMatch.second);
	    //cout<<"pt:"<<track->pt()<<"vs"<<tpr->pt()<<endl<<"eta:"<<track->eta()<<"vs"<<tpr->eta()<<endl<<"phi:"<<track->phi()<<"vs"<<tpr->phi()<<endl<<"chi2:"<<-tp.begin()->second<<endl<<"-------------Next------------------"<<endl;//check it is doing correct things
	    
	    //Get the decay chain of this track
	    MaskOut.clear();
	    for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )  {//g4Track loop begin
	      SimTrack *thisTrk;
	      trackIdLink(g4Track_iter->trackId())
	      if (find(MaskOut.begin(),MaskOut.end(),thisTrk)==MaskOut.end()) {
		DChain.clear();SimChains.clear();
		SimTrackDaughtersTree(thisTrk);
		SimVertex thisVtx;
		do {
		  if (!thisTrk->noVertex()) {
		    thisVtx=SVC[thisTrk->vertIndex()];
		    if (!thisVtx.noParent()) {
		      trackIdLink(thisVtx.parentIndex())
		      //add parent particle to each Chain
		      vector<SimTrack *>::iterator ParToSim_iter=find(ParToSim.begin(),ParToSim.end(),thisTrk);
		      if (ParToSim_iter==ParToSim.end()) {
			for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
			  SimChains_iter->insert(SimChains_iter->begin(),IsParInHep->size());
			RecordSimTrack(thisTrk)
		      }
		      else { Int_t pos=FindSimTrackRecordingPosition(ParToSim_iter-ParToSim.begin());
			for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
			  SimChains_iter->insert(SimChains_iter->begin(),pos);
		      }
		    }
		    else break;
		  }
		  else break;
		}
		while(true);
		//HepMC Particles
		HepMCChains.clear();
		if (!thisTrk->noGenpart()) {
		  HepMC::GenEvent::particle_iterator genPar = HepGenEvent->particles_begin();
		  for (int count=1; count<thisTrk->genpartIndex()&&genPar != HepGenEvent->particles_end(); count++ )
		    genPar++;
		  if (genPar != HepGenEvent->particles_end()) {
		    HepMCParentTree(*genPar);
		  }
		  else {ReportError(6,"genpartIndex() Error or HepMC is empty")}
		}
		//merge the HepMC and SimTrack Decay Chains
		for (vector< vector<Int_t> >::iterator SimChains_iter = SimChains.begin(); SimChains_iter !=  SimChains.end(); ++SimChains_iter )
		  for (vector< vector<Int_t> >::iterator HepMCChains_iter = HepMCChains.begin(); HepMCChains_iter !=  HepMCChains.end(); ++HepMCChains_iter ) {
		    vector<Int_t> thisChain(HepMCChains_iter->rbegin(),HepMCChains_iter->rend());
		    thisChain.insert(thisChain.end(),SimChains_iter->begin(),SimChains_iter->end());
		    //see if thisChain is the same with previous muons
		    int Muref=-1;
		    for (vector<int>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end()&&Muref<(int) eta->size()-1; DChain_iter++ ) {
		      if (*DChain_iter==-1) {
			DChain_iter++;
			vector<int>::iterator DChain_begin=DChain_iter;
			for (; DChain_iter != DChains->end()&&*DChain_iter!=-2&&*DChain_iter!=-1; DChain_iter++ ) ;
			if (DChain_begin!=DChain_iter) 
			  if (IstheSameDChain(thisChain,* new vector<int> (DChain_begin,DChain_iter))) {
			    theSameWithMuon->back()=Muref;
			    break;
			  }
		      }
		      if (*DChain_iter==-2) Muref++;
		    }
		    DChains->push_back(-1);
		    DChains->insert(DChains->end(),thisChain.begin(),thisChain.end());
		    Byte_t newtype=classification(thisChain);
		    if (find(type.begin(),type.end(),newtype)==type.end()) type.push_back(newtype);
		  }
	      }
	    }//g4Track loop end
	    Long64_t TypeRecord=0;
	    /*for (vector <Byte_t>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	      if (*type_iter/10==5) 
	      for (vector <Byte_t>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	      if (*type_iter/10==3||*type_iter/10==4) *type_iter=0;*/
	    for (vector <Byte_t>::iterator type_iter=type.begin();type_iter!=type.end();type_iter++)
	      if (*type_iter) TypeRecord=TypeRecord*100+(*type_iter);
	    MuonType->push_back(TypeRecord);
	    break;
	  }//SimTrk End
	  //if no matched simulated track found
	  NumMisMatch++;
	  break;
	}//recotrk loop end
      if (TrkParticles_eta->size()<eta->size()) {
	if (glbtrack_pos==tC.size()) {ReportError(4,"TrkRef Error")}
	if (TrkParticles_eta->size()!=eta->size()-1) {ReportError(5,"MC Match Coding Error")}
	nValidTrackerHits->push_back(0);
	nValidPixelHits->push_back(0);
	nLostPixelHits->push_back(0);
	nLostTrackerHits->push_back(0);
	nBadMuonHits->push_back(0);
	TrkParticles_pt->push_back(0);
	TrkParticles_eta->push_back(0);
	TrkParticles_phi->push_back(0);
	TrkParticles_pdgId->push_back(0);
	TrkParticles_charge->push_back(-100.0);
	MatchQuality->push_back(-100.0);
	MuonType->push_back(0);
      }
    }//muon loop end
#ifdef LocalRun
  bool Print=false;
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
    printf("\n\n");
    unsigned int Muref=0;
    for (vector<Int_t>::iterator DChain_iter = DChains->begin(); DChain_iter != DChains->end(); DChain_iter++ )
      {
	if (*DChain_iter==-1) printf("\n");
	else if (*DChain_iter==-2) {
	  printf("\n Mu (%d): (%3.2fcm,%3.2fcm,%3.2fcm) --- (%3.2fGeV/c,%3.2f,%3.2f) ---> Type",Muref,Vertex_x->at(Muref),Vertex_y->at(Muref),Vertex_z->at(Muref),pt->at(Muref)*cosh(eta->at(Muref)),eta->at(Muref),phi->at(Muref));
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
bool CSCPriEff::IstheSameDChain(vector<int> &ThisChain,vector<int> &AnotherChain)
{
  bool ChainIncluded=false;
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

//classification tracker muons
//21-prompt Mu from LightMeson
//22-prompt Mu from HeavyMeson
//23-prompt Mu from LightBaryon
//24-prompt Mu from HeavyBaryon
//25-prompt Mu from WZ
//26-Mu from WZ but not prompt (should not happen)
//27-prompt Mu from Others (should not happen)
//30-punch through: it is a hadron and doesn't decay in the detector region
//40-punch through and decay in flight: it is a hadron and decays outside of HCAL (including HO)
//51-Decay in flight from LightMeson
//52-Decay in flight from HeavyMeson
//53-Decay in flight from LightBaryon
//54-Decay in flight from HeavyBaryon
//60-others
//0-no matched simulated tracks
//number combination: duplicate decay chain types(e.g. 23 means prompt Mu from jets and a punch through )
Byte_t CSCPriEff::classification(vector<Int_t> &Chain)
{
  int MuPos=-1;
  for (vector<Int_t>::reverse_iterator iter = Chain.rbegin(); iter != Chain.rend(); iter++) {
    int pid=(*Gen_pdgId)[*iter];
    ParticleType ParType=ParticleCata(pid);
    if (abs(pid)==13) MuPos=*iter;
    if (MuPos>=0&&(ParType==W||ParType==Z)) {
      if ((*IsParInHep)[MuPos]) return 25;
      else return 26;
    }
    if (MuPos<0&&ParType>=LightMeson&&ParType<=BottomBaryon) return 30;
    if (MuPos>=0&&ParType>=LightMeson&&ParType<=BottomBaryon) {
      float deltaRPhi2=((*Gen_vx)[MuPos]-HepMCFourVec[0])*((*Gen_vx)[MuPos]-HepMCFourVec[0])+((*Gen_vy)[MuPos]-HepMCFourVec[1])*((*Gen_vy)[MuPos]-HepMCFourVec[1]);
      float deltaZ=abs((*Gen_vz)[MuPos]-HepMCFourVec[2]);
      float RPhi2=((*Gen_vx)[MuPos])*((*Gen_vx)[MuPos])+((*Gen_vy)[MuPos])*((*Gen_vy)[MuPos]);
      float Z=abs((*Gen_vz)[MuPos]);
      if ((*IsParInHep)[MuPos]||(deltaRPhi2<9&&deltaZ<30)) {//vertex region
	if (ParType==LightMeson) return 21;
	if (ParType>LightMeson&&ParType<LightBaryon) return 22;
	if (ParType==LightBaryon) return 23;
	if (ParType>LightBaryon&&ParType<=BottomBaryon) return 24;
      }
      if ((RPhi2<161604.&&Z<568.)||(RPhi2>82024.96&&RPhi2<161604.&&Z<666.)) {//inside HCAL
	if (ParType==LightMeson) return 51;
	if (ParType>LightMeson&&ParType<LightBaryon) return 52;
	if (ParType==LightBaryon) return 53;
	if (ParType>LightBaryon&&ParType<=BottomBaryon) return 54;
      }
      return 40;
    }
  }
  if (MuPos>=0) return 26;
  return 60;
}


void CSCPriEff::SimTrackDaughtersTree(SimTrack * thisTrk)
{
  //To avoid duplicate particle saving
  vector<SimTrack *>::iterator ParToSim_iter=find(ParToSim.begin(),ParToSim.end(),thisTrk);
  if (ParToSim_iter==ParToSim.end())
    {
      DChain.push_back(IsParInHep->size());
      RecordSimTrack(thisTrk)
    }
  else DChain.push_back(FindSimTrackRecordingPosition(ParToSim_iter-ParToSim.begin()));
  MaskOut.push_back(thisTrk);
  if (Daughters[thisTrk].size()>0) 
    for (vector<SimTrack *>::iterator Daughter=Daughters[thisTrk].begin();Daughter!=Daughters[thisTrk].end();Daughter++)
      SimTrackDaughtersTree(*Daughter);
  else SimChains.push_back(DChain);
  DChain.pop_back();
}


void CSCPriEff::HepMCParentTree(HepMC::GenParticle *genPar) {
  HepMC::GenVertex *thisVtx = genPar->production_vertex();
  bool ChainEnd=true;
  if (thisVtx) {
      for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	if ((*pgenD)->pdg_id()!=92)  {//Pythia special code for string, we only care about the particles after hadronization
	  ChainEnd=false;
	  vector<HepMC::GenParticle *>::iterator ParToHep_iter=find(ParToHep.begin(),ParToHep.end(),*pgenD);
	  if (ParToHep_iter==ParToHep.end())
	    {
	      DChain.push_back(IsParInHep->size());
	      RecordHepMC((*pgenD))
	    }
	  else DChain.push_back(FindHepMCRecordingPosition(ParToHep_iter-ParToHep.begin()));
	  HepMCParentTree(*pgenD);
	  DChain.pop_back();
	}
  }
  if (ChainEnd) HepMCChains.push_back(DChain);
}

void CSCPriEff::ClearVecs() {
//clear all the other vectors except pt
  eta->clear();  phi->clear();  chargeMinus->clear();
  isGlobalMu->clear();  isTrackerMu->clear();
  Vertex_x->clear();  Vertex_y->clear();  Vertex_z->clear();
  isoR03sumPt->clear();  isoR03emEt->clear();  isoR03hadEt->clear();
  isoR03hoEt->clear();  isoR03nJets->clear();  isoR03nTracks->clear();
  isoR05sumPt->clear();  isoR05emEt->clear();  isoR05hadEt->clear();
  isoR05hoEt->clear();  isoR05nJets->clear();  isoR05nTracks->clear();
  isoemVetoEt->clear();  isohadVetoEt->clear();  isohoVetoEt->clear();
  TrkKink->clear();  GlbKink->clear();  TrkRecoChi2->clear();
  CaloE_emMax->clear();  CaloE_emS9->clear();  CaloE_emS25->clear();
  CaloE_hadMax->clear();  CaloE_hadS9->clear();
  Calo_emPos_R->clear();  Calo_emPos_eta->clear();  Calo_emPos_phi->clear();
  Calo_hadPos_R->clear();  Calo_hadPos_eta->clear();  Calo_hadPos_phi->clear();

  dEdx->clear();  dEdxError->clear();
  dEdx_numberOfSaturatedMeasurements->clear();  dEdx_numberOfMeasurements->clear();

  nValidTrackerHits->clear(); nValidPixelHits->clear();  
  nLostTrackerHits->clear(); nLostPixelHits->clear(); 
  nBadMuonHits->clear();
  
  StationMask->clear();RequiredStationMask->clear();
  for(unsigned int stationIdx = 0; stationIdx <4; ++stationIdx) {
    TrackDist[stationIdx]->clear();
    TrackDistErr[stationIdx]->clear();
    DXTrackToSegment[stationIdx]->clear();
    DYTrackToSegment[stationIdx]->clear();
    DXErrTrackToSegment[stationIdx]->clear();
    DYErrTrackToSegment[stationIdx]->clear();
    DRTrackToSegment[stationIdx]->clear();
    DRErrTrackToSegment[stationIdx]->clear();
    IsSegmentBelongsToTrackByDR[stationIdx]->clear();
    IsSegmentBelongsToTrackByCleaning[stationIdx]->clear();
  }

  vx->clear();   vxError->clear();
  vy->clear();   vyError->clear();
  vz->clear();   vzError->clear();

  Gen_pt->clear();  Gen_eta->clear();  Gen_phi->clear(); Gen_pdgId->clear();
  Gen_vx->clear();  Gen_vy->clear();  Gen_vz->clear();  Gen_vt->clear();
  
  MySelector->clear();
  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++)
    SelectorPassed[whichcut]->clear(); 

  TrkParticles_pt->clear();  TrkParticles_eta->clear(); TrkParticles_phi->clear();
  TrkParticles_pdgId->clear(); TrkParticles_charge->clear(); MatchQuality->clear();

  ParToSim.clear(); ParToHep.clear(); IsParInHep->clear();
  DChains->clear(); theSameWithMuon->clear(); MuonType->clear();
}

#define MakeVecBranch(Name,Var,Type) Var=new vector<Type>();Tracks_Tree->Branch(Name,&Var)

CSCPriEff::CSCPriEff(const edm::ParameterSet& pset) {
//---- Get the input parameters
  FileName = pset.getUntrackedParameter<string>("FileName","MuonTracks.root");
  maxChamberDist = pset.getUntrackedParameter<double>("maxChamberDist",-3.);
  maxChamberDistPull = pset.getUntrackedParameter<double>("maxChamberDistPull",-3.);
  Cuts = pset.getUntrackedParameter< vector<string> >("StandardMuonCuts",Cuts);//default is empty
  PrimaryVerticesTag = pset.getUntrackedParameter<string>("PrimaryVertices","offlinePrimaryVerticesWithBS");//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
  const InputTag tracksTag_default("generalTracks");
  tracksTag = pset.getUntrackedParameter<InputTag>("tracksTag",tracksTag_default);
  dEdxTag = pset.getUntrackedParameter<string>("dEdxTag","dedxHarmonic2");//Other options already available in RECO files are dedxMedian and dedxTruncated40. 
  HepMCTag = pset.getUntrackedParameter<string>("HepMCTag","generator");

//Build a root file and a TTree
  file=new TFile(FileName.c_str(),"RECREATE");
  Tracks_Tree = new TTree ("MuTrkCand","TrksCandidates") ;
  Summarization_Tree = new TTree ("Summary","Summary") ;
  ErrorMsg_Tree = new TTree ("Errors","Errors") ;
  Tracks_Tree->SetCircular(500000);//the max events in a single root file
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
  MakeVecBranch("TrkKink",TrkKink,Float_t); MakeVecBranch("GlbKink",GlbKink,Float_t); MakeVecBranch("TrkRecoChi2",TrkRecoChi2,Float_t);
  MakeVecBranch("CaloE_emMax",CaloE_emMax,Float_t);  MakeVecBranch("CaloE_emS9",CaloE_emS9,Float_t);  MakeVecBranch("CaloE_emS25",CaloE_emS25,Float_t);
  MakeVecBranch("CaloE_hadMax",CaloE_hadMax,Float_t);  MakeVecBranch("CaloE_hadS9",CaloE_hadS9,Float_t);
  MakeVecBranch("Calo_emPos_R",Calo_emPos_R,Float_t);  MakeVecBranch("Calo_emPos_eta",Calo_emPos_eta,Float_t);  MakeVecBranch("Calo_emPos_phi",Calo_emPos_phi,Float_t);
  MakeVecBranch("Calo_hadPos_R",Calo_hadPos_R,Float_t);  MakeVecBranch("Calo_hadPos_eta",Calo_hadPos_eta,Float_t);  MakeVecBranch("Calo_hadPos_phi",Calo_hadPos_phi,Float_t);
  //InnerTrack
  MakeVecBranch("dEdx",dEdx,Float_t);
  MakeVecBranch("dEdxError",dEdxError,Float_t);
  MakeVecBranch("dEdx_numberOfSaturatedMeasurements",dEdx_numberOfSaturatedMeasurements,Int_t);
  MakeVecBranch("dEdx_numberOfMeasurements",dEdx_numberOfMeasurements,Int_t);
  MakeVecBranch("nValidTrackerHits",nValidTrackerHits,UInt_t);
  MakeVecBranch("nValidPixelHits",nValidPixelHits,UInt_t);
  MakeVecBranch("nLostTrackerHits",nLostTrackerHits,UInt_t);
  MakeVecBranch("nLostPixelHits",nLostPixelHits,UInt_t);
  MakeVecBranch("nBadMuonHits",nBadMuonHits,UInt_t);
  //Muon System
  char TrackDistBranchName[14]="TrackDistCSC1",TrackDistErrBranchName[17]="TrackDistErrCSC1",
    DXTrackToSegmentName[21]="DXTrackToSegmentCSC1",DYTrackToSegmentName[21]="DYTrackToSegmentCSC1",
    DXErrTrackToSegmentName[24]="DXErrTrackToSegmentCSC1",DYErrTrackToSegmentName[24]="DYErrTrackToSegmentCSC1",
    DRTrackToSegmentName[21]="DRTrackToSegmentCSC1",DRErrTrackToSegmentName[24]="DRErrTrackToSegmentCSC1",
    IsSegmentBelongsToTrackByDRName[32]="IsSegmentBelongsToTrackByDRCSC1",IsSegmentBelongsToTrackByCleaningName[38]="IsSegmentBelongsToTrackByCleaningCSC1";
  for(unsigned int stationIdx = 0; stationIdx <4; ++stationIdx) {
    MakeVecBranch(TrackDistBranchName,TrackDist[stationIdx],Float_t);
    MakeVecBranch(TrackDistErrBranchName,TrackDistErr[stationIdx],Float_t);
    MakeVecBranch(DXTrackToSegmentName,DXTrackToSegment[stationIdx],Float_t);
    MakeVecBranch(DYTrackToSegmentName,DYTrackToSegment[stationIdx],Float_t);
    MakeVecBranch(DXErrTrackToSegmentName,DXErrTrackToSegment[stationIdx],Float_t);
    MakeVecBranch(DYErrTrackToSegmentName,DYErrTrackToSegment[stationIdx],Float_t);
    MakeVecBranch(DRTrackToSegmentName,DRTrackToSegment[stationIdx],Float_t);
    MakeVecBranch(DRErrTrackToSegmentName,DRErrTrackToSegment[stationIdx],Float_t);
    MakeVecBranch(IsSegmentBelongsToTrackByDRName,IsSegmentBelongsToTrackByDR[stationIdx],Bool_t);
    MakeVecBranch(IsSegmentBelongsToTrackByCleaningName,IsSegmentBelongsToTrackByCleaning[stationIdx],Bool_t);
    TrackDistBranchName[12]++; TrackDistErrBranchName[15]++;
    DXTrackToSegmentName[19]++;DYTrackToSegmentName[19]++;
    DXErrTrackToSegmentName[22]++;DYErrTrackToSegmentName[22]++;
    DRTrackToSegmentName[19]++;DRErrTrackToSegmentName[22]++;
    IsSegmentBelongsToTrackByDRName[30]++;IsSegmentBelongsToTrackByCleaningName[36]++;
  }
  MakeVecBranch("StationMask",StationMask,UInt_t);MakeVecBranch("RequiredStationMask",RequiredStationMask,UInt_t); 
  //PV
  MakeVecBranch("PVx",vx,Float_t);  MakeVecBranch("PVy",vy,Float_t);  MakeVecBranch("PVz",vz,Float_t);
  MakeVecBranch("PVxError",vxError,Float_t);  MakeVecBranch("PVyError",vyError,Float_t);  MakeVecBranch("PVzError",vzError,Float_t);
  //Generation Level Muons
  MakeVecBranch("Gen_pt",Gen_pt,Float_t);  MakeVecBranch("Gen_eta",Gen_eta,Float_t);  MakeVecBranch("Gen_phi",Gen_phi,Float_t);   MakeVecBranch("Gen_pdgId",Gen_pdgId,Int_t);
  MakeVecBranch("Gen_vx",Gen_vx,Float_t);  MakeVecBranch("Gen_vy",Gen_vy,Float_t);  MakeVecBranch("Gen_vz",Gen_vz,Float_t);  MakeVecBranch("Gen_vt",Gen_vt,Float_t);
  MakeVecBranch("IsParInHep",IsParInHep,Bool_t);
  Tracks_Tree->Branch("HepMCVertex",HepMCFourVec,"HepMCVertex[4]/F");

  //Simulated Tracks
  MakeVecBranch("TrkParticles_pt",TrkParticles_pt,Float_t);  MakeVecBranch("TrkParticles_eta",TrkParticles_eta,Float_t);  MakeVecBranch("TrkParticles_phi",TrkParticles_phi,Float_t); MakeVecBranch("TrkParticles_charge",TrkParticles_charge,Int_t);
  MakeVecBranch("TrkParticles_pdgId",TrkParticles_pdgId,Int_t);  MakeVecBranch("MatchQuality",MatchQuality,Double_t);
  MakeVecBranch("DChains",DChains,Int_t); MakeVecBranch("theSameWithMuon",theSameWithMuon,Int_t);
  MakeVecBranch("MuonType",MuonType,Long64_t);
  //Cuts
  num_Cuts=Cuts.size();
  for (unsigned int whichcut=0;whichcut<num_Cuts;whichcut++) 
    if (muon::selectionTypeFromString(Cuts[whichcut])!=(muon::SelectionType)-1) {
      MakeVecBranch(Cuts[whichcut].c_str(),SelectorPassed[whichcut],Bool_t); 
    }
    else {
      Cuts.erase(Cuts.begin()+whichcut);
      num_Cuts--;
    }
  MakeVecBranch("MySelector",MySelector,Bool_t);
  //Summary
  ErrorMsg_Tree->Branch("ErrorMsg",&Error.ErrorCode,"ErrorCode/b:RunNum/l:EventNum");
  Summarization_Tree->Branch("Summarization",&Summarization.Total_Events,"Total_Events/l:Total_TrackerMuons:Total_GlobalMuon:Total_GlobalnotTrackerMuon");
  Summarization.Total_Events=0;Summarization.Total_TrackerMuons=0;Summarization.Total_GlobalMuon=0;Summarization.Total_GlobalnotTrackerMuon=0;
  FirstEntry=true;
  NumMisMatch=0;
}

CSCPriEff::~CSCPriEff()
{
   printf("MisRate:%f%%\n",NumMisMatch/(float) Summarization.Total_TrackerMuons*100);
   //Summarizations
   Summarization_Tree->Fill();
   file->Write();
   file->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCPriEff);
