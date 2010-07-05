// -*- C++ -*-
//
// Package:    CSCPriEff
// Class:      CSCPriEff
// 
/**\class CSCPriEff CSCPriEff.cc UserCode/CSCPriEff/src/CSCPriEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhang Jinzhong
//         Created:  Mon Jun  7 22:19:50 CEST 2010
// $Id$
//
//

#include "../interface/CSCPriEff.h"

using namespace edm;
using namespace std;
using namespace reco;

#define ReportError(code,MSG) Error.ErrorCode=code;    \
  Error.Run=Info.RUN;			      \
  Error.Event=Info.EVENT;		      \
  ErrorMsg_Tree->Fill();		      \
  if (LocalRun) cerr<<MSG<<endl;

#define trackIdLink(TrkID)  unsigned int pos=0;				\
  vector<SimTrack>::const_iterator LinkToTrk;				\
  for (LinkToTrk = STC.begin(); LinkToTrk != STC.end(); ++LinkToTrk )	\
    if ((int) LinkToTrk->trackId()== (int) TrkID)				\
      {pos=LinkToTrk-STC.begin();break;}				\
  if (LinkToTrk==STC.end()) {ReportError(8,"parentIndex/trackId Error")}

#define GenSimMomentum(ppt,peta,pphi,ppdg,pcharge) Gen_pdgId->push_back(ppdg);Gen_pt->push_back(ppt);Gen_eta->push_back(peta);Gen_phi->push_back(pphi);Gen_charge->push_back(pcharge);
#define GenSimVertex(vtx,vty,vtz,vtt) Gen_vx->push_back(vtx);Gen_vy->push_back(vty);Gen_vz->push_back(vtz);Gen_vt->push_back(vtt);


// ------------ method called on each new Event  ------------
bool
CSCPriEff::filter(edm::Event& event, const edm::EventSetup& iSetup)
{
  Summarization.Total_Events++;

//Determine if keeps the event
  pt->clear();
  //Set Muons Handle
  Handle<reco::MuonCollection> Muon;
  event.getByLabel("muons", Muon);
  if (!Muon.isValid()) return false;
  reco::MuonCollection const & muons = *Muon;
  
  for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
    if (iter->isTrackerMuon()) pt->push_back(iter->pt());
  if (pt->size()==0) return false;

//general event information
  Info.RUN   = event.id ().run ();
  Info.EVENT = event.id ().event();
  Info.LS    = event.luminosityBlock ();
  Info.ORBIT = event.orbitNumber ();
  Info.BX = event.bunchCrossing ();

//clear all the other vectors except pt
  eta->clear();  phi->clear();
  whichGenPar->clear(); whichGenId->clear();
  chargeMinus->clear();  isGlobalMu->clear();
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
  
  vx->clear();   vxError->clear();
  vy->clear();   vyError->clear();
  vz->clear();   vzError->clear();

  Gen_pt->clear();  Gen_eta->clear();  Gen_phi->clear(); Gen_pdgId->clear();
  Gen_vx->clear();  Gen_vy->clear();  Gen_vz->clear();  Gen_vt->clear();
  Gen_charge->clear();

  TrkParticles_pt->clear();  TrkParticles_eta->clear(); TrkParticles_phi->clear();
  TrkParticles_pdgId->clear(); DChains->clear();

  //dEdx
  Handle< ValueMap<reco::DeDxData> > dEdxTrackHandle;
  event.getByLabel(dEdxTag, dEdxTrackHandle);
  const ValueMap<reco::DeDxData> dEdxTrack = *(dEdxTrackHandle.product());

  //HepMC Particles
  Handle<edm::HepMCProduct> HepMCH;
  event.getByLabel(HepMCTag, HepMCH);
  HepMC::GenEvent * HepGenEvent = new HepMC::GenEvent(*(HepMCH->GetEvent()));
  
  //Reco Tracks
  Handle<reco::TrackCollection> trackCollectionH;
  event.getByLabel(tracksTag,trackCollectionH);
  reco::TrackCollection tC = *trackCollectionH.product();
  Handle< View<Track> > trackCollectionHV;
  event.getByLabel(tracksTag,trackCollectionHV);

  //Tracking Particles (collection of SimTracks and Hits)
  Handle<TrackingParticleCollection> TPCollectionH ;
  event.getByType(TPCollectionH);

  //Simulated Tracks: the trackId() is not the position
  Handle<SimTrackContainer> STCollectionH;
  event.getByLabel("g4SimHits", STCollectionH);
  SimTrackContainer STC = *STCollectionH.product();

  //Simulated Vertices: the vertexId() is just the position
  Handle<SimVertexContainer> SVCollectionH;
  event.getByLabel("g4SimHits", SVCollectionH);
  SimVertexContainer SVC = *SVCollectionH.product();

  //RecoToSim Tracks Association
  ESHandle<TrackAssociatorBase> Associator;
  iSetup.get<TrackAssociatorRecord>().get(MatchAlgorithm, Associator);
  RecoToSimCollection RecoToSim = Associator->associateRecoToSim(trackCollectionHV,TPCollectionH,&event);
  particlenum=0;
  for ( reco::MuonCollection::const_iterator iter = muons.begin(); iter != muons.end() ; ++iter)
    //muon loop begin
    if (iter->isTrackerMuon())
      {
	//muon information (pt has been recorded in the beginning of this function)
	eta->push_back(iter->eta());
	phi->push_back(iter->phi());
	if (iter->charge()==-1) chargeMinus->push_back(true);
	else chargeMinus->push_back(false);
	isGlobalMu->push_back(iter->isGlobalMuon());
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
	reco::TrackRef track = iter->innerTrack();
	dEdx->push_back(dEdxTrack[track].dEdx());
	dEdxError->push_back(dEdxTrack[track].dEdxError());
	dEdx_numberOfSaturatedMeasurements->push_back(dEdxTrack[track].numberOfSaturatedMeasurements());
	dEdx_numberOfMeasurements->push_back(dEdxTrack[track].numberOfMeasurements());

	DChains->push_back(-2);
	unsigned int i;
	for(i=0; i<tC.size(); i++)
	  //inner recotrk loop begin
	  if (track == reco::TrackRef(trackCollectionH,i) )
	    {
	      RefToBase<Track> trk(trackCollectionHV, i);
	      vector<pair<TrackingParticleRef, double> > tp;
	      if(RecoToSim.find(trk) != RecoToSim.end())
		{
		  tp = RecoToSim[trk];
		  if (tp.size()!=0)
		    {
		      TrackingParticleRef tpr = tp.begin()->first;
      		      //Simulated Tracks
		      TrkParticles_pt->push_back(tpr->pt());
		      TrkParticles_eta->push_back(tpr->eta());
		      TrkParticles_phi->push_back(tpr->phi());
		      TrkParticles_pdgId->push_back(tpr->pdgId());
		      //save the decay chain from the final product to the beginning
		      vector<bool> MaskOut(tpr->g4Tracks().size(),false);
		      //if a SimTrack in g4Tracks() has a parent SimTrack in g4Tracks(), mask the daughter out (to prevent dupicate tracks saving)
		      for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )
			{
			  SimTrack thisTrk=*g4Track_iter;
			  do
			    {
			      if (!thisTrk.noVertex())
				{
				  SimVertex thisVtx=SVC[thisTrk.vertIndex()];
				  if (!thisVtx.noParent())
				    {
				      vector<SimTrack>::const_iterator FindTrack;
				      for (FindTrack = tpr->g4Track_begin(); FindTrack != tpr->g4Track_end(); ++FindTrack )
					{
					  trackIdLink(thisVtx.parentIndex())
					  thisTrk=STC[pos];
					  if (thisTrk.trackId()==FindTrack->trackId()) 
					    {
					      MaskOut[FindTrack-tpr->g4Track_begin()]==true;
					      break;
					    }
					}
				    }
				  else break;
				}
			      else break;
			    }
			  while(true);
			}
		      //find daughters of each SimTrack
		      vector< vector<int> > daughters(STC.size(),* new vector<int>);
		      for (vector<SimTrack>::const_iterator FindTrack = STC.begin(); FindTrack != STC.end(); ++FindTrack )
			if (!FindTrack->noVertex())
			  {
			    SimVertex thisVtx=SVC[FindTrack->vertIndex()];
			    if (!thisVtx.noParent())
			      {
				trackIdLink(thisVtx.parentIndex())
				daughters[pos].push_back(FindTrack-STC.begin());
			      }
			  }
		      //Get the decay chain of this track
		      printf("\n------*Muon (%d) related TrackingParticles*------",eta.size()-1);
		      for (vector<SimTrack>::const_iterator g4Track_iter = tpr->g4Track_begin(); g4Track_iter != tpr->g4Track_end(); ++g4Track_iter )
			if (!MaskOut[g4Track_iter-tpr->g4Track_begin()])
			  {
			    //cout<<"pt:"<<track->pt()<<"vs"<<tpr->pt()<<endl<<"eta:"<<track->eta()<<"vs"<<tpr->eta()<<endl<<"phi:"<<track->phi()<<"vs"<<tpr->phi()<<endl<<"chi2:"<<-tp.begin()->second<<endl<<"-------------Next------------------"<<endl;//check it is doing correct things
			    SimTrack thisTrk=*g4Track_iter;
			    trackIdLink(thisTrk.trackId())
			    printf("\n---!SimTracks!---\n<daughters>\n");
			    DChain.clear();
			    SimTrackDaughtersTree(SVC,STC,daughters,pos);
			    printf("<parents>\n");
			    SimVertex thisVtx;
			    do
			      {
				printf("%d ",thisTrk.type());
				if (!thisTrk.noVertex())
				  {
				    thisVtx=SVC[thisTrk.vertIndex()];
				    printf("(%5.5f,%5.5f,%5.5f,%fus) --> ",thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9);
				    if (!thisVtx.noParent()) 
				      {
					thisTrk=STC[thisVtx.parentIndex()];
					//add parent particle to each Chain
					if (DChains->size()>2)
					  {
					    for (vector<Int_t>::iterator DChains_iter = DChains->end(); (DChains_iter != DChains->begin()+1)&&((*DChains_iter)!=-2) ; --DChains_iter )
					      if ((*DChains_iter)>=0&&(*(DChains_iter-1))<0) DChains->insert(DChains_iter,particlenum);
					  }
					particlenum++;
					GenSimMomentum(thisTrk.momentum().pt(),thisTrk.momentum().eta(),thisTrk.momentum().phi(),thisTrk.type(),(Byte_t) thisTrk.charge())
					if (!thisTrk.noVertex())
					  {
					    thisVtx=SVC[thisTrk.vertIndex()];
					    GenSimVertex(thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9)
					  }
					else {GenSimVertex(0,0,0,0)}
				      }
				    else break;
				  }
				else break;
			      }
			    while(true);
			    if (!thisTrk.noGenpart())
			      {
				HepMC::GenEvent::particle_iterator genPar = HepGenEvent->particles_begin();
				for (int count=1; count<thisTrk.genpartIndex(); count++ )
				  genPar++;
				printf("\n---!HepMC!---\n<Parents>\n");
				HepMCParentTree(*genPar);
				//printf("\n<Daughters>\n");
				//HepMCDaughtersTree(*genPar);//not necessary for simulated tracks
			      }
			    printf("\n--------------------------------------------");
			  }
		    }
		}
	      else
		{
		  ReportError(6,"No matched simulated track found")
		  TrkParticles_pt->push_back(0);
		  TrkParticles_eta->push_back(0);
		  TrkParticles_phi->push_back(0);
		  TrkParticles_pdgId->push_back(0);
		}
	      break;
	    }//inner recotrk loop end
	if (i==tC.size()) {ReportError(7,"TrkRef Error")}
      }//muon loop end
  printf("\n");
  for (unsigned int parref = 0; parref<Gen_pdgId->size(); parref++)
    printf("%d - %d; ",parref,Gen_pdgId->at(parref));
  printf("\n");
  for (vector<Int_t>::const_iterator DChains_iter = DChains->begin(); DChains_iter != DChains->end(); ++DChains_iter )
    {
      printf("%d --> ", *DChains_iter);
    }

//Primary Vertex
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(PrimaryVerticesTag.c_str(),recVtxs);
  if (recVtxs.isValid())
    for(reco::VertexCollection::const_iterator v=recVtxs->begin(); v!=recVtxs->end(); ++v)
      {
	vx->push_back(v->x());	   vxError->push_back(v->xError());
	vy->push_back(v->y());	   vyError->push_back(v->yError());
	vz->push_back(v->z());	   vzError->push_back(v->zError());
      }
  else if (FirstEntry) {ReportError(2,PrimaryVerticesTag.c_str()<<" information is not valid.")}

  FirstEntry=false;
  Tracks_Tree->Fill();
  return true;
}

CSCPriEff::CSCPriEff(const edm::ParameterSet& pset)
{
//---- Get the input parameters
  FileName = pset.getUntrackedParameter<string>("FileName","MuonTracks.root");
  minP = pset.getUntrackedParameter<double>("minP",20.);//
  maxP = pset.getUntrackedParameter<double>("maxP",100.);//
  LocalRun = pset.getUntrackedParameter<bool>("LocalRun",false);//
  MaxDR = pset.getUntrackedParameter<double>("MaxDR",0.15);//
  MaxRelpT = pset.getUntrackedParameter<double>("MaxRelpT",0.2);//
  MaxDR2 = MaxDR*MaxDR;
  MatchAlgorithm = pset.getUntrackedParameter<string>("RecoToSimAlgorithm","TrackAssociatorByChi2");
  minTrackHits = pset.getUntrackedParameter<uint>("minTrackHits",10);//
  PrimaryVerticesTag = pset.getUntrackedParameter<string>("PrimaryVertices","offlinePrimaryVerticesWithBS");//"offlinePrimaryVerticesWithBS": Primary vertex reconstructed using the tracks taken from the generalTracks collection, and imposing the offline beam spot as a constraint in the fit of the vertex position. Another possible tag is "offlinePrimaryVertices", which is Primary vertex reconstructed using the tracks taken from the generalTracks collection
  tracksTag = pset.getUntrackedParameter<string>("tracksTag","generalTracks");
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
  pt =  new vector<Float_t>();  eta = new vector<Float_t>();  phi = new vector<Float_t>();
  chargeMinus = new vector<Bool_t>();  isGlobalMu = new vector<Bool_t>();
  whichGenPar = new vector<Int_t>(); whichGenId = new vector<Int_t>();
  Vertex_x = new vector<Float_t>();  Vertex_y = new vector<Float_t>();  Vertex_z = new vector<Float_t>();
  isoR03sumPt = new vector<Float_t>();  isoR03emEt = new vector<Float_t>();  isoR03hadEt = new vector<Float_t>();
  isoR03hoEt = new vector<Float_t>();  isoR03nJets = new vector<Float_t>();  isoR03nTracks = new vector<Float_t>();
  isoR05sumPt = new vector<Float_t>();  isoR05emEt = new vector<Float_t>();  isoR05hadEt = new vector<Float_t>();
  isoR05hoEt = new vector<Float_t>();  isoR05nJets = new vector<Float_t>();  isoR05nTracks = new vector<Float_t>();
  isoemVetoEt = new vector<Float_t>();  isohadVetoEt = new vector<Float_t>();  isohoVetoEt = new vector<Float_t>();
  TrkKink =  new vector<Float_t>(); GlbKink =  new vector<Float_t>(); TrkRecoChi2 = new vector<Float_t>();
  CaloE_emMax = new vector<Float_t>();  CaloE_emS9 = new vector<Float_t>();  CaloE_emS25 = new vector<Float_t>();
  CaloE_hadMax = new vector<Float_t>();  CaloE_hadS9 = new vector<Float_t>();
  Calo_emPos_R = new vector<Float_t>();  Calo_emPos_eta = new vector<Float_t>();  Calo_emPos_phi = new vector<Float_t>();
  Calo_hadPos_R = new vector<Float_t>();  Calo_hadPos_eta = new vector<Float_t>();  Calo_hadPos_phi = new vector<Float_t>();
  dEdx = new vector<Float_t>();  dEdxError = new vector<Float_t>(); 
  dEdx_numberOfSaturatedMeasurements = new vector<Int_t>();  dEdx_numberOfMeasurements = new vector<Int_t>();

  Tracks_Tree->Branch("pt",&pt);  Tracks_Tree->Branch("eta",&eta);  Tracks_Tree->Branch("phi",&phi);
  Tracks_Tree->Branch("chargeMinus",&chargeMinus); Tracks_Tree->Branch("isGlobalMu",&isGlobalMu);
  Tracks_Tree->Branch("whichGenPar",&whichGenPar); Tracks_Tree->Branch("whichGenId",&whichGenId);
  Tracks_Tree->Branch("Vertex_x",&Vertex_x);  Tracks_Tree->Branch("Vertex_y",&Vertex_y);  Tracks_Tree->Branch("Vertex_z",&Vertex_z);
  Tracks_Tree->Branch("isoR03sumPt",&isoR03sumPt);  Tracks_Tree->Branch("isoR03emEt",&isoR03emEt);  Tracks_Tree->Branch("isoR03hadEt",&isoR03hadEt);
  Tracks_Tree->Branch("isoR03hoEt",&isoR03hoEt);  Tracks_Tree->Branch("isoR03nJets",&isoR03nJets);  Tracks_Tree->Branch("isoR03nTracks",&isoR03nTracks);
  Tracks_Tree->Branch("isoR05sumPt",&isoR05sumPt);  Tracks_Tree->Branch("isoR05emEt",&isoR05emEt);  Tracks_Tree->Branch("isoR05hadEt",&isoR05hadEt);
  Tracks_Tree->Branch("isoR05hoEt",&isoR05hoEt);  Tracks_Tree->Branch("isoR05nJets",&isoR05nJets);  Tracks_Tree->Branch("isoR05nTracks",&isoR05nTracks);
  Tracks_Tree->Branch("isoemVetoEt",&isoemVetoEt);  Tracks_Tree->Branch("isohadVetoEt",&isohadVetoEt);  Tracks_Tree->Branch("isohoVetoEt",&isohoVetoEt);
  Tracks_Tree->Branch("TrkKink",&TrkKink); Tracks_Tree->Branch("GlbKink",&GlbKink); Tracks_Tree->Branch("TrkRecoChi2",&TrkRecoChi2);
  Tracks_Tree->Branch("CaloE_emMax",&CaloE_emMax);  Tracks_Tree->Branch("CaloE_emS9",&CaloE_emS9);  Tracks_Tree->Branch("CaloE_emS25",&CaloE_emS25);
  Tracks_Tree->Branch("CaloE_hadMax",&CaloE_hadMax);  Tracks_Tree->Branch("CaloE_hadS9",&CaloE_hadS9);
  Tracks_Tree->Branch("Calo_emPos_R",&Calo_emPos_R);  Tracks_Tree->Branch("Calo_emPos_eta",&Calo_emPos_eta);  Tracks_Tree->Branch("Calo_emPos_phi",&Calo_emPos_phi);
  Tracks_Tree->Branch("Calo_hadPos_R",&Calo_hadPos_R);  Tracks_Tree->Branch("Calo_hadPos_eta",&Calo_hadPos_eta);  Tracks_Tree->Branch("Calo_hadPos_phi",&Calo_hadPos_phi);
  Tracks_Tree->Branch("dEdx",&dEdx);  Tracks_Tree->Branch("dEdxError",&dEdxError);
  Tracks_Tree->Branch("dEdx_numberOfSaturatedMeasurements",dEdx_numberOfSaturatedMeasurements);  Tracks_Tree->Branch("dEdx_numberOfMeasurements",&dEdx_numberOfMeasurements);

  //PV
  vx = new vector<Float_t>(); vy = new vector<Float_t>(); vz = new vector<Float_t>();
  vxError = new vector<Float_t>(); vyError =  new vector<Float_t>(); vzError = new vector<Float_t>();
  Tracks_Tree->Branch("PVx",&vx);  Tracks_Tree->Branch("PVy",&vy);  Tracks_Tree->Branch("PVz",&vz);
  Tracks_Tree->Branch("PVxError",&vxError);  Tracks_Tree->Branch("PVyError",&vyError);  Tracks_Tree->Branch("PVzError",&vzError);
  //Generation Level Muons
  Gen_pt = new vector<Float_t>();  Gen_eta = new vector<Float_t>();  Gen_phi = new vector<Float_t>();    Gen_pdgId = new vector<Int_t>();
  Gen_vx = new vector<Float_t>();  Gen_vy = new vector<Float_t>();  Gen_vz = new vector<Float_t>();   Gen_vt = new vector<Float_t>();
  Gen_charge=new vector<Byte_t>();
  Tracks_Tree->Branch("Gen_pt",&Gen_pt);  Tracks_Tree->Branch("Gen_eta",&Gen_eta);  Tracks_Tree->Branch("Gen_phi",&Gen_phi);   Tracks_Tree->Branch("Gen_pdgId",&Gen_pdgId);
  Tracks_Tree->Branch("Gen_vx",&Gen_vx);  Tracks_Tree->Branch("Gen_vy",&Gen_vy);  Tracks_Tree->Branch("Gen_vz",&Gen_vz);  Tracks_Tree->Branch("Gen_vt",&Gen_vt);
  Tracks_Tree->Branch("Gen_charge",&Gen_charge);
  //Simulated Tracks
  TrkParticles_pt = new vector<Float_t>();  TrkParticles_eta = new vector<Float_t>(); TrkParticles_phi = new vector<Float_t>();
  TrkParticles_pdgId = new vector<Int_t>();  DChains = new vector<Int_t>();
  Tracks_Tree->Branch("TrkParticles_pt",&TrkParticles_pt);  Tracks_Tree->Branch("TrkParticles_eta",&TrkParticles_eta);  Tracks_Tree->Branch("TrkParticles_phi",&TrkParticles_phi);
  Tracks_Tree->Branch("TrkParticles_pdgId",&TrkParticles_pdgId);
  Tracks_Tree->Branch("DChains",&DChains);

  //Summary
  ErrorMsg_Tree->Branch("ErrorMsg",&Error.ErrorCode,"ErrorCode/b:RunNum/l:EventNum");
  Summarization_Tree->Branch("Summarization",&Summarization.Total_Events,"Total_Events/l:Total_RecoMuons:Total_GenMuons:Total_GenElectrons:Total_GenPhotons");
  Summarization.Total_Events=0;Summarization.Total_RecoMuons=0;Summarization.Total_GenMuons=0;Summarization.Total_GenElectrons=0;Summarization.Total_GenPhotons=0;
  FirstEntry=true;
}


CSCPriEff::~CSCPriEff()
{
   //Summarizations
   Summarization_Tree->Fill();
   file->Write();
   file->Close();
}

void CSCPriEff::SimTrackDaughtersTree(edm::SimVertexContainer & SVC, edm::SimTrackContainer & STC, vector< vector<int> > & daughters, unsigned int trackpos)
{
  DChain.push_back(particlenum);
  particlenum++;
  GenSimMomentum(STC[trackpos].momentum().pt(),STC[trackpos].momentum().eta(),STC[trackpos].momentum().phi(),STC[trackpos].type(),(Byte_t) STC[trackpos].charge())
  if (!STC[trackpos].noVertex())
    {
      SimVertex thisVtx=SVC[STC[trackpos].vertIndex()];
      printf("%d (%5.5f,%5.5f,%5.5f,%fus) --> ",STC[trackpos].type(),thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9);
      GenSimVertex(thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9)     
    }
  else 
    {
      printf("%d --> ",STC[trackpos].type());
      GenSimVertex(0,0,0,0)
    }

  if (daughters[trackpos].size()>0)
    {
      SimTrackDaughtersTree(SVC,STC,daughters,daughters[trackpos].at(0));
      for (unsigned int i=1;i<daughters[trackpos].size();i++)
	{
	  SimVertex thisVtx=SVC[STC[trackpos].vertIndex()];
	  printf("%d (%5.5f,%5.5f,%5.5f,%fus) --> ",STC[trackpos].type(),thisVtx.position().x(),thisVtx.position().y(),thisVtx.position().z(),thisVtx.position().t()*1e9);
	  SimTrackDaughtersTree(SVC,STC,daughters,daughters[trackpos].at(i));
	}
    }
  else 
    {
      DChains->push_back(-1);
      DChains->insert(DChains->end(),DChain.begin(),DChain.end());
      printf("\n");
    }
  DChain.pop_back();
}

void CSCPriEff::HepMCParentTree(HepMC::GenParticle *genPar)
{
  printf("%d ",genPar->pdg_id());
  HepMC::GenVertex *thisVtx = genPar->production_vertex();
  if (thisVtx) 
    {
      printf("(%5.5f,%5.5f,%5.5f,%fus) <--| ",thisVtx->position().x()/10.,thisVtx->position().y()/10.,thisVtx->position().z()/10.,thisVtx->position().t()/299.792458);//change unit from mm to cm, mm to us
      for (HepMC::GenVertex::particles_in_const_iterator pgenD = thisVtx->particles_in_const_begin(); pgenD != thisVtx->particles_in_const_end(); ++pgenD)
	if (genPar->pdg_id()!=92)//Pythia special code for string, we only care about the particles after hadronization
	  {
	    if (pgenD == thisVtx->particles_in_const_begin()) HepMCParentTree(*pgenD); 
	    else
	      {
		printf("\n %d <--|",genPar->pdg_id());
		HepMCParentTree(*pgenD);
	      }
	  }
    }
}

void CSCPriEff::HepMCDaughtersTree(HepMC::GenParticle *genPar)
{
  cout<<genPar->pdg_id()<<" |--> ";
  HepMC::GenVertex *thisVtx = genPar->end_vertex();
  if (thisVtx) 
    {
      for (HepMC::GenVertex::particles_out_const_iterator pgenD = thisVtx->particles_out_const_begin(); pgenD != thisVtx->particles_out_const_end(); ++pgenD)
	if (pgenD == thisVtx->particles_out_const_begin()) HepMCDaughtersTree(*pgenD); 
	else
	  {
	    printf("\n %d |--> ",genPar->pdg_id());
	    HepMCDaughtersTree(*pgenD);
	  }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCPriEff);
