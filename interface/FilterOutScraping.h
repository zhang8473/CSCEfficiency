// Modified from
// CMSSW/DPGAnalysis/Skims/interface/FilterOutScraping.h
// CMSSW/DPGAnalysis/Skims/src/FilterOutScraping.cc

// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//


class FilterOutScraping {
public:
  explicit FilterOutScraping( const edm::ParameterSet & );
  virtual bool filter ( const edm::Event & );
  ~FilterOutScraping();
  
private:
  double thresh;
  unsigned int numtrack;
  edm::InputTag tracks_;

  reco::TrackBase::TrackQuality _trackQuality;
};

FilterOutScraping::FilterOutScraping(const edm::ParameterSet& iConfig)
{
  
  thresh =  iConfig.getUntrackedParameter<double>("thresh",0.2);
  numtrack = iConfig.getUntrackedParameter<unsigned int>("numtrack",10);
  tracks_ = iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("generalTracks"));
}

FilterOutScraping::~FilterOutScraping()
{
}

bool FilterOutScraping::filter( const edm::Event& iEvent)
{
  bool accepted = false;
  float fraction = 0;  
  // get GeneralTracks collection

  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel(tracks_,tkRef);    
  const reco::TrackCollection* tkColl = tkRef.product();

  //std::cout << "Total Number of Tracks " << tkColl->size() << std::endl;
  
  int numhighpurity=0;
  _trackQuality = reco::TrackBase::qualityByName("highPurity");

  if(tkColl->size()>numtrack){ 
    reco::TrackCollection::const_iterator itk = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();
    for(;itk!=itk_e;++itk){
      // std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
      if(itk->quality(_trackQuality)) numhighpurity++;
    }
    fraction = (float)numhighpurity/(float)tkColl->size();
    if(fraction>thresh) accepted=true;
  }else{
    //if less than 10 Tracks accept the event anyway    
    accepted= true;
  }
  return accepted;
}
