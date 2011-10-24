// Modified from
// CMSSW/DPGAnalysis/Skims/src/GoodVertexFilter.cc

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// class declaration
//

class GoodVertexFilter {
   public:
      explicit GoodVertexFilter(const edm::ParameterSet&);
      virtual bool filter(const edm::Event&, const edm::EventSetup&);
      ~GoodVertexFilter();

   private:
      edm::InputTag vertexSrc;        
      unsigned int minNDOF;
      double maxAbsZ;
      double maxd0;
};

GoodVertexFilter::GoodVertexFilter(const edm::ParameterSet& iConfig)
{
  vertexSrc = iConfig.getParameter<edm::InputTag>("vertexCollection");
  minNDOF = iConfig.getParameter<unsigned int>("minimumNDOF");
  maxAbsZ = iConfig.getParameter<double>("maxAbsZ");
  maxd0 = iConfig.getParameter<double>("maxd0");

}


GoodVertexFilter::~GoodVertexFilter() {}

bool
GoodVertexFilter::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 bool result = false; 
 edm::Handle<reco::VertexCollection> pvHandle; 
 iEvent.getByLabel(vertexSrc,pvHandle);
 const reco::VertexCollection & vertices = *pvHandle.product();
 for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
  {
      if(it->ndof() > minNDOF && 
         ( (maxAbsZ <=0 ) || fabs(it->z()) <= maxAbsZ ) &&
         ( (maxd0 <=0 ) || fabs(it->position().rho()) <= maxd0 )
       ) result = true;
  }

 return result;
}
