// -*- C++ -*-
//
// Package:    PrecisionTiming/PTAnalysis
// Class:      TimePUJetIdAnalyzer
// 
/**\class TimePUJetIdAnalyzer TimePUJetIdAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/TimePUJetIdAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Martina Malberti
//         Created:  Mon, 10 Oct 2016 14:06:02 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

////
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <vector>

#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


using namespace std;
using namespace edm;
using namespace reco;


struct eventInfo
{

  int evId;
  int npu;
  float vtx4D_t;
  float vtx4D_z;
  float vtx_z;
  float genvtx_z;

  float jetPt;
  float jetEta;
  float jetPhi;
  float jetDR2;
  float jetDR2cleaned;
  float jetDT2n;
  float jetDT2c;
  float jetDT2;
  bool  jetIsMatchedToGen;
  float jetGenPt;
  float jetGenEta;
  float jetGenPhi;
  float jetDrToGen;

  vector<float> chTime;
  vector<float> chPt;
  vector<float> chEta;
  vector<float> chPhi;

  vector<float> clusTime;
  vector<float> clusPt;
  vector<float> clusEta;
  vector<float> clusPhi;

};


class TimePUJetIdAnalyzer : public edm::EDAnalyzer  
{
public:
  explicit TimePUJetIdAnalyzer(const edm::ParameterSet&);
  ~TimePUJetIdAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
    
  void initEventStructure();

      
  //---inputs  
  EDGetTokenT<vector<PileupSummaryInfo> > PileUpToken_;
  EDGetTokenT<View<reco::Vertex> > vertexToken_;
  EDGetTokenT<View<reco::Vertex> > vertex4DToken_;
  EDGetTokenT<View<reco::PFJet> > chsjetsToken_;
  EDGetTokenT<View<reco::GenJet> > genjetsToken_;
  EDGetTokenT<View<reco::Muon> > muonsToken_;
  EDGetTokenT<View<reco::GenParticle> > genparticlesToken_;
  EDGetTokenT<ValueMap<float> > trackTimeToken_;

  

  //--- outputs
  edm::Service<TFileService> fs_;
  TTree *eventTree;
  eventInfo evInfo;

};
