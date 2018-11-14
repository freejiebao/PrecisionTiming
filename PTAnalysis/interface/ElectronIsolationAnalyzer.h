// -*- C++ -*-
//
// Package:    PrecisionTiming/ElectronIsolationAnalyzer
// Class:      ElectronIsolationAnalyzer
//
/**\class ElectronIsolationAnalyzer ElectronIsolationAnalyzer.cc PrecisionTiming/ElectronIsolationAnalyzer/plugins/ElectronIsolationAnalyzer.cc

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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "TTree.h"
#include <vector>

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
using namespace math;

struct eventInfo {
    int           npu;
    vector<float> track_pt;
    vector<float> track_eta;
    vector<float> track_phi;
    vector<float> track_dz;
    vector<float> track_dz3D;
    vector<float> track_dz4D;
    vector<float> track_t;

    float vtxGen_z;
    float vtxGen_t;
    float vtx_z;
    float vtx_zErr;
    float vtx_t;
    float vtx_tErr;
    float vtx3D_z;
    float vtx3D_zErr;
    float vtx4D_z;
    float vtx4D_zErr;
    float vtx4D_t;
    float vtx4D_tErr;
    float vtx_isFake;
    float vtx3D_isFake;
    float vtx4D_isFake;

    // -- store the dr between the electron and pfCand, the use it to find the dr of veto cone.
    vector<float> drep;
    vector<float> electron_pt;
    vector<float> electron_eta;
    vector<float> electron_phi;
    vector<float> electron_sigmaIetaIeta;
    vector<float> electron_dz;
    vector<float> electron_dxy;
    vector<float> electron_dz3D;
    vector<float> electron_dxy3D;
    vector<float> electron_dz4D;
    vector<float> electron_dxy4D;
    vector<float> electron_t;
    vector<bool>  electron_isPrompt;
    vector<bool>  electron_isMatchedToGenJet;
    vector<bool>  electron_isMatchedToGenJet2;
    vector<float> electron_r9;
    vector<float> electron_chIso[10];
    vector<float> electron_chIso_dT[10][10];
    vector<float> electron_chIso_simVtx[10];
    vector<float> electron_chIso_dT_simVtx[10][10];
    vector<float> electron_chIso_dT_4D[10][10];
    vector<float> electron_chIso_reldZ[10];
    vector<float> electron_chIso_reldZ_dT[10][10];
    vector<float> electron_chIso_reldZ_dT_4D[10][10];
    vector<int>   track_elecIndex;

    /*vector<int> passVetoId;
    vector<int> passLooseId;
    vector<int> passMediumId;
    vector<int> passTightId;*/
    //vector<int> passHEEPId;
};

class ElectronIsolationAnalyzer : public edm::EDAnalyzer {
public:
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    explicit ElectronIsolationAnalyzer(const edm::ParameterSet&);
    ~ElectronIsolationAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob();

    void initEventStructure();

    //---inputs
    EDGetTokenT<genXYZ>                       genXYZToken_;
    EDGetTokenT<float>                        genT0Token_;
    EDGetTokenT<vector<PileupSummaryInfo>>    PileUpToken_;
    EDGetTokenT<View<reco::Vertex>>           vertexToken3D_;
    EDGetTokenT<View<reco::Vertex>>           vertexToken4D_;
    EDGetTokenT<edm::View<reco::PFCandidate>> pfcandToken_;
    EDGetTokenT<View<reco::GenParticle>>      genPartToken_;
    EDGetTokenT<vector<SimVertex>>            genVertexToken_;
    EDGetTokenT<View<reco::GenJet>>           genJetsToken_;
    EDGetTokenT<View<reco::GsfElectron>>      barrelElectronsToken_;
    //EDGetTokenT<View<reco::GsfElectron>>      endcapElectronsToken_;
    // ID decisions objects
    EDGetTokenT<edm::ValueMap<bool>> eleVetoIdMapToken_;
    EDGetTokenT<edm::ValueMap<bool>> eleLooseIdMapToken_;
    EDGetTokenT<edm::ValueMap<bool>> eleMediumIdMapToken_;
    EDGetTokenT<edm::ValueMap<bool>> eleTightIdMapToken_;
    //--- outputs
    edm::Service<TFileService> fs_;
    TTree*                     eventTree[10];
    eventInfo                  evInfo[10];

    //--- options
    vector<double> timeResolutions_;
    vector<double> isoConeDR_;
    bool           saveTracks_;
    float          maxDz_;
    float          minDr_;
    bool           isAOD_;
};

bool  isPromptElectron(const reco::GsfElectron& electron, const edm::View<reco::GenParticle>& genParticles);
bool  isMatchedToGenJet(const reco::GsfElectron& electron, const edm::View<reco::GenJet>& genJet);
bool  isMatchedToGenJet2(const reco::GsfElectron& electron, const edm::View<reco::GenJet>& genJets);
float Get_dEtaInSeed(const reco::GsfElectron& ele);
float Get_epCut(const reco::Candidate& ele);
