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

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
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
#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "TTree.h"
#include <TRandom.h>
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
    vector<float> track_dz4D;
    vector<float> track_dz3D;
    vector<float> track_dxy3D;
    vector<float> track_dxy4D;
    vector<float> track_t;
    vector<int>   track_elecIndex;

    float  vtxGen_z;
    float  vtxGen_t;
    float  vtx4D_z;
    float  vtx4D_zErr;
    float  vtx4D_t;
    float  vtx4D_tErr;
    float  vtx3D_z;
    float  vtx3D_zErr;
    int    vtx4D_isFake;
    int    vtx3D_isFake;
    double rho;
    double rho_calo;
    // -- info of electron work point
    vector<double> electron_sigmaIetaIeta;
    vector<double> electron_dEtaInSeed;
    vector<double> electron_dPhiIn;
    vector<double> electron_hoe;
    vector<double> electron_energy_sc;
    vector<double> electron_pf_isolation;
    vector<double> electron_pf_isolation_calo;
    vector<double> electron_ooEmooP;
    vector<double> electron_mHits;
    vector<bool>   electron_pass_conversion_veto;
    vector<double> electronsc_eta;
    // -- store the dr between the electron and pfCand, the use it to find the dr of veto cone.
    vector<float> drep;
    vector<float> electron_pt;
    vector<float> electron_eta;
    vector<float> electron_phi;
    vector<float> electron_dz4D;
    vector<float> electron_dxy4D;
    vector<float> electron_dz3D;
    vector<float> electron_dxy3D;
    vector<float> electron_t;
    vector<int>   electron_isPrompt;
    vector<int>   electron_isMatchedToGenJet;
    vector<int>   electron_isMatchedToGenJet2;
    vector<int>   electron_isFromTauDecay;
    vector<float> electron_r9;

    vector<float> electron_chIso_dZ05_simVtx;
    vector<float> electron_chIso_dZ05_dT_simVtx;

    vector<float> electron_chIso_dZ1_simVtx;
    vector<float> electron_chIso_dZ1_dT_simVtx;

    vector<float> electron_chIso_dZ2_simVtx;
    vector<float> electron_chIso_dZ2_dT_simVtx;

    vector<float> electron_chIso_dZ05;
    vector<float> electron_chIso_dZ05_dT;

    vector<float> electron_chIso_dZ1;
    vector<float> electron_chIso_dZ1_dT;

    vector<float> electron_chIso_dZ2;
    vector<float> electron_chIso_dZ2_dT;

    vector<float> electron_chIso_reldZ;
    vector<float> electron_chIso_reldZ_dT;

    vector<float> electron_chIso_dZele05;
    vector<float> electron_chIso_dZele05_dTele;

    vector<float> electron_chIso_dZele1;
    vector<float> electron_chIso_dZele1_dTele;

    vector<float> electron_chIso_dZele2;
    vector<float> electron_chIso_dZele2_dTele;
    /*
    vector<int> passVetoId;
    vector<int> passLooseId;
    vector<int> passMediumId;
    vector<int> passTightId;
    */
};

class ElectronIsolationAnalyzer : public edm::EDAnalyzer {
public:
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    explicit ElectronIsolationAnalyzer(const edm::ParameterSet&);
    ~ElectronIsolationAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
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
    EDGetTokenT<View<reco::GsfElectron>>      endcapElectronsToken_;
    EDGetTokenT<double>                       RhoToken_;
    EDGetTokenT<double>                       RhoCaloToken_;
    EffectiveAreas                            effectiveAreas_;
    //EDGetTokenT<View<reco::Conversion>>       convsToken_;
    //EDGetTokenT<reco::BeamSpot>               thebsToken_;
    edm::Handle<reco::ConversionCollection> _convs;
    edm::Handle<reco::BeamSpot>             _thebs;
    /*
    // ID decisions objects
    EDGetTokenT<edm::ValueMap<bool>> eleVetoIdMapToken_;
    EDGetTokenT<edm::ValueMap<bool>> eleLooseIdMapToken_;
    EDGetTokenT<edm::ValueMap<bool>> eleMediumIdMapToken_;
    EDGetTokenT<edm::ValueMap<bool>> eleTightIdMapToken_;
    */
    //--- outputs
    edm::Service<TFileService> fs_;
    TTree*                     eventTree[10];
    eventInfo                  evInfo[10];

    //--- options
    vector<double> timeResolutions_;
    double         isoConeDR_;
    bool           saveTracks_;
    double         maxDz_;
    double         minDr_;
    double         minTrackPt_;
    bool           useVertexClosestToGen_;
    bool           isAOD_;
};

bool  isPromptElectron(const reco::GsfElectron& electron, const edm::View<reco::GenParticle>& genParticles);
bool  isMatchedToGenJet(const reco::GsfElectron& electron, const edm::View<reco::GenJet>& genJet);
bool  isMatchedToGenJet2(const reco::GsfElectron& electron, const edm::View<reco::GenJet>& genJets);
bool  isFromTau(const reco::GsfElectron& electron, const edm::View<reco::GenParticle>& genParticles);
float Get_dEtaInSeed(const reco::GsfElectron& ele);
float Get_epCut(const reco::Candidate& ele);
