// -*- C++ -*-%
//
// Package:    PrecisionTiming/ElectronIsolationAnalyzer
// Class:      ElectronIsolationAnalyzer
//
/**\class ElectronIsolationAnalyzer ElectronIsolationAnalyzer.cc PrecisionTiming/PTAnalysis/plugins/ElectronIsolationAnalyzer.cc

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
#include <iostream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PrecisionTiming/PTAnalysis/interface/ElectronIsolationAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include <TMath.h>
#include <TRandom.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//
// constructors and destructor
//
ElectronIsolationAnalyzer::ElectronIsolationAnalyzer(const edm::ParameterSet& iConfig)
    : genXYZToken_(consumes<genXYZ>(iConfig.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
      genT0Token_(consumes<float>(iConfig.getUntrackedParameter<edm::InputTag>("genT0Tag"))), PileUpToken_(consumes<vector<PileupSummaryInfo>>(iConfig.getParameter<InputTag>("PileUpTag"))),
      vertexToken3D_(consumes<View<reco::Vertex>>(iConfig.getParameter<InputTag>("VertexTag3D"))),
      vertexToken4D_(consumes<View<reco::Vertex>>(iConfig.getParameter<InputTag>("VertexTag4D"))),
      pfcandToken_(consumes<View<reco::PFCandidate>>(iConfig.getParameter<InputTag>("PFCandidateTag"))),
      genPartToken_(consumes<View<reco::GenParticle>>(iConfig.getUntrackedParameter<InputTag>("genPartTag"))),
      genVertexToken_(consumes<vector<SimVertex>>(iConfig.getUntrackedParameter<InputTag>("genVtxTag"))),
      genJetsToken_(consumes<View<reco::GenJet>>(iConfig.getUntrackedParameter<InputTag>("genJetsTag"))),
      barrelElectronsToken_(consumes<View<reco::GsfElectron>>(iConfig.getUntrackedParameter<edm::InputTag>("barrelElectronsTag"))),
      endcapElectronsToken_(consumes<View<reco::GsfElectron>>(iConfig.getUntrackedParameter<edm::InputTag>("endcapElectronsTag"))),
      RhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho"))),
      RhoCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho_Calo"))),
      effectiveAreas_((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()) {
    timeResolutions_       = iConfig.getUntrackedParameter<vector<double>>("timeResolutions");
    isoConeDR_             = iConfig.getUntrackedParameter<double>("isoConeDR");
    saveTracks_            = iConfig.getUntrackedParameter<bool>("saveTracks");
    maxDz_                 = iConfig.getUntrackedParameter<double>("maxDz");
    minDr_                 = iConfig.getUntrackedParameter<double>("minDr");
    minTrackPt_            = iConfig.getUntrackedParameter<double>("minTrackPt");
    useVertexClosestToGen_ = iConfig.getUntrackedParameter<bool>("useVertexClosestToGen");
    isAOD_                 = iConfig.getUntrackedParameter<bool>("isAOD");
    //Now do what ever initialization is needed
    for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
        eventTree[iRes] = fs_->make<TTree>(Form("tree_%dps", int(timeResolutions_[iRes] * 1000)), Form("tree_%dps", int(timeResolutions_[iRes] * 1000)));
        cout << iRes << "  " << timeResolutions_[iRes] << "  " << eventTree[iRes]->GetName() << endl;
    }
    gRandom = new TRandom();
}

ElectronIsolationAnalyzer::~ElectronIsolationAnalyzer() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void ElectronIsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // -- get the vertex 3D collection
    Handle<View<reco::Vertex>> Vertex3DCollectionH;
    iEvent.getByToken(vertexToken3D_, Vertex3DCollectionH);
    const edm::View<reco::Vertex>& vertices3D = *Vertex3DCollectionH;

    // -- get the vertex 4D collection
    Handle<View<reco::Vertex>> Vertex4DCollectionH;
    iEvent.getByToken(vertexToken4D_, Vertex4DCollectionH);
    const edm::View<reco::Vertex>& vertices4D = *Vertex4DCollectionH;

    // -- get the PU
    Handle<vector<PileupSummaryInfo>> PileupInfos;
    if (!iEvent.isRealData()) {
        iEvent.getByToken(PileUpToken_, PileupInfos);
    }
    else
        return;

    // -- get the barrel electrons
    Handle<View<reco::GsfElectron>> BarrelElectronCollectionH;
    iEvent.getByToken(barrelElectronsToken_, BarrelElectronCollectionH);
    const edm::View<reco::GsfElectron>& barrelElectrons = *BarrelElectronCollectionH;

    // -- get the endcap electrons
    Handle<View<reco::GsfElectron>> EndcapElectronCollectionH;
    iEvent.getByToken(endcapElectronsToken_, EndcapElectronCollectionH);
    const edm::View<reco::GsfElectron>& endcapElectrons = *EndcapElectronCollectionH;

    // -- get the PFCandidate collection
    Handle<View<reco::PFCandidate>> PFCandidateCollectionH;
    iEvent.getByToken(pfcandToken_, PFCandidateCollectionH);
    const edm::View<reco::PFCandidate>& pfcands = *PFCandidateCollectionH;

    // -- get the gen particles collection
    Handle<View<reco::GenParticle>> GenParticleCollectionH;
    iEvent.getByToken(genPartToken_, GenParticleCollectionH);
    const edm::View<reco::GenParticle>& genParticles = *GenParticleCollectionH;

    // -- get the gen jets collection
    Handle<View<reco::GenJet>> GenJetCollectionH;
    iEvent.getByToken(genJetsToken_, GenJetCollectionH);
    const edm::View<reco::GenJet>& genJets = *GenJetCollectionH;

    // -- get the rho Value
    Handle<double> rhoH;
    iEvent.getByToken(RhoToken_, rhoH);
    double rho_ = rhoH.isValid() ? (float)(*rhoH) : 0;

    // -- get the rho_calo Value
    Handle<double> rhoCaloH;
    iEvent.getByToken(RhoCaloToken_, rhoCaloH);
    double rho_calo_ = rhoCaloH.isValid() ? (float)(*rhoCaloH) : 0;

    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    /*Handle<edm::ValueMap<bool>> veto_id_decisions;
    Handle<edm::ValueMap<bool>> loose_id_decisions;
    Handle<edm::ValueMap<bool>> medium_id_decisions;
    Handle<edm::ValueMap<bool>> tight_id_decisions;
    iEvent.getByToken(eleVetoIdMapToken_, veto_id_decisions);
    iEvent.getByToken(eleLooseIdMapToken_, loose_id_decisions);
    iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
    iEvent.getByToken(eleTightIdMapToken_, tight_id_decisions);*/

    SimVertex genPV;
    if (isAOD_) {
        // -- get the genXYZ
        Handle<genXYZ> genXYZH;
        iEvent.getByToken(genXYZToken_, genXYZH);
        auto xyz = genXYZH.product();  //

        // -- get the genT0
        Handle<float> genT0H;
        iEvent.getByToken(genT0Token_, genT0H);
        auto t = *genT0H.product();

        //---get truth PV from the genXYZ and genT0
        auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
        genPV  = SimVertex(v, t);
    }
    else {
        // -- get the gen vertex collection
        Handle<vector<SimVertex>> GenVertexCollectionH;
        iEvent.getByToken(genVertexToken_, GenVertexCollectionH);
        const vector<SimVertex>& genVertices = *GenVertexCollectionH;

        //-- get truth PV
        genPV = genVertices.at(0);
    }

    // -- initialize output tree
    initEventStructure();

    // -- number of pileup events
    int nPU = 0;
    if (!iEvent.isRealData()) {
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for (PVI = PileupInfos->begin(); PVI != PileupInfos->end(); ++PVI) {
            Int_t pu_bunchcrossing = PVI->getBunchCrossing();
            if (pu_bunchcrossing == 0) {
                nPU = PVI->getPU_NumInteractions();
            }
        }
    }

    double mindz       = 999999.;
    int    pv_index_3D = -1;
    int    pv_index_4D = -1;

    // -- find the reco vertex closest to the gen vertex (3D)
    for (unsigned int ivtx = 0; ivtx < vertices3D.size(); ivtx++) {
        const reco::Vertex& vtx = vertices3D[ivtx];
        float               dz  = std::abs(vtx.z() - genPV.position().z());
        if (dz < mindz) {
            mindz       = dz;
            pv_index_3D = ivtx;
        }
    }

    // -- find the reco vertex closest to the gen vertex (4D)
    mindz = 999999.;
    for (unsigned int ivtx = 0; ivtx < vertices4D.size(); ivtx++) {
        const reco::Vertex& vtx = vertices4D[ivtx];
        float               dz  = std::abs(vtx.z() - genPV.position().z());
        if (dz < mindz) {
            mindz       = dz;
            pv_index_4D = ivtx_4D;
        }
    }

    // -- if PV index == -1, use highest ranked vertex
    if (pv_index_3D == -1)
        pv_index_3D = 0;
    if (pv_index_4D == -1)
        pv_index_4D = 0;
    // -- get isolation around a candidate electron
    // --- using only vtx closest to gen vtx
    const reco::Vertex& vtx3D = vertices3D[pv_index_3D];
    const reco::Vertex& vtx4D = vertices4D[pv_index_4D];

    int       electronIndex = -1;
    const int nResol        = timeResolutions_.size();
    // -- start loop over barrel electrons
    for (unsigned int iele = 0; iele < barrelElectrons.size(); iele++) {

        const reco::GsfElectron& electron = barrelElectrons[iele];

        // -- minimal checks
        if (electron.pt() < 10.)
            continue;
        if (electron.gsfTrack().isNull())
            continue;
        //if (fabs(electron.eta()) > 1.5)
        //    continue;
        electronIndex++;
        // -- check if prompt or fake electron
        bool isPromptEle    = isPromptElectron(electron, genParticles);
        bool isMatchedJet   = isMatchedToGenJet(electron, genJets);
        bool isMatchedJet2  = isMatchedToGenJet2(electron, genJets);
        bool isFromTauDecay = isFromTau(electron, genParticles);
        // -- Look up and save the ID decisions
        /*float sieieCut          = 0.0112;
        float dEtaSeedCut       = 0.00377;
        float dPhiInCut         = 0.0884;
        float hoeCut            = Get_hoeCut();
        float relIsoWithEACut   = 0.133;
        float epCut             = 0.129;
        int   mHitsCut          = 1;
        bool  conversionVetoCut = TRUE;

        float          sieie        = electron.full5x5_sigmaIetaIeta();
        float          dEtaSeed     = Get_dEtaSeed(electron);
        float          dPhiIn       = electron.deltaPhiSuperClusterTrackAtVtx();
        float          hoe          = Get_hoe(electron);
        float          relIsoWithEA = 0.112 + 0.506 / electron.pt();
        float          ep           = Get_epCut(electron);
        int            mHits        = Get_mHitsCut(electron);
        bool           conversion float if (< 0.0105 &&);*/
        /*auto el            = barrelElectrons.ptrAt(iele);
        bool isPassVeto_   = (*veto_id_decisions)[el];
        bool isPassLoose_  = (*loose_id_decisions)[el];
        bool isPassMedium_ = (*medium_id_decisions)[el];
        bool isPassTight_  = (*tight_id_decisions)[el];*/

        // -- compute charged isolations
        float chIso_dZ05_simVtx = 0.;
        float chIso_dZ1_simVtx  = 0.;
        float chIso_dZ2_simVtx  = 0.;
        float chIso_dZ05        = 0.;
        float chIso_dZ1         = 0.;
        float chIso_dZ2         = 0.;
        float chIso_reldZ       = 0.;
        float chIso_dZele05     = 0.;
        float chIso_dZele1      = 0.;
        float chIso_dZele2      = 0.;

        float chIso_dZ05_dT_simVtx[nResol];
        float chIso_dZ1_dT_simVtx[nResol];
        float chIso_dZ2_dT_simVtx[nResol];
        float chIso_dZ05_dT[nResol];
        float chIso_dZ1_dT[nResol];
        float chIso_dZ2_dT[nResol];
        float chIso_reldZ_dT[nResol];
        float chIso_dZele05_dTele[nResol];
        float chIso_dZele1_dTele[nResol];
        float chIso_dZele2_dTele[nResol];
        float time[nResol];

        // -- initialize
        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
            chIso_dZ05_dT_simVtx[iRes] = 0.;
            chIso_dZ1_dT_simVtx[iRes]  = 0.;
            chIso_dZ2_dT_simVtx[iRes]  = 0.;
            chIso_dZ05_dT[iRes]        = 0.;
            chIso_dZ1_dT[iRes]         = 0.;
            chIso_dZ2_dT[iRes]         = 0.;
            chIso_reldZ_dT[iRes]       = 0.;
            chIso_dZele05_dTele[iRes]  = 0.;
            chIso_dZele1_dTele[iRes]   = 0.;
            chIso_dZele2_dTele[iRes]   = 0.;
            time[iRes]                 = 0.;
        }

        // -- firstly, get the electron time
        float        elecTime = 0.;  // record time of electron
        unsigned int icandTag = -1;
        float        dr_      = 99999999.;
        for (unsigned icand = 0; icand < pfcands.size(); ++icand) {
            const reco::PFCandidate& pfcand = pfcands[icand];
            if (pfcand.charge() == 0)
                continue;

            // -- get the track ref
            auto           pfcandRef = pfcands.refAt(icand);
            reco::TrackRef trackRef  = pfcandRef->trackRef();
            if (trackRef.isNull())
                continue;
            if (!trackRef->quality(reco::TrackBase::highPurity))
                continue;
            float dr = deltaR(electron.eta(), electron.phi(), pfcand.eta(), pfcand.phi());
            if (dr < dr_) {
                dr_      = dr;
                icandTag = icand;
                //-- get electron time from pfcand
                elecTime = pfcand.time();
            }
            if (isoConeDR_ == 0.3 && saveTracks_) {
                for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                    evInfo[iRes].drep.push_back(dr);
                }
            }
        }
        // -- firstly, get the electron time

        // -- compute loose electron ID
        //
        // -- full5x5_sigmaIetaIeta
        double sigmaIEtaIEta = electron.full5x5_sigmaIetaIeta();
        // -- abs(dEtaSeed)
        double dEtaInSeed = fabs(electron.superCluster().isNonnull() && electron.superCluster()->seed().isNonnull() ? electron.deltaEtaSuperClusterTrackAtVtx() - electron.superCluster()->eta() + electron.superCluster()->seed()->eta() : std::numeric_limits<float>::max());
        // -- abs(dPhiIn)
        double dPhiIn = fabs(electron.deltaPhiSuperClusterTrackAtVtx());
        // -- H/E
        double hoe       = electron.hadronicOverEm();
        double energy_sc = electron.superCluster()->energy();
        // -- relIsoWithEA
        float eA                = effectiveAreas_.getEffectiveArea(fabs(electron.superCluster()->eta()));
        float pf_isolation      = ((electron.pfIsolationVariables().sumChargedHadronPt + std::max(0., electron.pfIsolationVariables().sumNeutralHadronEt + electron.pfIsolationVariables().sumPhotonEt - eA * rho_)) / electron.pt());
        float pf_isolation_calo = ((electron.pfIsolationVariables().sumChargedHadronPt + std::max(0., electron.pfIsolationVariables().sumNeutralHadronEt + electron.pfIsolationVariables().sumPhotonEt - eA * rho_calo_)) / electron.pt());
        // --  EleEInverseMinusPInverseCut
        double ooEmooP = fabs((1.0 / electron.ecalEnergy() - electron.eSuperClusterOverP() / electron.ecalEnergy()));
        // -- expected missing inner hits
        double mHits = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
        //pass conversion veto
        int pass_conversion_veto = electron.passConversionVeto();

        // -- loop over charged pf candidates
        for (unsigned icand = 0; icand < pfcands.size(); ++icand) {
            const reco::PFCandidate& pfcand = pfcands[icand];
            if (pfcand.charge() == 0)
                continue;
            // -- get the track ref
            auto           pfcandRef = pfcands.refAt(icand);
            reco::TrackRef trackRef  = pfcandRef->trackRef();
            if (trackRef.isNull())
                continue;
            if (!trackRef->quality(reco::TrackBase::highPurity))
                continue;
            //-- use tracks with pT above thtreshold
            if (pfcand.pt() < minTrackPt_)
                continue;
            // -- get dz, dxy
            float dz3D  = std::abs(trackRef->dz(vtx3D.position()));
            float dz4D  = std::abs(trackRef->dz(vtx4D.position()));
            float dxy3D = std::abs(trackRef->dxy(vtx3D.position()));
            float dxy4D = std::abs(trackRef->dxy(vtx4D.position()));

            float dz3Drel = std::abs(dz3D / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx3D.zError() * vtx3D.zError()));
            float dz4Drel = std::abs(dz4D / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx4D.zError() * vtx4D.zError()));

            float dzsim  = std::abs(trackRef->vz() - genPV.position().z());
            float dxysim = sqrt(pow(trackRef->vx() - genPV.position().x(), 2) + pow(trackRef->vy() - genPV.position().y(), 2));

            float dzele  = std::abs(trackRef->dz(vtx4D.position()) - electron.gsfTrack()->dz(vtx4D.position()));
            float dxyele = std::abs(trackRef->dxy(vtx4D.position()) - electron.gsfTrack()->dxy(vtx4D.position()));

            float dr = deltaR(electron.eta(), electron.phi(), pfcand.eta(), pfcand.phi());

            // --- no timing
            if (dr > minDr_ && dr < isoConeDR_) {
                // -- sim vertex
                if (dzsim < 0.05 && dxysim < 0.02) {
                    chIso_dZ05_simVtx += pfcand.pt();
                }
                if (dzsim < 0.1 && dxysim < 0.02) {
                    chIso_dZ1_simVtx += pfcand.pt();
                }
                if (dzsim < 0.2 && dxysim < 0.02) {
                    chIso_dZ2_simVtx += pfcand.pt();
                }

                // -- reco vtx closest to the sim one
                if (dz3D < 0.05 && dxy3D < 0.02) {
                    chIso_dZ05 += pfcand.pt();
                }
                if (dz3D < 0.1 && dxy3D < 0.02) {
                    chIso_dZ1 += pfcand.pt();
                }
                if (dz3D < 0.2 && dxy3D < 0.02) {
                    chIso_dZ2 += pfcand.pt();
                }

                // -- using reco vtx closest to the sim one and cut on relative dz
                if (dz3Drel < 3.0 && dxy3D < 0.02) {
                    chIso_reldZ += pfcand.pt();
                }

                // -- dz wrt to electron
                if (dzele < 0.05) {
                    chIso_dZele05 += pfcand.pt();
                }
                if (dzele < 0.1) {
                    chIso_dZele1 += pfcand.pt();
                }
                if (dzele < 0.2) {
                    chIso_dZele2 += pfcand.pt();
                }
            }

            // --- with timing
            if (dr > minDr_ && dr < isoConeDR_) {
                for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                    double targetTimeResol  = timeResolutions_[iRes];
                    double defaultTimeResol = 0.;
                    if (pfcand.isTimeValid())
                        defaultTimeResol = double(pfcand.timeError());
                    double extra_resol = sqrt(targetTimeResol * targetTimeResol - defaultTimeResol * defaultTimeResol);
                    double dtsim       = 0.;
                    double dt          = 0.;
                    double dtele       = 0.;
                    time[iRes]         = -999.;
                    if (pfcand.isTimeValid()) {
                        double rnd = gRandom->Gaus(0., extra_resol);
                        //cout << "target time resol = "<< targetTimeResol << "  extra_resol = "<< extra_resol << "  rnd = " << rnd <<endl;
                        time[iRes] = pfcand.time() + rnd;
                        dtsim      = std::abs(time[iRes] - genPV.position().t() * 1000000000.);
                        cout << "dasim:" << dtsim << endl;
                        dt    = std::abs(time[iRes] - vtx4D.t());
                        dtele = std::abs(time[iRes] - elecTime);
                    }
                    else {
                        dtsim = 0.;
                        dt    = 0.;
                        dtele = 0.;
                    }

                    // -- sim vertex
                    if (dtsim < 3 * targetTimeResol && dzsim < 0.05 && dxysim < 0.02) {
                        chIso_dZ05_dT_simVtx[iRes] += pfcand.pt();
                    }
                    if (dtsim < 3 * targetTimeResol && dzsim < 0.1 && dxysim < 0.02) {
                        chIso_dZ1_dT_simVtx[iRes] += pfcand.pt();
                    }
                    if (dtsim < 3 * targetTimeResol && dzsim < 0.2 && dxysim < 0.02) {
                        chIso_dZ2_dT_simVtx[iRes] += pfcand.pt();
                    }

                    // -- reco vtx closest to the sim vtx
                    if (dt < 3 * targetTimeResol && dz4D < 0.05 && dxy4D < 0.02) {
                        chIso_dZ05_dT[iRes] += pfcand.pt();
                    }
                    if (dt < 3 * targetTimeResol && dz4D < 0.1 && dxy4D < 0.02) {
                        chIso_dZ1_dT[iRes] += pfcand.pt();
                    }
                    if (dt < 3 * targetTimeResol && dz4D < 0.2 && dxy4D < 0.02) {
                        chIso_dZ2_dT[iRes] += pfcand.pt();
                    }

                    // -- using reco vtx closest to the sim vtx and cut on relative dz
                    if (dt < 3 * targetTimeResol && dz4Drel < 3.0 && dxy4D < 0.02) {
                        chIso_reldZ_dT[iRes] += pfcand.pt();
                    }

                    // -- dT wrt to electron
                    if (dtele < sqrt(2) * 3 * targetTimeResol && dzele < 0.05) {
                        chIso_dZele05_dTele[iRes] += pfcand.pt();
                    }
                    if (dtele < sqrt(2) * 3 * targetTimeResol && dzele < 0.1) {
                        chIso_dZele1_dTele[iRes] += pfcand.pt();
                    }
                    if (dtele < sqrt(2) * 3 * targetTimeResol && dzele < 0.2) {
                        chIso_dZele2_dTele[iRes] += pfcand.pt();
                    }

                    // -- save info for tracks in the isolation cone (only for DR = 0.3)
                    if (saveTracks_ && (dz4D < 1.0 || dz3D < 1. || dzele < 1.)) {  // save a subset of tracks with loose dz selection
                        evInfo[iRes].track_t.push_back(time[iRes]);
                        evInfo[iRes].track_dz4D.push_back(trackRef->dz(vtx4D.position()));
                        evInfo[iRes].track_dz3D.push_back(trackRef->dz(vtx3D.position()));
                        evInfo[iRes].track_dxy4D.push_back(trackRef->dxy(vtx4D.position()));
                        evInfo[iRes].track_dxy3D.push_back(trackRef->dxy(vtx3D.position()));
                        evInfo[iRes].track_pt.push_back(pfcand.pt());
                        evInfo[iRes].track_eta.push_back(pfcand.eta());
                        evInfo[iRes].track_phi.push_back(pfcand.phi());
                        evInfo[iRes].track_elecIndex.push_back(electronIndex);
                    }
                }  // end loop over time resolutions
            }
        }  // end loop over tracks

        // -- fill electron info for each resolution scenario
        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
            evInfo[iRes].electron_pt.push_back(electron.pt());
            evInfo[iRes].electron_eta.push_back(electron.eta());
            evInfo[iRes].electronsc_eta.push_back(electron.superCluster()->eta());
            evInfo[iRes].electron_phi.push_back(electron.phi());
            // -- info of electron work point
            evInfo[iRes].electron_sigmaIetaIeta.push_back(sigmaIEtaIEta);
            evInfo[iRes].electron_dEtaInSeed.push_back(dEtaInSeed);
            evInfo[iRes].electron_dPhiIn.push_back(dPhiIn);
            evInfo[iRes].electron_hoe.push_back(hoe);
            evInfo[iRes].electron_energy_sc.push_back(energy_sc);
            evInfo[iRes].electron_pf_isolation.push_back(pf_isolation);
            evInfo[iRes].electron_pf_isolation_calo.push_back(pf_isolation_calo);
            evInfo[iRes].electron_ooEmooP.push_back(ooEmooP);
            evInfo[iRes].electron_mHits.push_back(mHits);
            evInfo[iRes].electron_pass_conversion_veto.push_back(pass_conversion_veto);

            evInfo[iRes].electron_dz4D.push_back(electron.gsfTrack()->dz(vtx4D.position()));
            evInfo[iRes].electron_dxy4D.push_back(electron.gsfTrack()->dxy(vtx4D.position()));
            evInfo[iRes].electron_dz3D.push_back(electron.gsfTrack()->dz(vtx3D.position()));
            evInfo[iRes].electron_dxy3D.push_back(electron.gsfTrack()->dxy(vtx3D.position()));
            evInfo[iRes].electron_t.push_back(elecTime);
            evInfo[iRes].electron_isPrompt.push_back(isPromptEle);
            evInfo[iRes].electron_isMatchedToGenJet.push_back(isMatchedJet);
            evInfo[iRes].electron_isMatchedToGenJet2.push_back(isMatchedJet2);
            evInfo[iRes].electron_isFromTauDecay.push_back(isFromTauDecay);
            evInfo[iRes].electron_r9.push_back(electron.r9());
            evInfo[iRes].electron_chIso_dZ05_simVtx.push_back(chIso_dZ05_simVtx);
            evInfo[iRes].electron_chIso_dZ05_dT_simVtx.push_back(chIso_dZ05_dT_simVtx[iRes]);

            evInfo[iRes].electron_chIso_dZ1_simVtx.push_back(chIso_dZ1_simVtx);
            evInfo[iRes].electron_chIso_dZ1_dT_simVtx.push_back(chIso_dZ1_dT_simVtx[iRes]);

            evInfo[iRes].electron_chIso_dZ2_simVtx.push_back(chIso_dZ2_simVtx);
            evInfo[iRes].electron_chIso_dZ2_dT_simVtx.push_back(chIso_dZ2_dT_simVtx[iRes]);

            evInfo[iRes].electron_chIso_dZ05.push_back(chIso_dZ05);
            evInfo[iRes].electron_chIso_dZ05_dT.push_back(chIso_dZ05_dT[iRes]);

            evInfo[iRes].electron_chIso_dZ1.push_back(chIso_dZ1);
            evInfo[iRes].electron_chIso_dZ1_dT.push_back(chIso_dZ1_dT[iRes]);

            evInfo[iRes].electron_chIso_dZ2.push_back(chIso_dZ2);
            evInfo[iRes].electron_chIso_dZ2_dT.push_back(chIso_dZ2_dT[iRes]);

            evInfo[iRes].electron_chIso_reldZ.push_back(chIso_reldZ);
            evInfo[iRes].electron_chIso_reldZ_dT.push_back(chIso_reldZ_dT[iRes]);

            evInfo[iRes].electron_chIso_dZele05.push_back(chIso_dZele05);
            evInfo[iRes].electron_chIso_dZele05_dTele.push_back(chIso_dZele05_dTele[iRes]);

            evInfo[iRes].electron_chIso_dZele1.push_back(chIso_dZele1);
            evInfo[iRes].electron_chIso_dZele1_dTele.push_back(chIso_dZele1_dTele[iRes]);

            evInfo[iRes].electron_chIso_dZele2.push_back(chIso_dZele2);
            evInfo[iRes].electron_chIso_dZele2_dTele.push_back(chIso_dZele2_dTele[iRes]);
        }  // end loop over time resolutions
    }      // end loop over barrel electrons

    electronIndex = 0;
    // -- start loop over endcap electrons
    for (unsigned int iele = 0; iele < endcapElectrons.size(); iele++) {

        const reco::GsfElectron& electron = endcapElectrons[iele];

        // -- minimal checks
        if (electron.pt() < 10.)
            continue;
        if (electron.gsfTrack().isNull())
            continue;
        //if (fabs(electron.eta()) > 1.5)
        //    continue;
        electronIndex++;
        // -- check if prompt or fake electron
        bool isPromptEle    = isPromptElectron(electron, genParticles);
        bool isMatchedJet   = isMatchedToGenJet(electron, genJets);
        bool isMatchedJet2  = isMatchedToGenJet2(electron, genJets);
        bool isFromTauDecay = isFromTau(electron, genParticles);
        // -- Look up and save the ID decisions
        /*float sieieCut          = 0.0112;
        float dEtaSeedCut       = 0.00377;
        float dPhiInCut         = 0.0884;
        float hoeCut            = Get_hoeCut();
        float relIsoWithEACut   = 0.133;
        float epCut             = 0.129;
        int   mHitsCut          = 1;
        bool  conversionVetoCut = TRUE;

        float          sieie        = electron.full5x5_sigmaIetaIeta();
        float          dEtaSeed     = Get_dEtaSeed(electron);
        float          dPhiIn       = electron.deltaPhiSuperClusterTrackAtVtx();
        float          hoe          = Get_hoe(electron);
        float          relIsoWithEA = 0.112 + 0.506 / electron.pt();
        float          ep           = Get_epCut(electron);
        int            mHits        = Get_mHitsCut(electron);
        bool           conversion float if (< 0.0105 &&);*/
        /*auto el            = barrelElectrons.ptrAt(iele);
        bool isPassVeto_   = (*veto_id_decisions)[el];
        bool isPassLoose_  = (*loose_id_decisions)[el];
        bool isPassMedium_ = (*medium_id_decisions)[el];
        bool isPassTight_  = (*tight_id_decisions)[el];*/

        // -- compute charged isolations
        float chIso_dZ05_simVtx = 0.;
        float chIso_dZ1_simVtx  = 0.;
        float chIso_dZ2_simVtx  = 0.;
        float chIso_dZ05        = 0.;
        float chIso_dZ1         = 0.;
        float chIso_dZ2         = 0.;
        float chIso_reldZ       = 0.;
        float chIso_dZele05     = 0.;
        float chIso_dZele1      = 0.;
        float chIso_dZele2      = 0.;

        float chIso_dZ05_dT_simVtx[nResol];
        float chIso_dZ1_dT_simVtx[nResol];
        float chIso_dZ2_dT_simVtx[nResol];
        float chIso_dZ05_dT[nResol];
        float chIso_dZ1_dT[nResol];
        float chIso_dZ2_dT[nResol];
        float chIso_reldZ_dT[nResol];
        float chIso_dZele05_dTele[nResol];
        float chIso_dZele1_dTele[nResol];
        float chIso_dZele2_dTele[nResol];
        float time[nResol];

        // -- initialize
        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
            chIso_dZ05_dT_simVtx[iRes] = 0.;
            chIso_dZ1_dT_simVtx[iRes]  = 0.;
            chIso_dZ2_dT_simVtx[iRes]  = 0.;
            chIso_dZ05_dT[iRes]        = 0.;
            chIso_dZ1_dT[iRes]         = 0.;
            chIso_dZ2_dT[iRes]         = 0.;
            chIso_reldZ_dT[iRes]       = 0.;
            chIso_dZele05_dTele[iRes]  = 0.;
            chIso_dZele1_dTele[iRes]   = 0.;
            chIso_dZele2_dTele[iRes]   = 0.;
            time[iRes]                 = 0.;
        }

        // -- firstly, get the electron time
        float        elecTime = 0.;  // record time of electron
        unsigned int icandTag = -1;
        float        dr_      = 99999999.;
        for (unsigned icand = 0; icand < pfcands.size(); ++icand) {
            const reco::PFCandidate& pfcand = pfcands[icand];
            if (pfcand.charge() == 0)
                continue;

            // -- get the track ref
            auto           pfcandRef = pfcands.refAt(icand);
            reco::TrackRef trackRef  = pfcandRef->trackRef();
            if (trackRef.isNull())
                continue;
            if (!trackRef->quality(reco::TrackBase::highPurity))
                continue;
            float dr = deltaR(electron.eta(), electron.phi(), pfcand.eta(), pfcand.phi());
            if (dr < dr_) {
                dr_      = dr;
                icandTag = icand;
                //-- get electron time from pfcand
                elecTime = pfcand.time();
            }
            if (isoConeDR_ == 0.3 && saveTracks_) {
                for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                    evInfo[iRes].drep.push_back(dr);
                }
            }
        }
        // -- firstly, get the electron time

        // -- compute loose electron ID
        //
        // -- full5x5_sigmaIetaIeta
        double sigmaIEtaIEta = electron.full5x5_sigmaIetaIeta();
        // -- abs(dEtaSeed)
        double dEtaInSeed = fabs(electron.superCluster().isNonnull() && electron.superCluster()->seed().isNonnull() ? electron.deltaEtaSuperClusterTrackAtVtx() - electron.superCluster()->eta() + electron.superCluster()->seed()->eta() : std::numeric_limits<float>::max());
        // -- abs(dPhiIn)
        double dPhiIn = fabs(electron.deltaPhiSuperClusterTrackAtVtx());
        // -- H/E
        double hoe       = electron.hadronicOverEm();
        double energy_sc = electron.superCluster()->energy();
        // -- relIsoWithEA
        float eA                = effectiveAreas_.getEffectiveArea(fabs(electron.superCluster()->eta()));
        float pf_isolation      = ((electron.pfIsolationVariables().sumChargedHadronPt + std::max(0., electron.pfIsolationVariables().sumNeutralHadronEt + electron.pfIsolationVariables().sumPhotonEt - eA * rho_)) / electron.pt());
        float pf_isolation_calo = ((electron.pfIsolationVariables().sumChargedHadronPt + std::max(0., electron.pfIsolationVariables().sumNeutralHadronEt + electron.pfIsolationVariables().sumPhotonEt - eA * rho_calo_)) / electron.pt());
        // --  EleEInverseMinusPInverseCut
        double ooEmooP = fabs((1.0 / electron.ecalEnergy() - electron.eSuperClusterOverP() / electron.ecalEnergy()));
        // -- expected missing inner hits
        double mHits = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
        //pass conversion veto
        int pass_conversion_veto = electron.passConversionVeto();

        // -- loop over charged pf candidates
        for (unsigned icand = 0; icand < pfcands.size(); ++icand) {
            const reco::PFCandidate& pfcand = pfcands[icand];
            if (pfcand.charge() == 0)
                continue;
            // -- get the track ref
            auto           pfcandRef = pfcands.refAt(icand);
            reco::TrackRef trackRef  = pfcandRef->trackRef();
            if (trackRef.isNull())
                continue;
            if (!trackRef->quality(reco::TrackBase::highPurity))
                continue;
            //-- use tracks with pT above thtreshold
            if (pfcand.pt() < minTrackPt_)
                continue;
            // -- get dz, dxy
            float dz3D  = std::abs(trackRef->dz(vtx3D.position()));
            float dz4D  = std::abs(trackRef->dz(vtx4D.position()));
            float dxy3D = std::abs(trackRef->dxy(vtx3D.position()));
            float dxy4D = std::abs(trackRef->dxy(vtx4D.position()));

            float dz3Drel = std::abs(dz3D / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx3D.zError() * vtx3D.zError()));
            float dz4Drel = std::abs(dz4D / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx4D.zError() * vtx4D.zError()));

            float dzsim  = std::abs(trackRef->vz() - genPV.position().z());
            float dxysim = sqrt(pow(trackRef->vx() - genPV.position().x(), 2) + pow(trackRef->vy() - genPV.position().y(), 2));

            float dzele  = std::abs(trackRef->dz(vtx4D.position()) - electron.gsfTrack()->dz(vtx4D.position()));
            float dxyele = std::abs(trackRef->dxy(vtx4D.position()) - electron.gsfTrack()->dxy(vtx4D.position()));

            float dr = deltaR(electron.eta(), electron.phi(), pfcand.eta(), pfcand.phi());

            // --- no timing

            // --- no timing
            if (dr > minDr_ && dr < isoConeDR_) {
                // -- sim vertex
                if (dzsim < 0.05 && dxysim < 0.02) {
                    chIso_dZ05_simVtx += pfcand.pt();
                }
                if (dzsim < 0.1 && dxysim < 0.02) {
                    chIso_dZ1_simVtx += pfcand.pt();
                }
                if (dzsim < 0.2 && dxysim < 0.02) {
                    chIso_dZ2_simVtx += pfcand.pt();
                }

                // -- reco vtx closest to the sim one
                if (dz3D < 0.05 && dxy3D < 0.02) {
                    chIso_dZ05 += pfcand.pt();
                }
                if (dz3D < 0.1 && dxy3D < 0.02) {
                    chIso_dZ1 += pfcand.pt();
                }
                if (dz3D < 0.2 && dxy3D < 0.02) {
                    chIso_dZ2 += pfcand.pt();
                }

                // -- using reco vtx closest to the sim one and cut on relative dz
                if (dz3Drel < 3.0 && dxy3D < 0.02) {
                    chIso_reldZ += pfcand.pt();
                }

                // -- dz wrt to electron
                if (dzele < 0.05) {
                    chIso_dZele05 += pfcand.pt();
                }
                if (dzele < 0.1) {
                    chIso_dZele1 += pfcand.pt();
                }
                if (dzele < 0.2) {
                    chIso_dZele2 += pfcand.pt();
                }
            }

            // --- with timing
            if (dr > minDr_ && dr < isoConeDR_) {
                for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                    double targetTimeResol  = timeResolutions_[iRes];
                    double defaultTimeResol = 0.;
                    if (pfcand.isTimeValid())
                        defaultTimeResol = double(pfcand.timeError());
                    double extra_resol = sqrt(targetTimeResol * targetTimeResol - defaultTimeResol * defaultTimeResol);
                    double dtsim       = 0.;
                    double dt          = 0.;
                    double dtele       = 0.;
                    time[iRes]         = -999.;
                    if (pfcand.isTimeValid()) {
                        double rnd = gRandom->Gaus(0., extra_resol);
                        //cout << "target time resol = "<< targetTimeResol << "  extra_resol = "<< extra_resol << "  rnd = " << rnd <<endl;
                        time[iRes] = pfcand.time() + rnd;
                        dtsim      = std::abs(time[iRes] - genPV.position().t() * 1000000000.);
                        cout << "dasim:" << dtsim << endl;
                        dt    = std::abs(time[iRes] - vtx4D.t());
                        dtele = std::abs(time[iRes] - elecTime);
                    }
                    else {
                        dtsim = 0.;
                        dt    = 0.;
                        dtele = 0.;
                    }

                    // -- sim vertex
                    if (dtsim < 3 * targetTimeResol && dzsim < 0.05 && dxysim < 0.02) {
                        chIso_dZ05_dT_simVtx[iRes] += pfcand.pt();
                    }
                    if (dtsim < 3 * targetTimeResol && dzsim < 0.1 && dxysim < 0.02) {
                        chIso_dZ1_dT_simVtx[iRes] += pfcand.pt();
                    }
                    if (dtsim < 3 * targetTimeResol && dzsim < 0.2 && dxysim < 0.02) {
                        chIso_dZ2_dT_simVtx[iRes] += pfcand.pt();
                    }

                    // -- reco vtx closest to the sim vtx
                    if (dt < 3 * targetTimeResol && dz4D < 0.05 && dxy4D < 0.02) {
                        chIso_dZ05_dT[iRes] += pfcand.pt();
                    }
                    if (dt < 3 * targetTimeResol && dz4D < 0.1 && dxy4D < 0.02) {
                        chIso_dZ1_dT[iRes] += pfcand.pt();
                    }
                    if (dt < 3 * targetTimeResol && dz4D < 0.2 && dxy4D < 0.02) {
                        chIso_dZ2_dT[iRes] += pfcand.pt();
                    }

                    // -- using reco vtx closest to the sim vtx and cut on relative dz
                    if (dt < 3 * targetTimeResol && dz4Drel < 3.0 && dxy4D < 0.02) {
                        chIso_reldZ_dT[iRes] += pfcand.pt();
                    }

                    // -- dT wrt to electron
                    if (dtele < sqrt(2) * 3 * targetTimeResol && dzele < 0.05) {
                        chIso_dZele05_dTele[iRes] += pfcand.pt();
                    }
                    if (dtele < sqrt(2) * 3 * targetTimeResol && dzele < 0.1) {
                        chIso_dZele1_dTele[iRes] += pfcand.pt();
                    }
                    if (dtele < sqrt(2) * 3 * targetTimeResol && dzele < 0.2) {
                        chIso_dZele2_dTele[iRes] += pfcand.pt();
                    }

                    // -- save info for tracks in the isolation cone (only for DR = 0.3)
                    if (saveTracks_ && (dz4D < 1.0 || dz3D < 1. || dzele < 1.)) {  // save a subset of tracks with loose dz selection
                        evInfo[iRes].track_t.push_back(time[iRes]);
                        evInfo[iRes].track_dz4D.push_back(trackRef->dz(vtx4D.position()));
                        evInfo[iRes].track_dz3D.push_back(trackRef->dz(vtx3D.position()));
                        evInfo[iRes].track_dxy4D.push_back(trackRef->dxy(vtx4D.position()));
                        evInfo[iRes].track_dxy3D.push_back(trackRef->dxy(vtx3D.position()));
                        evInfo[iRes].track_pt.push_back(pfcand.pt());
                        evInfo[iRes].track_eta.push_back(pfcand.eta());
                        evInfo[iRes].track_phi.push_back(pfcand.phi());
                        evInfo[iRes].track_elecIndex.push_back(electronIndex);
                    }
                }  // end loop over time resolutions
            }
        }  // end loop over tracks

        // -- fill electron info for each resolution scenario
        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {

            evInfo[iRes].electron_pt.push_back(electron.pt());
            evInfo[iRes].electron_eta.push_back(electron.eta());
            evInfo[iRes].electronsc_eta.push_back(electron.superCluster()->eta());
            evInfo[iRes].electron_phi.push_back(electron.phi());
            // -- info of electron work point
            evInfo[iRes].electron_sigmaIetaIeta.push_back(sigmaIEtaIEta);
            evInfo[iRes].electron_dEtaInSeed.push_back(dEtaInSeed);
            evInfo[iRes].electron_dPhiIn.push_back(dPhiIn);
            evInfo[iRes].electron_hoe.push_back(hoe);
            evInfo[iRes].electron_energy_sc.push_back(energy_sc);
            evInfo[iRes].electron_pf_isolation.push_back(pf_isolation);
            evInfo[iRes].electron_pf_isolation_calo.push_back(pf_isolation_calo);
            evInfo[iRes].electron_ooEmooP.push_back(ooEmooP);
            evInfo[iRes].electron_mHits.push_back(mHits);
            evInfo[iRes].electron_pass_conversion_veto.push_back(pass_conversion_veto);

            evInfo[iRes].electron_dz4D.push_back(electron.gsfTrack()->dz(vtx4D.position()));
            evInfo[iRes].electron_dxy4D.push_back(electron.gsfTrack()->dxy(vtx4D.position()));
            evInfo[iRes].electron_dz3D.push_back(electron.gsfTrack()->dz(vtx3D.position()));
            evInfo[iRes].electron_dxy3D.push_back(electron.gsfTrack()->dxy(vtx3D.position()));
            evInfo[iRes].electron_t.push_back(elecTime);
            evInfo[iRes].electron_isPrompt.push_back(isPromptEle);
            evInfo[iRes].electron_isMatchedToGenJet.push_back(isMatchedJet);
            evInfo[iRes].electron_isMatchedToGenJet2.push_back(isMatchedJet2);
            evInfo[iRes].electron_isFromTauDecay.push_back(isFromTauDecay);
            evInfo[iRes].electron_r9.push_back(electron.r9());
            evInfo[iRes].electron_chIso_dZ05_simVtx.push_back(chIso_dZ05_simVtx);
            evInfo[iRes].electron_chIso_dZ05_dT_simVtx.push_back(chIso_dZ05_dT_simVtx[iRes]);

            evInfo[iRes].electron_chIso_dZ1_simVtx.push_back(chIso_dZ1_simVtx);
            evInfo[iRes].electron_chIso_dZ1_dT_simVtx.push_back(chIso_dZ1_dT_simVtx[iRes]);

            evInfo[iRes].electron_chIso_dZ2_simVtx.push_back(chIso_dZ2_simVtx);
            evInfo[iRes].electron_chIso_dZ2_dT_simVtx.push_back(chIso_dZ2_dT_simVtx[iRes]);

            evInfo[iRes].electron_chIso_dZ05.push_back(chIso_dZ05);
            evInfo[iRes].electron_chIso_dZ05_dT.push_back(chIso_dZ05_dT[iRes]);

            evInfo[iRes].electron_chIso_dZ1.push_back(chIso_dZ1);
            evInfo[iRes].electron_chIso_dZ1_dT.push_back(chIso_dZ1_dT[iRes]);

            evInfo[iRes].electron_chIso_dZ2.push_back(chIso_dZ2);
            evInfo[iRes].electron_chIso_dZ2_dT.push_back(chIso_dZ2_dT[iRes]);

            evInfo[iRes].electron_chIso_reldZ.push_back(chIso_reldZ);
            evInfo[iRes].electron_chIso_reldZ_dT.push_back(chIso_reldZ_dT[iRes]);

            evInfo[iRes].electron_chIso_dZele05.push_back(chIso_dZele05);
            evInfo[iRes].electron_chIso_dZele05_dTele.push_back(chIso_dZele05_dTele[iRes]);

            evInfo[iRes].electron_chIso_dZele1.push_back(chIso_dZele1);
            evInfo[iRes].electron_chIso_dZele1_dTele.push_back(chIso_dZele1_dTele[iRes]);

            evInfo[iRes].electron_chIso_dZele2.push_back(chIso_dZele2);
            evInfo[iRes].electron_chIso_dZele2_dTele.push_back(chIso_dZele2_dTele[iRes]);
        }  // end loop over time resolutions
    }      // end loop over barrel electrons
    // -- fill info per event
    for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
        evInfo[iRes].npu          = nPU;
        evInfo[iRes].vtx4D_t      = vtx4D.t();
        evInfo[iRes].vtx4D_tErr   = vtx4D.tError();
        evInfo[iRes].vtx4D_z      = vtx4D.z();
        evInfo[iRes].vtx4D_zErr   = vtx4D.zError();
        evInfo[iRes].vtx3D_z      = vtx3D.z();
        evInfo[iRes].vtx3D_zErr   = vtx3D.zError();
        evInfo[iRes].vtxGen_z     = genPV.position().z();
        evInfo[iRes].vtxGen_t     = genPV.position().t();
        evInfo[iRes].vtx4D_isFake = vtx4D.isFake();
        evInfo[iRes].vtx3D_isFake = vtx3D.isFake();
        evInfo[iRes].rho          = rho_;
        evInfo[iRes].rho_calo     = rho_calo_;
    }

    // --- fill the tree
    for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
        eventTree[iRes]->Fill();
    }

    /*
    for (unsigned int iele = 0; iele < endcapElectrons.size(); iele++) {

        const reco::GsfElectron& electron = endcapElectrons[iele];

        // -- minimal checks
        if (electron.pt() < 15.)
            continue;
        if (electron.gsfTrack().isNull())
            continue;
        if (fabs(electron.eta()) > 1.5)
            continue;
        electronIndex++;
        // -- check if prompt or fake electron
        bool isPromptEle   = isPromptElectron(electron, genParticles);
        bool isMatchedJet  = isMatchedToGenJet(electron, genJets);
        bool isMatchedJet2 = isMatchedToGenJet2(electron, genJets);

        // -- Look up and save the ID decisions
        auto el            = endcapElectrons.ptrAt(iele);
        bool isPassVeto_   = (*veto_id_decisions)[el];
        bool isPassLoose_  = (*loose_id_decisions)[el];
        bool isPassMedium_ = (*medium_id_decisions)[el];
        bool isPassTight_  = (*tight_id_decisions)[el];

        // -- compute charged isolations
        float chIso[nCones];
        float chIso_dT[nCones][nResol];
        float chIso_dT_4D[nCones][nResol];
        float chIso_reldZ[nCones];
        float chIso_reldZ_dT[nCones][nResol];
        float chIso_reldZ_dT_4D[nCones][nResol];
        float chIso_simVtx[nCones];
        float chIso_dT_simVtx[nCones][nResol];
        float time[nCones][nResol];

        // -- initialize
        for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
            chIso_simVtx[iCone] = 0.;
            chIso[iCone]        = 0.;
            chIso_reldZ[iCone]  = 0.;
            for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                chIso_dT_simVtx[iCone][iRes] = 0.;
                chIso_dT[iCone][iRes]        = 0.;
                chIso_reldZ_dT[iCone][iRes]  = 0.;
                time[iCone][iRes]            = 0.;
            }
        }

        // -- firstly, get the electron time
        float        elecTime = 0.;  // record time of electron
        unsigned int icandTag = -1;
        float        dr_      = 99999999.;
        for (unsigned icand = 0; icand < pfcands.size(); ++icand) {
            const reco::PFCandidate& pfcand = pfcands[icand];
            if (pfcand.charge() == 0)
                continue;

            // -- get the track ref
            auto           pfcandRef = pfcands.refAt(icand);
            reco::TrackRef trackRef  = pfcandRef->trackRef();
            if (trackRef.isNull())
                continue;
            if (!trackRef->quality(reco::TrackBase::highPurity))
                continue;
            float dr = deltaR(electron.eta(), electron.phi(), pfcand.eta(), pfcand.phi());
            if (dr < dr_) {
                dr_      = dr;
                icandTag = icand;
                //-- get electron time from pfcand
                elecTime = pfcand.time();
            }
            for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                if (isoConeDR_[iCone] == 0.3 && saveTracks_) {
                    for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                        evInfo[iRes].drep.push_back(dr);
                    }
                }
            }
        }
        // -- firstly, get the electron time

        // -- loop over charged pf candidates
        for (unsigned icand = 0; icand < pfcands.size(); ++icand) {
            if (icand == icandTag)
                continue;
            const reco::PFCandidate& pfcand = pfcands[icand];
            if (pfcand.charge() == 0)
                continue;
            // -- get the track ref
            auto           pfcandRef = pfcands.refAt(icand);
            reco::TrackRef trackRef  = pfcandRef->trackRef();
            if (trackRef.isNull())
                continue;
            if (!trackRef->quality(reco::TrackBase::highPurity))
                continue;
            // -- get dz, dxy
            float dz    = std::abs(trackRef->dz(vtx.position()));
            float dz3D  = std::abs(trackRef->dz(vtx3D.position()));
            float dz4D  = std::abs(trackRef->dz(vtx4D.position()));
            float dxy   = std::abs(trackRef->dxy(vtx.position()));
            float dxy3D = std::abs(trackRef->dxy(vtx3D.position()));
            float dxy4D = std::abs(trackRef->dxy(vtx4D.position()));

            float dz3Drel = std::abs(dz3D / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx3D.zError() * vtx3D.zError()));
            float dzrel   = std::abs(dz / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx.zError() * vtx.zError()));
            float dz4Drel = std::abs(dz / std::sqrt(trackRef->dzError() * trackRef->dzError() + vtx4D.zError() * vtx4D.zError()));

            float dzsim  = std::abs(trackRef->vz() - genPV.position().z());
            float dxysim = sqrt(pow(trackRef->vx() - genPV.position().x(), 2) + pow(trackRef->vy() - genPV.position().y(), 2));

            float dr = deltaR(electron.eta(), electron.phi(), pfcand.eta(), pfcand.phi());

            // --- no timing

            // -- using sim vertex
            if (dzsim < maxDz_ && dxysim < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        chIso_simVtx[iCone] += pfcand.pt();
                    }
                }
            }

            // -- using reco vtx closest to the sim one
            if (dz3D < maxDz_ && dxy3D < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        chIso[iCone] += pfcand.pt();
                    }
                }
            }

            // -- using reco vtx closest to the sim one and cut on relative dz
            if (dz3Drel < 3.0 && dxy3D < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        chIso_reldZ[iCone] += pfcand.pt();
                    }
                }
            }

            // --- with timing

            // -- using sim vtx
            if (dzsim < maxDz_ && dxysim < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                            double time_resol  = timeResolutions_[iRes];
                            double extra_resol = sqrt(time_resol * time_resol - 0.03 * 0.03);
                            cout << iRes << "  " << time_resol << "  " << extra_resol << endl;
                            double dt         = 0.;
                            time[iCone][iRes] = -999.;
                            if (pfcand.isTimeValid()) {
                                time[iCone][iRes] = pfcand.time() + gRandom->Gaus(0., extra_resol);
                                dt                = std::abs(time[iCone][iRes] - genPV.position().t() * 1000000000.);
                                cout << pfcand.time() << "   " << time[iCone][iRes] << "   " << genPV.position().t() * 1000000000. << endl;
                            }
                            else {
                                dt = 0;
                            }
                            if (dt < 3 * time_resol) {
                                chIso_dT_simVtx[iCone][iRes] += pfcand.pt();
                            }
                        }  // end loop over time resolutions
                    }
                }  // end loop over cone sizes
            }

            // -- using reco vertex closest in dz to the gen vtx
            if (dz < maxDz_ && dxy < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                            double time_resol  = timeResolutions_[iRes];
                            double extra_resol = sqrt(time_resol * time_resol - 0.03 * 0.03);
                            double dt          = 0.;
                            if (pfcand.isTimeValid()) {
                                time[iCone][iRes] = pfcand.time() + gRandom->Gaus(0., extra_resol);
                                //elecTime          = elecTime + gRandom->Gaus(0., extra_resol);
                                //dt = std::abs(time[iCone][iRes] - vtx.t());
                                dt = std::abs(time[iCone][iRes] - vtx.t());
                                //dt                = std::abs(time[iCone][iRes] - elecTime);
                            }
                            else {
                                dt = 0;
                            }
                            //if (dt < 3 * sqrt(2) * time_resol) {
                            if (dt < 3 * time_resol) {
                                chIso_dT[iCone][iRes] += pfcand.pt();
                            }
                        }  // end loop over time resolutions

                        // save info for tracks in the isolation cone
                        if (isoConeDR_[iCone] == 0.3 && saveTracks_) {
                            for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                                evInfo[iRes].track_dz.push_back(trackRef->dz(vtx.position()));
                            }
                        }
                    }
                }  // end loop over cone sizes
            }
            // -- using reco vertex closest in dz to the gen vtx with timing info: vtx4D
            if (dz4D < maxDz_ && dxy4D < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                            double time_resol  = timeResolutions_[iRes];
                            double extra_resol = sqrt(time_resol * time_resol - 0.03 * 0.03);
                            double dt          = 0.;
                            if (pfcand.isTimeValid()) {
                                time[iCone][iRes] = pfcand.time() + gRandom->Gaus(0., extra_resol);
                                //elecTime          = elecTime + gRandom->Gaus(0., extra_resol);
                                //dt = std::abs(time[iCone][iRes] - vtx.t());
                                dt = std::abs(time[iCone][iRes] - vtx4D.t());
                                //dt                = std::abs(time[iCone][iRes] - elecTime);
                            }
                            else {
                                dt = 0;
                            }
                            //if (dt < 3 * sqrt(2) * time_resol) {
                            if (dt < 3 * time_resol) {
                                chIso_dT_4D[iCone][iRes] += pfcand.pt();
                            }
                        }  // end loop over time resolutions

                        // save info for tracks in the isolation cone
                        if (isoConeDR_[iCone] == 0.3 && saveTracks_) {
                            for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                                evInfo[iRes].track_t.push_back(time[iCone][iRes]);
                                evInfo[iRes].track_dz.push_back(trackRef->dz(vtx.position()));
                                evInfo[iRes].track_dz3D.push_back(trackRef->dz(vtx3D.position()));
                                evInfo[iRes].track_dz4D.push_back(trackRef->dz(vtx4D.position()));
                                evInfo[iRes].track_pt.push_back(pfcand.pt());
                                evInfo[iRes].track_eta.push_back(pfcand.eta());
                                evInfo[iRes].track_phi.push_back(pfcand.phi());
                                evInfo[iRes].track_elecIndex.push_back(electronIndex);
                            }
                        }
                    }
                }  // end loop over cone sizes
            }

            // -- using reco vertex closest in dz to the sim vertex and use a relative dz cut
            if (dzrel < 3.0 && dxy < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                            double time_resol  = timeResolutions_[iRes];
                            double extra_resol = sqrt(time_resol * time_resol - 0.03 * 0.03);
                            double dt          = 0.;
                            time[iCone][iRes]  = -999.;
                            if (pfcand.isTimeValid()) {
                                time[iCone][iRes] = pfcand.time() + gRandom->Gaus(0., extra_resol);
                                dt                = std::abs(time[iCone][iRes] - vtx.t());
                            }
                            else {
                                dt = 0;
                            }
                            if (dt < 3 * time_resol) {
                                chIso_reldZ_dT[iCone][iRes] += pfcand.pt();
                            }
                        }  // end loop over time resolutions
                    }
                }  // end loop over cone sizes
            }

            // -- using reco vertex closest in dz to the sim vertex and use a relative dz cut
            if (dz4Drel < 3.0 && dxy4D < 0.02) {
                for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                    if (dr > minDr_ && dr < isoConeDR_[iCone]) {
                        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
                            double time_resol  = timeResolutions_[iRes];
                            double extra_resol = sqrt(time_resol * time_resol - 0.03 * 0.03);
                            double dt          = 0.;
                            time[iCone][iRes]  = -999.;
                            if (pfcand.isTimeValid()) {
                                time[iCone][iRes] = pfcand.time() + gRandom->Gaus(0., extra_resol);
                                dt                = std::abs(time[iCone][iRes] - vtx4D.t());
                            }
                            else {
                                dt = 0;
                            }
                            if (dt < 3 * time_resol) {
                                chIso_reldZ_dT_4D[iCone][iRes] += pfcand.pt();
                            }
                        }  // end loop over time resolutions
                    }
                }  // end loop over cone sizes
            }

        }  // end loop over tracks

        // fill electron info
        for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
            evInfo[iRes].electron_pt.push_back(electron.pt());
            evInfo[iRes].electron_eta.push_back(electron.eta());
            evInfo[iRes].electron_phi.push_back(electron.phi());
            evInfo[iRes].electron_dz.push_back(electron.gsfTrack()->dz(vtx.position()));
            evInfo[iRes].electron_dxy.push_back(electron.gsfTrack()->dxy(vtx.position()));
            evInfo[iRes].electron_dz3D.push_back(electron.gsfTrack()->dz(vtx3D.position()));
            evInfo[iRes].electron_dxy3D.push_back(electron.gsfTrack()->dxy(vtx3D.position()));
            evInfo[iRes].electron_dz4D.push_back(electron.gsfTrack()->dz(vtx4D.position()));
            evInfo[iRes].electron_dxy4D.push_back(electron.gsfTrack()->dxy(vtx4D.position()));
            evInfo[iRes].electron_t.push_back(elecTime);
            evInfo[iRes].electron_isPrompt.push_back(isPromptEle);
            evInfo[iRes].electron_isMatchedToGenJet.push_back(isMatchedJet);
            evInfo[iRes].electron_isMatchedToGenJet2.push_back(isMatchedJet2);
            evInfo[iRes].electron_r9.push_back(electron.r9());
            evInfo[iRes].electron_sigmaIetaIeta.push_back(electron.sigmaIetaIeta());
            evInfo[iRes].passVetoId.push_back((int)isPassVeto_);
            evInfo[iRes].passLooseId.push_back((int)isPassLoose_);
            evInfo[iRes].passMediumId.push_back((int)isPassMedium_);
            evInfo[iRes].passTightId.push_back((int)isPassTight_);

            for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
                (evInfo[iRes].electron_chIso_simVtx[iCone]).push_back(chIso_simVtx[iCone]);
                (evInfo[iRes].electron_chIso_dT_simVtx[iCone][iRes]).push_back(chIso_dT_simVtx[iCone][iRes]);
                (evInfo[iRes].electron_chIso[iCone]).push_back(chIso[iCone]);
                (evInfo[iRes].electron_chIso_dT[iCone][iRes]).push_back(chIso_dT[iCone][iRes]);
                (evInfo[iRes].electron_chIso_dT_4D[iCone][iRes]).push_back(chIso_dT_4D[iCone][iRes]);
                (evInfo[iRes].electron_chIso_reldZ[iCone]).push_back(chIso_reldZ[iCone]);
                (evInfo[iRes].electron_chIso_reldZ_dT[iCone][iRes]).push_back(chIso_reldZ_dT[iCone][iRes]);
                (evInfo[iRes].electron_chIso_reldZ_dT_4D[iCone][iRes]).push_back(chIso_reldZ_dT_4D[iCone][iRes]);
            }
        }
    }  // end loop over endcap electrons
   */
}  // end analyze event

// ------------ method called once each job just before starting event loop  ------------
void ElectronIsolationAnalyzer::beginJob() {
    for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
        eventTree[iRes]->Branch("drep", &evInfo[iRes].drep);  // add for veto dr

        eventTree[iRes]->Branch("npu", &evInfo[iRes].npu);
        eventTree[iRes]->Branch("vtxGen_z", &evInfo[iRes].vtxGen_z);
        eventTree[iRes]->Branch("vtxGen_t", &evInfo[iRes].vtxGen_t);
        eventTree[iRes]->Branch("vtx3D_z", &evInfo[iRes].vtx3D_z);
        eventTree[iRes]->Branch("vtx3D_zErr", &evInfo[iRes].vtx3D_zErr);
        eventTree[iRes]->Branch("vtx4D_z", &evInfo[iRes].vtx4D_z);
        eventTree[iRes]->Branch("vtx4D_zErr", &evInfo[iRes].vtx4D_zErr);
        eventTree[iRes]->Branch("vtx4D_t", &evInfo[iRes].vtx4D_t);
        eventTree[iRes]->Branch("vtx4D_tErr", &evInfo[iRes].vtx4D_tErr);
        eventTree[iRes]->Branch("vtx4D_isFake", &evInfo[iRes].vtx4D_isFake);
        eventTree[iRes]->Branch("vtx3D_isFake", &evInfo[iRes].vtx3D_isFake);
        eventTree[iRes]->Branch("electron_pt", &evInfo[iRes].electron_pt);
        eventTree[iRes]->Branch("electron_eta", &evInfo[iRes].electron_eta);
        eventTree[iRes]->Branch("electronsc_eta", &evInfo[iRes].electronsc_eta);
        eventTree[iRes]->Branch("electron_phi", &evInfo[iRes].electron_phi);
        eventTree[iRes]->Branch("electron_dz3D", &evInfo[iRes].electron_dz3D);
        eventTree[iRes]->Branch("electron_dz4D", &evInfo[iRes].electron_dz4D);
        eventTree[iRes]->Branch("electron_dxy3D", &evInfo[iRes].electron_dxy3D);
        eventTree[iRes]->Branch("electron_dxy4D", &evInfo[iRes].electron_dxy4D);
        eventTree[iRes]->Branch("electron_t", &evInfo[iRes].electron_t);
        eventTree[iRes]->Branch("electron_isPrompt", &evInfo[iRes].electron_isPrompt);
        eventTree[iRes]->Branch("electron_isMatchedToGenJet", &evInfo[iRes].electron_isMatchedToGenJet);
        eventTree[iRes]->Branch("electron_isFromTauDecay", &evInfo[iRes].electron_isFromTauDecay);

        // -- info of electron work point
        eventTree[iRes]->Branch("rho", &evInfo[iRes].rho);
        eventTree[iRes]->Branch("rho_calo", &evInfo[iRes].rho_calo);
        eventTree[iRes]->Branch("electron_sigmaIetaIeta", &evInfo[iRes].electron_sigmaIetaIeta);
        eventTree[iRes]->Branch("electron_dEtaInSeed", &evInfo[iRes].electron_dEtaInSeed);
        eventTree[iRes]->Branch("electron_dPhiIn", &evInfo[iRes].electron_dPhiIn);
        eventTree[iRes]->Branch("electron_hoe", &evInfo[iRes].electron_hoe);
        eventTree[iRes]->Branch("electron_energy_sc", &evInfo[iRes].electron_energy_sc);
        eventTree[iRes]->Branch("electron_pf_isolation", &evInfo[iRes].electron_pf_isolation);
        eventTree[iRes]->Branch("electron_pf_isolation_calo", &evInfo[iRes].electron_pf_isolation_calo);
        eventTree[iRes]->Branch("electron_ooEmooP", &evInfo[iRes].electron_ooEmooP);
        eventTree[iRes]->Branch("electron_mHits", &evInfo[iRes].electron_mHits);
        eventTree[iRes]->Branch("electron_pass_conversion_veto", &evInfo[iRes].electron_pass_conversion_veto);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ05_simVtx", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ05_simVtx);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ05_dT_simVtx", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ05_dT_simVtx);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ1_simVtx", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ1_simVtx);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ1_dT_simVtx", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ1_dT_simVtx);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ2_simVtx", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ2_simVtx);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ2_dT_simVtx", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ2_dT_simVtx);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ05", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ05);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ05_dT", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ05_dT);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ1", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ1);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ1_dT", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ1_dT);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ2", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ2);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZ2_dT", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZ2_dT);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_reldZ", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_reldZ);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_reldZ_dT", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_reldZ_dT);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZele05", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZele05);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZele05_dTele", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZele05_dTele);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZele1", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZele1);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZele1_dTele", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZele1_dTele);

        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZele2", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZele2);
        eventTree[iRes]->Branch(Form("electron_chIso%.2d_dZele2_dTele", int(isoConeDR_ * 10)), &evInfo[iRes].electron_chIso_dZele2_dTele);

        if (saveTracks_) {
            eventTree[iRes]->Branch("track_t", &evInfo[iRes].track_t);
            eventTree[iRes]->Branch("track_dz4D", &evInfo[iRes].track_dz4D);
            eventTree[iRes]->Branch("track_dz3D", &evInfo[iRes].track_dz3D);
            eventTree[iRes]->Branch("track_dxy4D", &evInfo[iRes].track_dxy4D);
            eventTree[iRes]->Branch("track_dxy3D", &evInfo[iRes].track_dxy3D);
            eventTree[iRes]->Branch("track_pt", &evInfo[iRes].track_pt);
            eventTree[iRes]->Branch("track_eta", &evInfo[iRes].track_eta);
            eventTree[iRes]->Branch("track_phi", &evInfo[iRes].track_phi);
            eventTree[iRes]->Branch("track_elecIndex", &evInfo[iRes].track_elecIndex);
        }
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void ElectronIsolationAnalyzer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ElectronIsolationAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// ------------ method initialize tree structure ----------------------------------------------------
void ElectronIsolationAnalyzer::initEventStructure() {
    // per-event trees:
    for (unsigned int iRes = 0; iRes < timeResolutions_.size(); iRes++) {
        evInfo[iRes].drep.clear();
        evInfo[iRes].npu          = -1;
        evInfo[iRes].vtxGen_t     = -999;
        evInfo[iRes].vtxGen_z     = -999;
        evInfo[iRes].vtx4D_z      = -999;
        evInfo[iRes].vtx4D_zErr   = -999;
        evInfo[iRes].vtx4D_t      = -999;
        evInfo[iRes].vtx4D_tErr   = -999;
        evInfo[iRes].vtx3D_z      = -999;
        evInfo[iRes].vtx3D_zErr   = -999;
        evInfo[iRes].vtx4D_isFake = -999;
        evInfo[iRes].vtx3D_isFake = -999;
        evInfo[iRes].rho          = -999;
        evInfo[iRes].rho_calo     = -999;

        evInfo[iRes].electron_pt.clear();
        evInfo[iRes].electron_eta.clear();
        evInfo[iRes].electronsc_eta.clear();
        evInfo[iRes].electron_phi.clear();

        // -- info of electron work point
        evInfo[iRes].electron_sigmaIetaIeta.clear();
        evInfo[iRes].electron_dEtaInSeed.clear();
        evInfo[iRes].electron_dPhiIn.clear();
        evInfo[iRes].electron_hoe.clear();
        evInfo[iRes].electron_energy_sc.clear();
        evInfo[iRes].electron_pf_isolation.clear();
        evInfo[iRes].electron_pf_isolation_calo.clear();
        evInfo[iRes].electron_ooEmooP.clear();
        evInfo[iRes].electron_mHits.clear();
        evInfo[iRes].electron_pass_conversion_veto.clear();

        evInfo[iRes].electron_dz3D.clear();
        evInfo[iRes].electron_dxy3D.clear();
        evInfo[iRes].electron_dz4D.clear();
        evInfo[iRes].electron_dxy4D.clear();
        evInfo[iRes].electron_t.clear();
        evInfo[iRes].electron_isPrompt.clear();
        evInfo[iRes].electron_isMatchedToGenJet.clear();
        evInfo[iRes].electron_isMatchedToGenJet2.clear();
        evInfo[iRes].electron_r9.clear();
        evInfo[iRes].electron_isFromTauDecay.clear();
        /*
        // clear IDs
        evInfo[iRes].passVetoId.clear();
        evInfo[iRes].passLooseId.clear();
        evInfo[iRes].passMediumId.clear();
        evInfo[iRes].passTightId.clear();

        for (unsigned int iCone = 0; iCone < isoConeDR_.size(); iCone++) {
            evInfo[iRes].electron_chIso_simVtx[iCone].clear();
            evInfo[iRes].electron_chIso_dT_simVtx[iCone][iRes].clear();
            evInfo[iRes].electron_chIso[iCone].clear();
            evInfo[iRes].electron_chIso_dT[iCone][iRes].clear();
            evInfo[iRes].electron_chIso_reldZ[iCone].clear();
            evInfo[iRes].electron_chIso_reldZ_dT[iCone][iRes].clear();
        }
*/
        evInfo[iRes].electron_chIso_dZ05_simVtx.clear();
        evInfo[iRes].electron_chIso_dZ05_dT_simVtx.clear();
        evInfo[iRes].electron_chIso_dZ1_simVtx.clear();
        evInfo[iRes].electron_chIso_dZ1_dT_simVtx.clear();
        evInfo[iRes].electron_chIso_dZ2_simVtx.clear();
        evInfo[iRes].electron_chIso_dZ2_dT_simVtx.clear();

        evInfo[iRes].electron_chIso_dZ05.clear();
        evInfo[iRes].electron_chIso_dZ05_dT.clear();
        evInfo[iRes].electron_chIso_dZ1.clear();
        evInfo[iRes].electron_chIso_dZ1_dT.clear();
        evInfo[iRes].electron_chIso_dZ2.clear();
        evInfo[iRes].electron_chIso_dZ2_dT.clear();

        evInfo[iRes].electron_chIso_reldZ.clear();
        evInfo[iRes].electron_chIso_reldZ_dT.clear();

        evInfo[iRes].electron_chIso_dZele05.clear();
        evInfo[iRes].electron_chIso_dZele05_dTele.clear();
        evInfo[iRes].electron_chIso_dZele1.clear();
        evInfo[iRes].electron_chIso_dZele1_dTele.clear();
        evInfo[iRes].electron_chIso_dZele2.clear();
        evInfo[iRes].electron_chIso_dZele2_dTele.clear();
        if (saveTracks_) {
            evInfo[iRes].track_t.clear();
            evInfo[iRes].track_dz4D.clear();
            evInfo[iRes].track_dz3D.clear();
            evInfo[iRes].track_dxy4D.clear();
            evInfo[iRes].track_dxy3D.clear();
            evInfo[iRes].track_pt.clear();
            evInfo[iRes].track_eta.clear();
            evInfo[iRes].track_phi.clear();
            evInfo[iRes].track_elecIndex.clear();
        }
    }
}

// --- Electrons matching to gen level -------- ----------------------------------------------------
bool isPromptElectron(const reco::GsfElectron& electron, const edm::View<reco::GenParticle>& genParticles) {
    bool isPrompt = false;

    for (unsigned int ip = 0; ip < genParticles.size(); ip++) {
        const reco::GenParticle& genp = genParticles[ip];
        if (std::abs(genp.pdgId()) != 11)
            continue;
        //if (genp.status() != 1 || !genp.isLastCopy())
        //  continue;  // -- from Simone
        if (!genp.isPromptFinalState())
            continue;
        if (genp.pt() < 5.0)
            continue;
        double dr = deltaR(electron, genp);
        if (dr > 0.1) {
            continue;
        }
        else {
            isPrompt = true;
            break;
        }
    }

    return isPrompt;
}

// --- matching to gen jet
bool isMatchedToGenJet(const reco::GsfElectron& electron, const edm::View<reco::GenJet>& genJets) {
    bool isMatched = false;

    for (unsigned int ip = 0; ip < genJets.size(); ip++) {
        const reco::GenJet& genj = genJets[ip];
        if (genj.pt() < 15.0 || genj.hadEnergy() / genj.energy() < 0.3)  // add information for hoe
            continue;
        double dr = deltaR(electron, genj);
        if (dr > 0.1) {
            continue;
        }
        else {
            isMatched = true;
            break;
        }
    }

    return isMatched;
}

bool isMatchedToGenJet2(const reco::GsfElectron& electron, const edm::View<reco::GenJet>& genJets) {
    bool isMatched = false;

    for (unsigned int ip = 0; ip < genJets.size(); ip++) {
        const reco::GenJet& genj = genJets[ip];
        if (genj.pt() < 15.0)
            continue;
        double dr = deltaR(electron, genj);
        if (dr > 0.1) {
            continue;
        }
        else {
            isMatched = true;
            break;
        }
    }

    return isMatched;
}

// --- matching to electrons from tau decays
bool isFromTau(const reco::GsfElectron& electron, const edm::View<reco::GenParticle>& genParticles) {

    bool fromTau = false;

    for (unsigned int ip = 0; ip < genParticles.size(); ip++) {
        const reco::GenParticle& genp = genParticles[ip];
        if (std::abs(genp.pdgId()) != 11)
            continue;
        if (!genp.isDirectPromptTauDecayProductFinalState())
            continue;
        double dr = deltaR(electron, genp);
        if (dr > 0.1) {
            continue;
        }
        else {
            fromTau = true;
            break;
        }
    }

    return fromTau;
}

/*float Get_dEtaInSeed(const reco::GsfElectron& ele) {
    return ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() ? ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

float Get_epCut(const reco::Candidate& ele) {
    const float ecal_energy_inverse = 1.0 / ele.ecalEnergy();
    const float eSCoverP            = ele.eSuperClusterOverP();
    return std::abs(1.0 - eSCoverP) * ecal_energy_inverse;
}

int Get_mHitsCut(const reco::Candidate& ele) {
    constexpr reco::HitPattern::HitCategory missingHitType = reco::HitPattern::MISSING_INNER_HITS;

    int mHits = electron.gsfTrack()->hitPattern().numberOfHits(missingHitType);
    return mHits;
}
float Get_hoe(const reco::Candidate& ele) {
    return 0.05 + 1.16 / ESC + 0.0324ρ / ele.superCluster()->energy();
}*/
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronIsolationAnalyzer);
