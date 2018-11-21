// per compilare: g++ -Wall -o EleIsoAnalysis `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFit -lRooFitCore -lFoam -lHtml -lMinuit eleiso.cc
/*
./EleIsoAnalysis /home/pku/xiaoj/eleiso/withTrack/DY PU200_DY 30 1 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/withTrack/DY-noPU noPU_DY 30 1 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/withTrack/QCD PU200_QCD 30 0 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/withTrack/QCD-noPU noPU_QCD 30 0 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/newplace/DY PU200_DY 30 1 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/newplace/DY_noPU noPU_DY 30 1 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/newplace/QCD PU200_QCD 30 0 0
./EleIsoAnalysis /home/pku/xiaoj/eleiso/newplace/QCD_noPU noPU_QCD 30 0 0
*/
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVirtualFitter.h"

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
    // Check the number of parameters
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " [folder] [samplename] [timeResol]  [usePromptelectrons] [ptReweighting]" << std::endl;
        return 1;
    }

    string  folderName(argv[1]);
    TString folder(argv[1]);
    TString samplename(argv[2]);
    TString timeResolution(argv[3]);
    int     usePromptelectrons = atoi(argv[4]);
    int     applyPtReweighting = atoi(argv[5]);

    TH2F* hratio2D;
    if (applyPtReweighting) {
        TFile* fWeights;
        if (folderName.find("noPU") != std::string::npos) {
            fWeights = TFile::Open("ptRatio_mu_noPU_TTbar.root");
            if (folderName.find("QCD") != std::string::npos)
                fWeights = TFile::Open("ptRatio_mu_noPU_QCD.root");
        }
        else {
            fWeights = TFile::Open("ptRatio_mu_TTbar.root");
            if (folderName.find("QCD") != std::string::npos)
                fWeights = TFile::Open("ptRatio_mu_QCD.root");
        }
        std::cout << "Applying pT/eta weights from " << fWeights->GetName() << std::endl;
        hratio2D = (TH2F*)fWeights->Get("hratio2D");
    }

    cout << endl;
    cout << "Start analyzing " << folderName << endl;

    // -- get TChain
    //TChain* chain = new TChain(("analysis/tree_" + timeResolution + "ps").c_str(), "tree");
    //chain->Add((folderName + "/eleiso_*.root").c_str());
    TFile* ipt = TFile::Open(folder + "/" + samplename);
    // -- tree vars
    TTreeReader fReader("analysis/tree_" + timeResolution + "ps", ipt);  //!the tree reader
    //TTree*      fChain = 0;  //!pointer to the analyzed TTree or TChain

    // Readers to access the data (delete the ones you do not need).
    TTreeReaderValue<Int_t>             npu                        = {fReader, "npu"};
    TTreeReaderValue<Float_t>           vtxGen_z                   = {fReader, "vtxGen_z"};
    TTreeReaderValue<Float_t>           vtxGen_t                   = {fReader, "vtxGen_t"};
    TTreeReaderValue<Float_t>           vtx3D_z                    = {fReader, "vtx3D_z"};
    TTreeReaderValue<Float_t>           vtx_z                      = {fReader, "vtx_z"};
    TTreeReaderValue<Float_t>           vtx_t                      = {fReader, "vtx_t"};
    TTreeReaderValue<Float_t>           vtx4D_z                    = {fReader, "vtx4D_z"};
    TTreeReaderValue<Float_t>           vtx4D_t                    = {fReader, "vtx4D_t"};
    TTreeReaderArray<float>             drep                       = {fReader, "drep"};
    TTreeReaderArray<float>             electron_pt                = {fReader, "electron_pt"};
    TTreeReaderArray<float>             electron_eta               = {fReader, "electron_eta"};
    TTreeReaderArray<float>             electron_phi               = {fReader, "electron_phi"};
    TTreeReaderArray<float>             electron_dz3D              = {fReader, "electron_dz3D"};
    TTreeReaderArray<float>             electron_dz4D              = {fReader, "electron_dz4D"};
    TTreeReaderArray<float>             electron_dxy3D             = {fReader, "electron_dxy3D"};
    TTreeReaderArray<float>             electron_dxy4D             = {fReader, "electron_dxy4D"};
    TTreeReaderValue<std::vector<bool>> electron_isPrompt          = {fReader, "electron_isPrompt"};
    TTreeReaderValue<std::vector<bool>> electron_isMatchedToGenJet = {fReader, "electron_isMatchedToGenJet"};
    TTreeReaderArray<float>             electron_r9                = {fReader, "electron_r9"};
    TTreeReaderArray<float>             electron_sigmaIetaIeta     = {fReader, "electron_sigmaIetaIeta"};
    TTreeReaderArray<float>             electron_chIso02           = {fReader, "electron_chIso02"};
    TTreeReaderArray<float>             electron_chIso02_dT        = {fReader, "electron_chIso02_dT"};
    TTreeReaderArray<float>             electron_chIso03           = {fReader, "electron_chIso03"};
    TTreeReaderArray<float>             electron_chIso03_dT        = {fReader, "electron_chIso03_dT"};
    TTreeReaderArray<float>             electron_chIso04           = {fReader, "electron_chIso04"};
    TTreeReaderArray<float>             electron_chIso04_dT        = {fReader, "electron_chIso04_dT"};
    TTreeReaderArray<float>             electron_chIso05           = {fReader, "electron_chIso05"};
    TTreeReaderArray<float>             electron_chIso05_dT        = {fReader, "electron_chIso05_dT"};
    TTreeReaderArray<float>             track_t                    = {fReader, "track_t"};
    TTreeReaderArray<float>             track_dz                   = {fReader, "track_dz"};
    TTreeReaderArray<float>             track_dz3D                 = {fReader, "track_dz3D"};
    TTreeReaderArray<float>             track_pt                   = {fReader, "track_pt"};
    TTreeReaderArray<float>             track_eta                  = {fReader, "track_eta"};
    TTreeReaderArray<float>             track_phi                  = {fReader, "track_phi"};
    // -- book histograms
    TH1F* h_npu = new TH1F("h_npu", "h_npu", 200, 140, 260);
    h_npu->GetXaxis()->SetTitle("number of pileup vertices");

    TH1F* h_vtx_dz3D = new TH1F("h_vtx_dz3D", "h_vtx_dz3D", 1000, -0.2, 0.2);
    h_vtx_dz3D->GetXaxis()->SetTitle("z_{3Dvtx} - z_{gen} (mm)");

    TH1F* h_vtx_dz = new TH1F("h_vtx_dz", "h_vtx_dz", 1000, -0.2, 0.2);
    h_vtx_dz->GetXaxis()->SetTitle("z_{vtx} - z_{gen} (mm)");

    TH1F* h_vtx_dt = new TH1F("h_vtx_dt", "h_vtx_dt", 1000, -0.5, 0.5);
    h_vtx_dt->GetXaxis()->SetTitle("t_{vtx} - t_{gen} (ns)");

    TH1F* h_vtx_dz4D = new TH1F("h_vtx_dz4D", "h_vtx_dz4D", 1000, -0.2, 0.2);
    h_vtx_dz4D->GetXaxis()->SetTitle("z_{4Dvtx} - z_{gen} (mm)");

    TH1F* h_vtx_dt4D = new TH1F("h_vtx_dt4D", "h_vtx_dt4D", 1000, -0.5, 0.5);
    h_vtx_dt4D->GetXaxis()->SetTitle("t_{4Dvtx} - t_{gen} (ns)");

    TH1F* h_vtx_dz3D_pull = new TH1F("h_vtx_dz3D_pull", "h_vtx_dz3D_pull", 200, -10, 10);
    h_vtx_dz3D_pull->GetXaxis()->SetTitle("(z_{3Dvtx} - z_{gen})/#sigma_{z}");

    TH1F* h_vtx_dz4D_pull = new TH1F("h_vtx_dz4D_pull", "h_vtx_dz4D_pull", 200, -10, 10);
    h_vtx_dz4D_pull->GetXaxis()->SetTitle("(z_{4Dvtx} - z_{gen})/#sigma_{z}");

    TH1F* h_vtx_dt4D_pull = new TH1F("h_vtx_dt4D_pull", "h_vtx_dt4D_pull", 200, -10, 10);
    h_vtx_dt4D_pull->GetXaxis()->SetTitle("(t_{4Dvtx} - t_{gen})/#sigma_{t}");

    TH1F* h_vtx_dz_pull = new TH1F("h_vtx_dz_pull", "h_vtx_dz_pull", 200, -10, 10);
    h_vtx_dz_pull->GetXaxis()->SetTitle("(z_{vtx} - z_{gen})/#sigma_{z}");

    TH1F* h_vtx_dt_pull = new TH1F("h_vtx_dt_pull", "h_vtx_dt_pull", 200, -10, 10);
    h_vtx_dt_pull->GetXaxis()->SetTitle("(t_{vtx} - t_{gen})/#sigma_{t}");

    TH1F* h_electron_pt = new TH1F("h_electron_pt", "h_electron_pt", 200, 0, 200);
    h_electron_pt->GetXaxis()->SetTitle("electron p_{T} (GeV)");

    TH1F* h_electron_eta = new TH1F("h_electron_eta", "h_electron_eta", 100, -3, 3);
    h_electron_eta->GetXaxis()->SetTitle("electron eta");

    TH1F* h_electron_phi = new TH1F("h_electron_phi", "h_electron_phi", 100, -4, 4);
    h_electron_phi->GetXaxis()->SetTitle("electron phi");

    TH2F* h2_electron_pt_vs_eta = new TH2F("h2_electron_pt_vs_eta", "h2_electron_pt_vs_eta", 150, -3, 3, 500, 0, 1000);

    // dR distribution
    TH1F* h_drep = new TH1F("h_drep", "h_drep", 3000, 0, 0.3);
    h_drep->GetXaxis()->SetTitle("dR_{(e,p)}");

    // -- relative isolation
    TH1F* h_electron_relChIso02 = new TH1F("h_electron_relChIso02", "h_electron_relChIso02", 5000, 0, 5);
    h_electron_relChIso02->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT = new TH1F("h_electron_relChIso02_dT", "h_electron_relChIso02_dT", 5000, 0, 5);
    h_electron_relChIso02_dT->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT_4D = new TH1F("h_electron_relChIso02_dT_4D", "h_electron_relChIso02_dT_4D", 5000, 0, 5);
    h_electron_relChIso02_dT_4D->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03 = new TH1F("h_electron_relChIso03", "h_electron_relChIso03", 5000, 0, 5);
    h_electron_relChIso03->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT = new TH1F("h_electron_relChIso03_dT", "h_electron_relChIso03_dT", 5000, 0, 5);
    h_electron_relChIso03_dT->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT_4D = new TH1F("h_electron_relChIso03_dT_4D", "h_electron_relChIso03_dT_4D", 5000, 0, 5);
    h_electron_relChIso03_dT_4D->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04 = new TH1F("h_electron_relChIso04", "h_electron_relChIso04", 5000, 0, 5);
    h_electron_relChIso04->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT = new TH1F("h_electron_relChIso04_dT", "h_electron_relChIso04_dT", 5000, 0, 5);
    h_electron_relChIso04_dT->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT_4D = new TH1F("h_electron_relChIso04_dT_4D", "h_electron_relChIso04_dT_4D", 5000, 0, 5);
    h_electron_relChIso04_dT_4D->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05 = new TH1F("h_electron_relChIso05", "h_electron_relChIso05", 5000, 0, 5);
    h_electron_relChIso05->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT = new TH1F("h_electron_relChIso05_dT", "h_electron_relChIso05_dT", 5000, 0, 5);
    h_electron_relChIso05_dT->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT_4D = new TH1F("h_electron_relChIso05_dT_4D", "h_electron_relChIso05_dT_4D", 5000, 0, 5);
    h_electron_relChIso05_dT_4D->GetXaxis()->SetTitle("relative charged isolation");

    // -- relative isolation barrel
    TH1F* h_electron_relChIso02_barrel = new TH1F("h_electron_relChIso02_barrel", "h_electron_relChIso02_barrel", 5000, 0, 5);
    h_electron_relChIso02_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT_barrel = new TH1F("h_electron_relChIso02_dT_barrel", "h_electron_relChIso02_dT_barrel", 5000, 0, 5);
    h_electron_relChIso02_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_barrel = new TH1F("h_electron_relChIso03_barrel", "h_electron_relChIso03_barrel", 5000, 0, 5);
    h_electron_relChIso03_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT_barrel = new TH1F("h_electron_relChIso03_dT_barrel", "h_electron_relChIso03_dT_barrel", 5000, 0, 5);
    h_electron_relChIso03_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_barrel = new TH1F("h_electron_relChIso04_barrel", "h_electron_relChIso04_barrel", 5000, 0, 5);
    h_electron_relChIso04_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT_barrel = new TH1F("h_electron_relChIso04_dT_barrel", "h_electron_relChIso04_dT_barrel", 5000, 0, 5);
    h_electron_relChIso04_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_barrel = new TH1F("h_electron_relChIso05_barrel", "h_electron_relChIso05_barrel", 5000, 0, 5);
    h_electron_relChIso05_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT_barrel = new TH1F("h_electron_relChIso05_dT_barrel", "h_electron_relChIso05_dT_barrel", 5000, 0, 5);
    h_electron_relChIso05_dT_barrel->GetXaxis()->SetTitle("relative charged isolation");

    // -- relative isolation endcap
    TH1F* h_electron_relChIso02_endcap = new TH1F("h_electron_relChIso02_endcap", "h_electron_relChIso02_endcap", 5000, 0, 5);
    h_electron_relChIso02_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT_endcap = new TH1F("h_electron_relChIso02_dT_endcap", "h_electron_relChIso02_dT_endcap", 5000, 0, 5);
    h_electron_relChIso02_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_endcap = new TH1F("h_electron_relChIso03_endcap", "h_electron_relChIso03_endcap", 5000, 0, 5);
    h_electron_relChIso03_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT_endcap = new TH1F("h_electron_relChIso03_dT_endcap", "h_electron_relChIso03_dT_endcap", 5000, 0, 5);
    h_electron_relChIso03_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_endcap = new TH1F("h_electron_relChIso04_endcap", "h_electron_relChIso04_endcap", 5000, 0, 5);
    h_electron_relChIso04_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT_endcap = new TH1F("h_electron_relChIso04_dT_endcap", "h_electron_relChIso04_dT_endcap", 5000, 0, 5);
    h_electron_relChIso04_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_endcap = new TH1F("h_electron_relChIso05_endcap", "h_electron_relChIso05_endcap", 5000, 0, 5);
    h_electron_relChIso05_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT_endcap = new TH1F("h_electron_relChIso05_dT_endcap", "h_electron_relChIso05_dT_endcap", 5000, 0, 5);
    h_electron_relChIso05_dT_endcap->GetXaxis()->SetTitle("relative charged isolation");

    // chIso wrt to sim vertex
    TH1F* h_electron_relChIso02_simVtx = new TH1F("h_electron_relChIso02_simVtx", "h_electron_relChIso02_simVtx", 5000, 0, 5);
    h_electron_relChIso02_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT_simVtx = new TH1F("h_electron_relChIso02_dT_simVtx", "h_electron_relChIso02_dT_simVtx", 5000, 0, 5);
    h_electron_relChIso02_dT_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_simVtx = new TH1F("h_electron_relChIso03_simVtx", "h_electron_relChIso03_simVtx", 5000, 0, 5);
    h_electron_relChIso03_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT_simVtx = new TH1F("h_electron_relChIso03_dT_simVtx", "h_electron_relChIso03_dT_simVtx", 5000, 0, 5);
    h_electron_relChIso03_dT_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_simVtx = new TH1F("h_electron_relChIso04_simVtx", "h_electron_relChIso04_simVtx", 5000, 0, 5);
    h_electron_relChIso04_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT_simVtx = new TH1F("h_electron_relChIso04_dT_simVtx", "h_electron_relChIso04_dT_simVtx", 5000, 0, 5);
    h_electron_relChIso04_dT_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_simVtx = new TH1F("h_electron_relChIso05_simVtx", "h_electron_relChIso05_simVtx", 5000, 0, 5);
    h_electron_relChIso05_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT_simVtx = new TH1F("h_electron_relChIso05_dT_simVtx", "h_electron_relChIso05_dT_simVtx", 5000, 0, 5);
    h_electron_relChIso05_dT_simVtx->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_simVtx_barrel = new TH1F("h_electron_relChIso02_simVtx_barrel", "h_electron_relChIso02_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso02_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT_simVtx_barrel = new TH1F("h_electron_relChIso02_dT_simVtx_barrel", "h_electron_relChIso02_dT_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso02_dT_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_simVtx_barrel = new TH1F("h_electron_relChIso03_simVtx_barrel", "h_electron_relChIso03_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso03_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT_simVtx_barrel = new TH1F("h_electron_relChIso03_dT_simVtx_barrel", "h_electron_relChIso03_dT_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso03_dT_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_simVtx_barrel = new TH1F("h_electron_relChIso04_simVtx_barrel", "h_electron_relChIso04_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso04_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT_simVtx_barrel = new TH1F("h_electron_relChIso04_dT_simVtx_barrel", "h_electron_relChIso04_dT_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso04_dT_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_simVtx_barrel = new TH1F("h_electron_relChIso05_simVtx_barrel", "h_electron_relChIso05_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso05_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT_simVtx_barrel = new TH1F("h_electron_relChIso05_dT_simVtx_barrel", "h_electron_relChIso05_dT_simVtx_barrel", 5000, 0, 5);
    h_electron_relChIso05_dT_simVtx_barrel->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_simVtx_endcap = new TH1F("h_electron_relChIso02_simVtx_endcap", "h_electron_relChIso02_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso02_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso02_dT_simVtx_endcap = new TH1F("h_electron_relChIso02_dT_simVtx_endcap", "h_electron_relChIso02_dT_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso02_dT_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_simVtx_endcap = new TH1F("h_electron_relChIso03_simVtx_endcap", "h_electron_relChIso03_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso03_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT_simVtx_endcap = new TH1F("h_electron_relChIso03_dT_simVtx_endcap", "h_electron_relChIso03_dT_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso03_dT_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_simVtx_endcap = new TH1F("h_electron_relChIso04_simVtx_endcap", "h_electron_relChIso04_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso04_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT_simVtx_endcap = new TH1F("h_electron_relChIso04_dT_simVtx_endcap", "h_electron_relChIso04_dT_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso04_dT_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_simVtx_endcap = new TH1F("h_electron_relChIso05_simVtx_endcap", "h_electron_relChIso05_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso05_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT_simVtx_endcap = new TH1F("h_electron_relChIso05_dT_simVtx_endcap", "h_electron_relChIso05_dT_simVtx_endcap", 5000, 0, 5);
    h_electron_relChIso05_dT_simVtx_endcap->GetXaxis()->SetTitle("relative charged isolation");

    // ratio iso
    TH1F* h_electron_relChIso03_ratio = new TH1F("h_electron_relChIso03_ratio", "h_electron_relChIso03_ratio", 1000, 0, 2);
    h_electron_relChIso03_ratio->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

    TH1F* h_electron_relChIso03_ratio_barrel = new TH1F("h_electron_relChIso03_ratio_barrel", "h_electron_relChIso03_ratio_barrel", 1000, 0, 2);
    h_electron_relChIso03_ratio_barrel->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

    TH1F* h_electron_relChIso03_ratio_endcap = new TH1F("h_electron_relChIso03_ratio_endcap", "h_electron_relChIso03_ratio_endcap", 1000, 0, 2);
    h_electron_relChIso03_ratio_endcap->GetXaxis()->SetTitle("relChIso_{ZTcut}/relChIso_{Zcut}");

    /*
    TH1F* h_electron_relChIso03_diff = new TH1F("h_electron_relChIso03_diff", "h_electron_relChIso03_diff", 1000, 0, 1);
    h_electron_relChIso03_diff->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

    TH1F* h_electron_relChIso03_diff_barrel = new TH1F("h_electron_relChIso03_diff_barrel", "h_electron_relChIso03_diff_barrel", 1000, 0, 1);
    h_electron_relChIso03_diff_barrel->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

    TH1F* h_electron_relChIso03_diff_endcap = new TH1F("h_electron_relChIso03_diff_endcap", "h_electron_relChIso03_diff_endcap", 1000, 0, 1);
    h_electron_relChIso03_diff_endcap->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");*/

    // -- loop over events
    //int maxEntries = std::min(int(chain ->GetEntries()),1000000);
    Long64_t maxEntries = fReader.GetEntries(false);
    cout << "Number of events to be analyzed : " << maxEntries << endl;
    float wgt     = 1;
    float pt      = 0;
    float maxdz   = 0.5;
    float maxdxy  = 0.2;
    int   i       = 0;
    long  total   = 0;
    int   notzero = 0;
    while (fReader.Next()) {
        if (i % 10000 == 0)
            cout << "Analyzing event " << i << "\r" << flush;
        h_npu->Fill(*npu);
        int nElectronsInEvent = 0;
        for (unsigned int iele = 0; iele < electron_pt.GetSize(); iele++) {
            pt = electron_pt[iele];
            if (pt < 15.)
                continue;
            //if (electron_isLoose[iele] == 0)
            //    continue;
            total++;
            if (abs(electron_chIso03[iele] - 0.) < 0.00000001)
                continue;
            notzero++;
            // -- pt/eta reweighting
            if (applyPtReweighting) {
                int bin = hratio2D->FindBin(electron_eta[iele], pt);
                wgt     = hratio2D->GetBinContent(bin);
            }
            else {
                wgt = 1;
            }

            // -- prompt electrons or non-prompt electrons
            bool pass = false;
            if ((bool)usePromptelectrons && electron_isPrompt->at(iele))
                pass = true;
            if (!(bool)usePromptelectrons && (!electron_isPrompt->at(iele) && electron_isMatchedToGenJet->at(iele)))
                pass = true;

            if (!pass)
                continue;
            //cout << "pass:" << pass << endl;
            nElectronsInEvent++;

            bool pass3D = fabs(electron_dz3D[iele]) < maxdz && fabs(electron_dxy3D[iele]) < maxdxy && !(*vtx3D_isFake);
            bool pass4D = fabs(electron_dz4D[iele]) < maxdz && fabs(electron_dxy4D[iele]) < maxdxy && !(*vtx4D_isFake);
            bool pass4d = fabs(electron_dz[iele]) < maxdz && fabs(electron_dxy[iele]) < maxdxy && !(*vtx_isFake);  // pass4d used for the event not apply the dt cut
            if (nElectronsInEvent == 1) {
                if (pass3D)
                    h_vtx_dz3D->Fill(*vtx3D_z - *vtxGen_z);
                if (pass4D)
                    h_vtx_dz4D->Fill(*vtx4D_z - *vtxGen_z);
                if (pass4D)
                    h_vtx_dt4D->Fill(*vtx4D_t - (*vtxGen_t) * 1000000000.);
                if (pass4d)
                    h_vtx_dz->Fill(*vtx_z - *vtxGen_z);
                if (pass4d)
                    h_vtx_dt->Fill(*vtx_t - (*vtxGen_t) * 1000000000.);
                if (pass3D)
                    h_vtx_dz3D_pull->Fill((*vtx3D_z - *vtxGen_z) / *vtx3D_zErr);
                if (pass4D)
                    h_vtx_dz4D_pull->Fill((*vtx4D_z - *vtxGen_z) / *vtx4D_zErr);
                if (pass4D)
                    h_vtx_dt4D_pull->Fill((*vtx4D_t - (*vtxGen_t) * 1000000000.) / *vtx4D_tErr);
                if (pass4d)
                    h_vtx_dz_pull->Fill((*vtx_z - *vtxGen_z) / *vtx_zErr);
                if (pass4d)
                    h_vtx_dt_pull->Fill((*vtx_t - (*vtxGen_t) * 1000000000.) / *vtx_tErr);
            }

            h_electron_pt->Fill(pt, wgt);
            h_electron_eta->Fill(electron_eta[iele], wgt);
            h_electron_phi->Fill(electron_phi[iele], wgt);
            h2_electron_pt_vs_eta->Fill(electron_eta[iele], pt, wgt);

            // control plots
            if (pass3D && pass4d) {
                float chIsoRatio = (electron_chIso03_dT[iele]) / electron_chIso03[iele];
                if (electron_chIso03[iele] == 0.)
                    chIsoRatio = 1.;
                h_electron_relChIso03_ratio->Fill(chIsoRatio, wgt);
                if (fabs(electron_eta[iele]) < 1.479) {
                    h_electron_relChIso03_ratio_barrel->Fill(chIsoRatio, wgt);
                }
                else {
                    h_electron_relChIso03_ratio_endcap->Fill(chIsoRatio, wgt);
                }
            }

            // chiIso plots - all
            if (pass3D) {
                if (electron_chIso02[iele] / pt < 5)
                    h_electron_relChIso02->Fill(electron_chIso02[iele] / pt, wgt);
                else
                    h_electron_relChIso02->(4.999, wgt);

                if (electron_chIso03[iele] / pt < 5)
                    h_electron_relChIso03->Fill(electron_chIso03[iele] / pt, wgt);
                else
                    h_electron_relChIso03->(4.999, wgt);

                if (electron_chIso04[iele] / pt < 5)
                    h_electron_relChIso04->Fill(electron_chIso04[iele] / pt, wgt);
                else
                    h_electron_relChIso04->(4.999, wgt);

                if (electron_chIso05[iele] / pt < 5)
                    h_electron_relChIso05->Fill(electron_chIso05[iele] / pt, wgt);
                else
                    h_electron_relChIso05->(4.999, wgt);

                if (electron_chIso02_simVtx[iele] / pt < 5)
                    h_electron_relChIso02_simVtx->Fill(electron_chIso02_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso02_simVtx->(4.999, wgt);

                if (electron_chIso03_simVtx[iele] / pt < 5)
                    h_electron_relChIso03_simVtx->Fill(electron_chIso03_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso03_simVtx->(4.999, wgt);

                if (electron_chIso04_simVtx[iele] / pt < 5)
                    h_electron_relChIso04_simVtx->Fill(electron_chIso04_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso04_simVtx->(4.999, wgt);

                if (electron_chIso05_simVtx[iele] / pt < 5)
                    h_electron_relChIso05_simVtx->Fill(electron_chIso05_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso05_simVtx->(4.999, wgt);
            }
            if (pass4d) {
                /*h_electron_relChIso02_dT->Fill(electron_chIso02_dT[iele] / pt, wgt);
                h_electron_relChIso03_dT->Fill(electron_chIso03_dT[iele] / pt, wgt);
                h_electron_relChIso04_dT->Fill(electron_chIso04_dT[iele] / pt, wgt);
                h_electron_relChIso05_dT->Fill(electron_chIso05_dT[iele] / pt, wgt);*/
                if (electron_chIso02_dT[iele] / pt < 5)
                    h_electron_relChIso02_dT->Fill(electron_chIso02_dT[iele] / pt, wgt);
                else
                    h_electron_relChIso02_dT->(4.999, wgt);

                if (electron_chIso03_dT[iele] / pt < 5)
                    h_electron_relChIso03_dT->Fill(electron_chIso03_dT[iele] / pt, wgt);
                else
                    h_electron_relChIso03_dT->(4.999, wgt);

                if (electron_chIso04_dT[iele] / pt < 5)
                    h_electron_relChIso04_dT->Fill(electron_chIso04_dT[iele] / pt, wgt);
                else
                    h_electron_relChIso04_dT->(4.999, wgt);

                if (electron_chIso05_dT[iele] / pt < 5)
                    h_electron_relChIso05_dT->Fill(electron_chIso05_dT[iele] / pt, wgt);
                else
                    h_electron_relChIso05_dT->(4.999, wgt);

                if (electron_chIso02_dT_simVtx[iele] / pt < 5)
                    h_electron_relChIso02_dT_simVtx->Fill(electron_chIso02_dT_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso02_dT_simVtx->(4.999, wgt);

                if (electron_chIso03_dT_simVtx[iele] / pt < 5)
                    h_electron_relChIso03_dT_simVtx->Fill(electron_chIso03_dT_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso03_dT_simVtx->(4.999, wgt);

                if (electron_chIso04_dT_simVtx[iele] / pt < 5)
                    h_electron_relChIso04_dT_simVtx->Fill(electron_chIso04_dT_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso04_dT_simVtx->(4.999, wgt);

                if (electron_chIso05_dT_simVtx[iele] / pt < 5)
                    h_electron_relChIso05_dT_simVtx->Fill(electron_chIso05_dT_simVtx[iele] / pt, wgt);
                else
                    h_electron_relChIso05_dT_simVtx->(4.999, wgt);
            }
            if (pass4D) {
                /*h_electron_relChIso02_dT->Fill(electron_chIso02_dT[iele] / pt, wgt);
                h_electron_relChIso03_dT->Fill(electron_chIso03_dT[iele] / pt, wgt);
                h_electron_relChIso04_dT->Fill(electron_chIso04_dT[iele] / pt, wgt);
                h_electron_relChIso05_dT->Fill(electron_chIso05_dT[iele] / pt, wgt);*/
                if (electron_chIso02_dT_4D[iele] / pt < 5)
                    h_electron_relChIso02_dT_4D->Fill(electron_chIso02_dT_4D[iele] / pt, wgt);
                else
                    h_electron_relChIso02_dT_4D->(4.999, wgt);

                if (electron_chIso03_dT_4D[iele] / pt < 5)
                    h_electron_relChIso03_dT_4D->Fill(electron_chIso03_dT_4D[iele] / pt, wgt);
                else
                    h_electron_relChIso03_dT_4D->(4.999, wgt);

                if (electron_chIso04_dT_4D[iele] / pt < 5)
                    h_electron_relChIso04_dT_4D->Fill(electron_chIso04_dT_4D[iele] / pt, wgt);
                else
                    h_electron_relChIso04_dT_4D->(4.999, wgt);

                if (electron_chIso05_dT_4D[iele] / pt < 5)
                    h_electron_relChIso05_dT_4D->Fill(electron_chIso05_dT_4D[iele] / pt, wgt);
                else
                    h_electron_relChIso05_dT_4D->(4.999, wgt);
            }
            // barrel
            if (fabs(electron_eta[iele]) < 1.479) {
                if (pass3D) {
                    /*h_electron_relChIso02_barrel->Fill(electron_chIso02[iele] / pt, wgt);
                    h_electron_relChIso03_barrel->Fill(electron_chIso03[iele] / pt, wgt);
                    h_electron_relChIso04_barrel->Fill(electron_chIso04[iele] / pt, wgt);
                    h_electron_relChIso05_barrel->Fill(electron_chIso05[iele] / pt, wgt);*/
                    if (electron_chIso02[iele] / pt < 5)
                        h_electron_relChIso02_barrel->Fill(electron_chIso02[iele] / pt, wgt);
                    else
                        h_electron_relChIso02_barrel->(4.999, wgt);

                    if (electron_chIso03[iele] / pt < 5)
                        h_electron_relChIso03_barrel->Fill(electron_chIso03[iele] / pt, wgt);
                    else
                        h_electron_relChIso03_barrel->(4.999, wgt);

                    if (electron_chIso04[iele] / pt < 5)
                        h_electron_relChIso04_barrel->Fill(electron_chIso04[iele] / pt, wgt);
                    else
                        h_electron_relChIso04_barrel->(4.999, wgt);

                    if (electron_chIso05[iele] / pt < 5)
                        h_electron_relChIso05_barrel->Fill(electron_chIso05[iele] / pt, wgt);
                    else
                        h_electron_relChIso05_barrel->(4.999, wgt);
                }
                if (pass4d) {
                    /*h_electron_relChIso02_dT_barrel->Fill(electron_chIso02_dT[iele] / pt, wgt);
                    h_electron_relChIso03_dT_barrel->Fill(electron_chIso03_dT[iele] / pt, wgt);
                    h_electron_relChIso04_dT_barrel->Fill(electron_chIso04_dT[iele] / pt, wgt);
                    h_electron_relChIso05_dT_barrel->Fill(electron_chIso05_dT[iele] / pt, wgt);*/
                    if (electron_chIso02_dT[iele] / pt < 5)
                        h_electron_relChIso02_dT_barrel->Fill(electron_chIso02_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso02_dT_barrel->(4.999, wgt);

                    if (electron_chIso03_dT[iele] / pt < 5)
                        h_electron_relChIso03_dT_barrel->Fill(electron_chIso03_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso03_dT_barrel->(4.999, wgt);

                    if (electron_chIso04_dT[iele] / pt < 5)
                        h_electron_relChIso04_dT_barrel->Fill(electron_chIso04_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso04_dT_barrel->(4.999, wgt);

                    if (electron_chIso05_dT[iele] / pt < 5)
                        h_electron_relChIso05_dT_barrel->Fill(electron_chIso05_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso05_dT_barrel->(4.999, wgt);
                }
            }
            // endcap
            else {
                if (pass3D) {
                    /*h_electron_relChIso02_endcap->Fill(electron_chIso02[iele] / pt, wgt);
                    h_electron_relChIso03_endcap->Fill(electron_chIso03[iele] / pt, wgt);
                    h_electron_relChIso04_endcap->Fill(electron_chIso04[iele] / pt, wgt);
                    h_electron_relChIso05_endcap->Fill(electron_chIso05[iele] / pt, wgt);*/
                    if (electron_chIso02[iele] / pt < 5)
                        h_electron_relChIso02_endcap->Fill(electron_chIso02[iele] / pt, wgt);
                    else
                        h_electron_relChIso02_endcap->(4.999, wgt);

                    if (electron_chIso03[iele] / pt < 5)
                        h_electron_relChIso03_endcap->Fill(electron_chIso03[iele] / pt, wgt);
                    else
                        h_electron_relChIso03_endcap->(4.999, wgt);

                    if (electron_chIso04[iele] / pt < 5)
                        h_electron_relChIso04_endcap->Fill(electron_chIso04[iele] / pt, wgt);
                    else
                        h_electron_relChIso04_endcap->(4.999, wgt);

                    if (electron_chIso05[iele] / pt < 5)
                        h_electron_relChIso05_endcap->Fill(electron_chIso05[iele] / pt, wgt);
                    else
                        h_electron_relChIso05_endcap->(4.999, wgt);
                }
                if (pass4d) {
                    /*h_electron_relChIso02_dT_endcap->Fill(electron_chIso02_dT[iele] / pt, wgt);
                    h_electron_relChIso03_dT_endcap->Fill(electron_chIso03_dT[iele] / pt, wgt);
                    h_electron_relChIso04_dT_endcap->Fill(electron_chIso04_dT[iele] / pt, wgt);
                    h_electron_relChIso05_dT_endcap->Fill(electron_chIso05_dT[iele] / pt, wgt);*/
                    if (electron_chIso02_dT[iele] / pt < 5)
                        h_electron_relChIso02_dT_endcap->Fill(electron_chIso02_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso02_dT_endcap->(4.999, wgt);

                    if (electron_chIso03_dT[iele] / pt < 5)
                        h_electron_relChIso03_dT_endcap->Fill(electron_chIso03_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso03_dT_endcap->(4.999, wgt);

                    if (electron_chIso04_dT[iele] / pt < 5)
                        h_electron_relChIso04_dT_endcap->Fill(electron_chIso04_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso04_dT_endcap->(4.999, wgt);

                    if (electron_chIso05_dT[iele] / pt < 5)
                        h_electron_relChIso05_dT_endcap->Fill(electron_chIso05_dT[iele] / pt, wgt);
                    else
                        h_electron_relChIso05_dT_endcap->(4.999, wgt);
                }
            }

        }  // end loop over electrons
        for (unsigned int idr = 0; idr < drep.GetSize(); idr++) {
            h_drep->Fill(drep[idr], wgt);
        }
        i++;
    }  // end loop over events
    cout << "total:" << total << endl;
    cout << "notzero:" << notzero << endl;
    cout << endl;

    cout << "Saving histograms in file "
         << "out_30ps_" << samplename << ".root" << endl;

    // -- save output file
    TFile* fout = new TFile("out_" + timeResolution + samplename, "recreate");
    h_npu->Write();
    h_electron_pt->Write();
    h_electron_eta->Write();
    h_electron_phi->Write();
    h2_electron_pt_vs_eta->Write();
    h_drep->Write();

    h_vtx_dz3D->Write();
    h_vtx_dz4D->Write();
    h_vtx_dt4D->Write();
    h_vtx_dz->Write();
    h_vtx_dt->Write();

    h_vtx_dz3D_pull->Write();
    h_vtx_dz4D_pull->Write();
    h_vtx_dt4D_pull->Write();
    h_vtx_dz_pull->Write();
    h_vtx_dt_pull->Write();
    /*h_electron_relChIso03_diff->Write();
    h_electron_relChIso03_diff_barrel->Write();
    h_electron_relChIso03_diff_endcap->Write();*/

    h_electron_relChIso03_ratio->Write();
    h_electron_relChIso03_ratio_barrel->Write();
    h_electron_relChIso03_ratio_endcap->Write();

    h_electron_relChIso02->Write();
    h_electron_relChIso02_barrel->Write();
    h_electron_relChIso02_endcap->Write();
    h_electron_relChIso02_dT->Write();
    h_electron_relChIso02_dT_barrel->Write();
    h_electron_relChIso02_dT_endcap->Write();
    h_electron_relChIso02_dT_4D->Write();

    h_electron_relChIso03->Write();
    h_electron_relChIso03_barrel->Write();
    h_electron_relChIso03_endcap->Write();
    h_electron_relChIso03_dT->Write();
    h_electron_relChIso03_dT_barrel->Write();
    h_electron_relChIso03_dT_endcap->Write();
    h_electron_relChIso03_dT_4D->Write();

    h_electron_relChIso04->Write();
    h_electron_relChIso04_barrel->Write();
    h_electron_relChIso04_endcap->Write();
    h_electron_relChIso04_dT->Write();
    h_electron_relChIso04_dT_barrel->Write();
    h_electron_relChIso04_dT_endcap->Write();
    h_electron_relChIso04_dT_4D->Write();

    h_electron_relChIso05->Write();
    h_electron_relChIso05_barrel->Write();
    h_electron_relChIso05_endcap->Write();
    h_electron_relChIso05_dT->Write();
    h_electron_relChIso05_dT_barrel->Write();
    h_electron_relChIso05_dT_endcap->Write();
    h_electron_relChIso05_dT_4D->Write();

    h_electron_relChIso02_simVtx->Write();
    h_electron_relChIso02_dT_simVtx->Write();
    /*h_electron_relChIso02_simVtx_barrel->Write();
    h_electron_relChIso02_dT_simVtx_barrel->Write();
    h_electron_relChIso02_simVtx_endcap->Write();
    h_electron_relChIso02_dT_simVtx_endcap->Write();*/

    h_electron_relChIso03_simVtx->Write();
    h_electron_relChIso03_dT_simVtx->Write();
    /*h_electron_relChIso03_simVtx_barrel->Write();
    h_electron_relChIso03_dT_simVtx_barrel->Write();
    h_electron_relChIso03_simVtx_endcap->Write();
    h_electron_relChIso03_dT_simVtx_endcap->Write();*/

    h_electron_relChIso04_simVtx->Write();
    h_electron_relChIso04_dT_simVtx->Write();
    /*h_electron_relChIso04_simVtx_barrel->Write();
    h_electron_relChIso04_dT_simVtx_barrel->Write();
    h_electron_relChIso04_simVtx_endcap->Write();
    h_electron_relChIso04_dT_simVtx_endcap->Write();*/

    h_electron_relChIso05_simVtx->Write();
    h_electron_relChIso05_dT_simVtx->Write();
    /*h_electron_relChIso05_simVtx_barrel->Write();
    h_electron_relChIso05_dT_simVtx_barrel->Write();
    h_electron_relChIso05_simVtx_endcap->Write();
    h_electron_relChIso05_dT_simVtx_endcap->Write();*/

    fout->Close();
}
