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
    TFile* ipt = TFile::Open(folder + "/eleiso_" + samplename + ".root");
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

    TH1F* h_electron_relChIso03 = new TH1F("h_electron_relChIso03", "h_electron_relChIso03", 5000, 0, 5);
    h_electron_relChIso03->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso03_dT = new TH1F("h_electron_relChIso03_dT", "h_electron_relChIso03_dT", 5000, 0, 5);
    h_electron_relChIso03_dT->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04 = new TH1F("h_electron_relChIso04", "h_electron_relChIso04", 5000, 0, 5);
    h_electron_relChIso04->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso04_dT = new TH1F("h_electron_relChIso04_dT", "h_electron_relChIso04_dT", 5000, 0, 5);
    h_electron_relChIso04_dT->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05 = new TH1F("h_electron_relChIso05", "h_electron_relChIso05", 5000, 0, 5);
    h_electron_relChIso05->GetXaxis()->SetTitle("relative charged isolation");

    TH1F* h_electron_relChIso05_dT = new TH1F("h_electron_relChIso05_dT", "h_electron_relChIso05_dT", 5000, 0, 5);
    h_electron_relChIso05_dT->GetXaxis()->SetTitle("relative charged isolation");

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

    //
    TH1F* h_electron_relChIso03_diff = new TH1F("h_electron_relChIso03_diff", "h_electron_relChIso03_diff", 1000, 0, 1);
    h_electron_relChIso03_diff->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

    TH1F* h_electron_relChIso03_diff_barrel = new TH1F("h_electron_relChIso03_diff_barrel", "h_electron_relChIso03_diff_barrel", 1000, 0, 1);
    h_electron_relChIso03_diff_barrel->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

    TH1F* h_electron_relChIso03_diff_endcap = new TH1F("h_electron_relChIso03_diff_endcap", "h_electron_relChIso03_diff_endcap", 1000, 0, 1);
    h_electron_relChIso03_diff_endcap->GetXaxis()->SetTitle("(relChIso_{Zcut}-relChIso_{ZTcut})/relChIso_{Zcut}");

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

            h_electron_pt->Fill(pt, wgt);
            h_electron_eta->Fill(electron_eta[iele], wgt);
            h_electron_phi->Fill(electron_phi[iele], wgt);
            h2_electron_pt_vs_eta->Fill(electron_eta[iele], pt, wgt);

            // control plots
            if (fabs(electron_dz3D[iele]) < maxdz && fabs(electron_dxy3D[iele]) < maxdxy && fabs(electron_dz4D[iele]) < maxdz && fabs(electron_dxy4D[iele]) < maxdxy) {
                float chIsoDiff = (electron_chIso03[iele] - electron_chIso03_dT[iele]) / electron_chIso03[iele];
                if (electron_chIso03[iele] == 0 && electron_chIso03[iele] == 0.)
                    chIsoDiff = 0.;
                h_electron_relChIso03_diff->Fill(chIsoDiff, wgt);
                if (fabs(electron_eta[iele]) < 1.4442) {
                    h_electron_relChIso03_diff_barrel->Fill(chIsoDiff, wgt);
                }
                else {
                    h_electron_relChIso03_diff_endcap->Fill(chIsoDiff, wgt);
                }
            }

            // all
            if (fabs(electron_dz3D[iele]) < maxdz && fabs(electron_dxy3D[iele]) < maxdxy) {
                h_electron_relChIso02->Fill(electron_chIso02[iele] / pt, wgt);
                h_electron_relChIso03->Fill(electron_chIso03[iele] / pt, wgt);
                h_electron_relChIso04->Fill(electron_chIso04[iele] / pt, wgt);
                h_electron_relChIso05->Fill(electron_chIso05[iele] / pt, wgt);
            }
            if (fabs(electron_dz4D[iele]) < maxdz && fabs(electron_dxy4D[iele]) < maxdxy) {
                h_electron_relChIso02_dT->Fill(electron_chIso02_dT[iele] / pt, wgt);
                h_electron_relChIso03_dT->Fill(electron_chIso03_dT[iele] / pt, wgt);
                h_electron_relChIso04_dT->Fill(electron_chIso04_dT[iele] / pt, wgt);
                h_electron_relChIso05_dT->Fill(electron_chIso05_dT[iele] / pt, wgt);
            }

            // barrel
            if (fabs(electron_eta[iele]) < 1.4442) {
                if (fabs(electron_dz3D[iele]) < maxdz && fabs(electron_dxy3D[iele]) < maxdxy) {
                    h_electron_relChIso02_barrel->Fill(electron_chIso02[iele] / pt, wgt);
                    h_electron_relChIso03_barrel->Fill(electron_chIso03[iele] / pt, wgt);
                    h_electron_relChIso04_barrel->Fill(electron_chIso04[iele] / pt, wgt);
                    h_electron_relChIso05_barrel->Fill(electron_chIso05[iele] / pt, wgt);
                }
                if (fabs(electron_dz4D[iele]) < maxdz && fabs(electron_dxy4D[iele]) < maxdxy) {
                    h_electron_relChIso02_dT_barrel->Fill(electron_chIso02_dT[iele] / pt, wgt);
                    h_electron_relChIso03_dT_barrel->Fill(electron_chIso03_dT[iele] / pt, wgt);
                    h_electron_relChIso04_dT_barrel->Fill(electron_chIso04_dT[iele] / pt, wgt);
                    h_electron_relChIso05_dT_barrel->Fill(electron_chIso05_dT[iele] / pt, wgt);
                }
            }
            // endcap
            else {
                if (fabs(electron_dz3D[iele]) < maxdz && fabs(electron_dxy3D[iele]) < maxdxy) {
                    h_electron_relChIso02_endcap->Fill(electron_chIso02[iele] / pt, wgt);
                    h_electron_relChIso03_endcap->Fill(electron_chIso03[iele] / pt, wgt);
                    h_electron_relChIso04_endcap->Fill(electron_chIso04[iele] / pt, wgt);
                    h_electron_relChIso05_endcap->Fill(electron_chIso05[iele] / pt, wgt);
                }
                if (fabs(electron_dz4D[iele]) < maxdz && fabs(electron_dxy4D[iele]) < maxdxy) {
                    h_electron_relChIso02_dT_endcap->Fill(electron_chIso02_dT[iele] / pt, wgt);
                    h_electron_relChIso03_dT_endcap->Fill(electron_chIso03_dT[iele] / pt, wgt);
                    h_electron_relChIso04_dT_endcap->Fill(electron_chIso04_dT[iele] / pt, wgt);
                    h_electron_relChIso05_dT_endcap->Fill(electron_chIso05_dT[iele] / pt, wgt);
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
    TFile* fout = new TFile("out_30ps_" + samplename + ".root", "recreate");
    h_npu->Write();
    h_electron_pt->Write();
    h_electron_eta->Write();
    h_electron_phi->Write();
    h2_electron_pt_vs_eta->Write();
    h_drep->Write();
    h_electron_relChIso03_diff->Write();
    h_electron_relChIso03_diff_barrel->Write();
    h_electron_relChIso03_diff_endcap->Write();

    h_electron_relChIso02->Write();
    h_electron_relChIso02_barrel->Write();
    h_electron_relChIso02_endcap->Write();
    h_electron_relChIso02_dT->Write();
    h_electron_relChIso02_dT_barrel->Write();
    h_electron_relChIso02_dT_endcap->Write();

    h_electron_relChIso03->Write();
    h_electron_relChIso03_barrel->Write();
    h_electron_relChIso03_endcap->Write();
    h_electron_relChIso03_dT->Write();
    h_electron_relChIso03_dT_barrel->Write();
    h_electron_relChIso03_dT_endcap->Write();

    h_electron_relChIso04->Write();
    h_electron_relChIso04_barrel->Write();
    h_electron_relChIso04_endcap->Write();
    h_electron_relChIso04_dT->Write();
    h_electron_relChIso04_dT_barrel->Write();
    h_electron_relChIso04_dT_endcap->Write();

    h_electron_relChIso05->Write();
    h_electron_relChIso05_barrel->Write();
    h_electron_relChIso05_endcap->Write();
    h_electron_relChIso05_dT->Write();
    h_electron_relChIso05_dT_barrel->Write();
    h_electron_relChIso05_dT_endcap->Write();

    fout->Close();
}
