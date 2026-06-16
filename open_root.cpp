#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include <TLegend.h>


void printBranches(TObjArray* branches, int depth = 0) {
    std::string indent(depth * 2, ' ');
    TIter next(branches);
    TBranch* branch;

    while ((branch = (TBranch*)next())) {
        std::cout << indent << branch->GetName()
                  << std::endl;

        // Recurse into sub-branches
        if (branch->GetListOfBranches()->GetEntries() > 0)
            printBranches(branch->GetListOfBranches(), depth + 1);
    }
}

int main() {

    double pi = TMath::Pi();
    double det_a = 4*pi*TMath::Power(10,4);

    TFile* file = TFile::Open("DATA/NuBeam_production_BooNE_50m_I174000A_0.dk2nu.root", "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: could not open file!" << std::endl;
        return 1;
    }

    // Print top-level keys
    TIter nextKey(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextKey())) {
        std::cout << "Name: " << key->GetName()
                  << "  Class: " << key->GetClassName()
                  << std::endl;
    }

    // Get the TTree and print its branches
    TTree* tree = (TTree*)file->Get("dk2nuTree");
    if (!tree) {
        std::cerr << "Error: could not find tree!" << std::endl;
        return 1;
    }

    printBranches(tree->GetListOfBranches());

    // Get the number of protons on target
    TTree* meta = (TTree*)file->Get("dkmetaTree");
    if (!tree) {
        std::cerr << "Error: could not find tree!" << std::endl;
        return 1;
    }

    // Simple Flux Histogram

    // Read Total Protons on Target from Metadata
    Double_t total_pots = 0;

    TTreeReader meta_reader(meta);
    TTreeReaderValue<Double_t> pots(meta_reader, "pots");

    while (meta_reader.Next()) {
    total_pots += *pots;
    }
    std::cout << "Total POTs: " << total_pots << std::endl;

    // Read Event Data
    TTreeReader data_reader(tree);
    TTreeReaderArray<Double_t> E(data_reader, "nuray.E");
    TTreeReaderArray<Double_t> wgts(data_reader, "nuray.wgt");
    TTreeReaderValue<Double_t> nimp_wgts(data_reader, "decay.nimpwt");
    TTreeReaderValue<Int_t> pdg(data_reader, "decay.ntype");


    // Define histograms for each neutrino type
    TH1F* hist_numu     = new TH1F("hist_numu",     "Flux by Particle Type; Energy (GeV); Flux (#nu / cm^{2} / GeV / POT)",          100, 0, 10);
    TH1F* hist_numubar  = new TH1F("hist_numubar",  "#bar{#nu}_{#mu}; Energy (GeV);  Flux (#nu / cm^{2} / GeV / POT)",    100, 0, 10);
    TH1F* hist_nue      = new TH1F("hist_nue",      "#nu_{e}; Energy (GeV);  Flux (#nu / cm^{2} / GeV / POT)",             100, 0, 10);
    TH1F* hist_nuebar   = new TH1F("hist_nuebar",   "#bar{#nu}_{e}; Energy (GeV);  Flux (#nu / cm^{2} / GeV / POT)",       100, 0, 10);

    hist_numu->Sumw2();
    hist_numubar->Sumw2();
    hist_nue->Sumw2();
    hist_nuebar->Sumw2();

    while (data_reader.Next()) {
        Double_t energy = E[5];
        Double_t weight = wgts[5] * (*nimp_wgts) / (det_a * total_pots);

        // PDG codes: 14 = numu, -14 = numubar, 12 = nue, -12 = nuebar
        switch (*pdg) {
            case  14: hist_numu->Fill(energy, weight);    break;
            case -14: hist_numubar->Fill(energy, weight); break;
            case  12: hist_nue->Fill(energy, weight);     break;
            case -12: hist_nuebar->Fill(energy, weight);  break;
        }
    }

    // Set colours
    hist_numu->SetLineColor(kBlue);
    hist_numubar->SetLineColor(kRed);
    hist_nue->SetLineColor(kGreen+2);
    hist_nuebar->SetLineColor(kMagenta);

    // Scale line width for visibility
    hist_numu->SetLineWidth(2);
    hist_numubar->SetLineWidth(2);
    hist_nue->SetLineWidth(2);
    hist_nuebar->SetLineWidth(2);

    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLogy();

    hist_numu->SetStats(0);
    hist_numubar->SetStats(0);
    hist_nue->SetStats(0);
    hist_nuebar->SetStats(0);

    // Draw on the same canvas with SAME after the first
    hist_numu->Draw("HIST");
    hist_numubar->Draw("HIST SAME");
    hist_nue->Draw("HIST SAME");
    hist_nuebar->Draw("HIST SAME");

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // (x1, y1, x2, y2) in pad coordinates
    legend->AddEntry(hist_numu,    "#nu_{#mu}",          "l");
    legend->AddEntry(hist_numubar, "#bar{#nu}_{#mu}",    "l");
    legend->AddEntry(hist_nue,     "#nu_{e}",             "l");
    legend->AddEntry(hist_nuebar,  "#bar{#nu}_{e}",       "l");
    legend->Draw();

    canvas->SaveAs("flux_by_flavour.png");

    file->Close();
    return 0;
}
