#include <iostream>
#include "TFile.h"
#include "TTree.h"

int main() {
    TFile* file = TFile::Open(
        "NuBeam_production_BooNE_50m_I174000A_0.dk2nu.root",
        "READ");

    TTree* tree = nullptr;
    file->GetObject("dk2nuTree", tree);

    std::cout << "Entries = "
              << tree->GetEntries()
              << std::endl;

    std::cout << "Maximum ancestor count = "
              << tree->GetMaximum("ancestor_")
              << std::endl;

    file->Close();
}