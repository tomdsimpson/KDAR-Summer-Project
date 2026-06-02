#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

int main() {
    // Open the ROOT file
    TFile* file = TFile::Open("NuBeam_production_BooNE_50m_I174000A_0.dk2nu.root", "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file!" << std::endl;
        return 1;
    }
}