
#include <iostream>
#include <memory>
#include "TFile.h"

using namespace std;

std::unique_ptr<TFile> myFile( TFile::Open("file.root", "RECREATE") );

