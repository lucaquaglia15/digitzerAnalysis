#ifndef UTILITIES_H
#define UTILITIES_H

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TH1F.h"
#include <stdlib.h>
#include <fstream>
#include <numeric>
#include <vector>

using namespace std;

struct chamber { //CMS digitizer
	std::string name = "ALICE-2-0"; 
	int stripNum = 16, numSamp = 1024, resolution = 12, offsetBits = 16; //resolution and offsetBits are in bits
	int offset = 0x7fff; //DC offeset in DAC (to be converted in mV)
	int strips[8] = {1,2,3,4,5,6,7,8};
	float sampFreq = 2.5, vpp = 1000.0; //Gs/s, mV
};

struct chamber1G { //CMS digitizer - 1 Gs/s
	std::string name = "ALICE-2-0"; 
	int stripNum = 16, numSamp = 1024, resolution = 12, offsetBits = 16; //resolution and offsetBits are in bits
	int offset = 0x7fff; //DC offeset in DAC (to be converted in mV)
	int strips[8] = {1,2,3,4,5,6,7,8};
	float sampFreq = 1., vpp = 1000.0; //Gs/s, mV
};


struct chamberEPDT { //EPDT digitizer
	std::string name = "ALICE-2-0";
	int stripNum, numSamp;
	int strips[7] = {1,2,3,4,5,6,7};
	float sampFreq = 2.5; //Gs/s
};

void treeProducer(int numFiles, string folderPath, string trigFile, string waveFile, int hv, bool verbose);

int hvCounter(string folderPath, string ext, bool verbose);

vector <double> hvReader(string scanPath, string ext, int numfiles, bool verbose);

double analyzer(int numFiles, string folderPath, int hv, bool verbose, struct chamber1G, int muonMin, int muonMax, int noiseMin, int noiseMax);
#endif