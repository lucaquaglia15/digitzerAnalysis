//Root & C++ classes
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"

#include "TAxis.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include <iterator>
#include <numeric>
#include <map>
#include "Fit/FitResult.h"
#include <array>

//My classes
#include "utilities.h"

using namespace std;

const string txtExt = ".txt"; //txt extension
const string rootExt = ".root"; //root extension
const string folder = "/home/luca/cernbox/PhD/digitizer904/000"; //folder with data
const string trig = "TR_0_"; //digitized trigger file name
const string wave = "wave_"; //file name + number from 0 to 15 for the strip number
const int muonMin = 150, muonMax = 250;
const int noiseMin = 400, noiseMax = 600;
bool verbose = false; //true = prints debug, false = no printout

void analyzeDigitizer(const int run) {

	cout << "Starting digitizer data analysis!" << endl;

	//From struct in .h file
	struct chamber1G testedDetector;

	//Info on the chamber
	cout << "Detector: " << testedDetector.name << endl;
	cout << "Number of strips: " << testedDetector.stripNum << endl;
	int readoutStrips = *(&testedDetector.strips + 1) - testedDetector.strips;
	cout << "Number of strips readout: " << readoutStrips << endl;
	cout << "Sampling frequency: " << testedDetector.sampFreq << "Gs/s" << endl;
	cout << "Number of samples: " << testedDetector.numSamp << endl;

	//Convert .txt files in .root
	cout << "Converting .txt data in .root format..." << endl;

	//Enter the scan folder
	string scanFolder = folder + to_string(run) + "/";
	gSystem->cd(scanFolder.c_str());

	//Count how many folders are there (# of HV points) -> proxy = how many .root files are there
	int file_count = hvCounter(scanFolder, rootExt, verbose);
    cout << "Number of HV points: " << file_count << endl;

    //Call tree producer function once per HV point
    for (int i = 0; i < file_count; i++) {
    	string hvPoint = scanFolder + "HV" + to_string(i+1) + "_DIGITIZER";
    	//Produce tree
    	//treeProducer(readoutStrips,hvPoint,trig,wave,i,verbose);
    }

    vector<double> eff, hvEff;

    for (int i = 0; i < file_count; i++) {
    	string hvPoint = scanFolder + "HV" + to_string(i+1) + "_DIGITIZER";
    	eff.push_back(analyzer(readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax));
    }

    //Read hv effective from .root CAEN file
    hvEff = hvReader(scanFolder, rootExt, file_count, verbose);


    //Plot efficiency(HV) curve
    TCanvas *cEff = new TCanvas();
    cEff->cd();

    TGraphErrors *effHv = new TGraphErrors(eff.size(),&hvEff[0],&eff[0],NULL,NULL);
    effHv->SetMarkerStyle(8);
    effHv->SetMarkerSize(1);
    effHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    effHv->GetYaxis()->SetTitle("Eff [%]");
    effHv->Draw("AP");
}