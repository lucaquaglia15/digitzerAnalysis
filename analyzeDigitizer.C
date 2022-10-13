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
const int muonMin = 150, muonMax = 250; //1 Gs/s
//const int muonMin = 350, muonMax = 500; //5 Gs/s - index of the vector (5Gs/s with 1024 samples = 204.8 ns window. 204.8/1024 (samples) = 0,2 ns/sample). 
//Time window at 5 Gs/s = 70-100 ns
const int noiseMin = 50, noiseMax = 100; //1 Gs/s
//const int noiseMin = 100, noiseMax = 200; //5 Gs/s - time for noise 140-180 ns -> now 20 - 40 ns

TF1* fit_eff;

bool verbose = false; //true = prints debug, false = no printout

void analyzeDigitizer(const int run, const bool makeTree) {

	cout << "Starting digitizer data analysis!" << endl;

	cout << "Run number: " << run << endl;

	//From struct in .h file
	struct chamber1G testedDetector; // 1 Gs/s
	//struct chamber testedDetector; // 5 Gs/s

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
    	//Produce tree if makeTree is true, otherwise tree already created
    	cout << "Converting hv point " << i+1 << endl;
    	if (makeTree) 
            treeProducer(readoutStrips,hvPoint,trig,wave,i,verbose);
    }

    //Compute time vector for each event
    vector<double> eff, hvEff;
    float time = 0;
    vector<float> times;
    for (int i = 0; i < testedDetector.numSamp; i++) {
        if (i == 0) times.push_back(0);
        else {
            time = time + (testedDetector.numSamp/testedDetector.sampFreq)/testedDetector.numSamp;
            times.push_back(time);
        }
    }

    //full analysis 
    /*for (int i = 0; i < file_count; i++) {
    	string hvPoint = scanFolder + "HV" + to_string(i+1) + "_DIGITIZER";
    	//eff.push_back(analyzer(readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax));
        eff.push_back(analyzer(times,readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax));
    }*/

    //analyze only hv 10 for a test
    for (int i = 9; i < file_count; i++) {
        string hvPoint = scanFolder + "HV" + to_string(i+1) + "_DIGITIZER";
        //eff.push_back(analyzer(readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax));
        eff.push_back(analyzer(times,readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax));
    }

    //Read hv effective from .root CAEN file
    /*hvEff = hvReader(scanFolder, run, file_count, verbose);

    //Fit function for eff(HV) curve + set param limits and names
   	fit_eff = new TF1("fit_eff","[0]/(1+TMath::Exp(-[1]*(x-[2])))",hvEff.front(),hvEff.back());

	fit_eff->SetParLimits(0,45.,110.);
	//fit_eff->SetParLimits(1,0.,0.5);
	//fit_eff->SetParLimits(2,9900.,11000.);

	double effmax = fit_eff->GetParameter(0); //Max efficiency
	double lambda = fit_eff->GetParameter(1); //Lambda parameter
	double HV50 = fit_eff->GetParameter(2); //HV for 50% eff

    //Plot efficiency(HV) curve
    TCanvas *cEff = new TCanvas();
    cEff->cd();

    TGraphErrors *effHv = new TGraphErrors(eff.size(),&hvEff[0],&eff[0],NULL,NULL);
    effHv->SetMarkerStyle(8);
    effHv->SetMarkerSize(1);
    effHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    effHv->GetYaxis()->SetTitle("Eff [%]");
    effHv->GetYaxis()->SetRangeUser(0,100);
    effHv->Fit(fit_eff,"RM+");
    effHv->Draw("AP");

    double WP = (TMath::Log(19)/lambda)+HV50+150; //WP with CMS definition
	double effwp = fit_eff->Eval(WP); //Efficiency at WP

    //Plot I(HV) curve
    TCanvas *cCurr = new TCanvas();
    cCurr->cd();

    TGraphErrors *currHv = new TGraphErrors(curr.size(),&hvEff[0],&curr[0],NULL,NULL);
    currHv->SetMarkerStyle(8);
    currHv->SetMarkerSize(1);
    currHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    currHv->GetYaxis()->SetTitle("I [#muA]");
    currHv->GetYaxis()->SetRangeUser(curr.front() - 0.1*curr.front(),curr.back() + 0.1*curr.back());
    currHv->Draw("AP");

    TFile *fAnalysisOut = new TFile(("Analysis_run_"+to_string(run)+".root").c_str(),"RECREATE");
    fAnalysisOut->cd();
    cEff->Write(("Efficiency_run_"+to_string(run)).c_str());
    //cCurr->Write("Current_run_"+to_string(run));
    fAnalysisOut->Close();*/
}
