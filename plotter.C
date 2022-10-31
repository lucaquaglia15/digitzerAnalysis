#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMUltiGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "analysis.h"
#include "Riostream.h"
#include "TSystem.h"
#include <fstream>
#include <vector>

using namespace std;

void plotter (const int runNumber) {
	
	int run, samples;
	string mix;
	double filter, hfo, co2, r134a, iso, sf6;
	vector<int> vRun, vSamples;
	vector<string> vMix;
	vector<double> vFilter, vHfo, vCo2, vR134a, vIso, vSf6;

	//open run file and get all info regarding the run
	ifstream hRuns;
	hRuns.open("scans.txt");

	while (hRuns >> run >> mix >> filter >> samples >> hfo >> co2 >> r134a >> iso >> sf6) {
		vRun.push_back(run);
		vMix.push_back(mix);
		vFilter.push_back(filter);
		vSamples.push_back(samples);
		vHfo.push_back(hfo);
		vCo2.push_back(co2);
		vR134a.push_back(r134a);
		vIso.push_back(iso);
		vSf6.push_back(sf6);
	}

	//Open path to root files from which to get the data
	string folder = "/home/luca/cernbox/PhD/digitizer904/000" + to_string(runNumber);
	gSystem->cd(folder.c_str());

	

}