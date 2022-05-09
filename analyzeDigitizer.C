#include "TTree.h"
#include "TBranch.h"

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
#include "TSystem.h"
#include "TAxis.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include <iterator>
#include <numeric>
#include <map>
#include "Fit/FitResult.h"

#include "utilities.h"

using namespace std;

void analyzeDigitizer() {

	cout << "Starting digitizer data analysis!" << endl;

	//From struct in .h file
	struct chamber testedDetector;

	cout << "Detector: " << testedDetector.name << endl;
	cout << "Number of strips: " << testedDetector.stripNum << endl;
	cout << "Sampling frequency: " << testedDetector.sampFreq << "Gs/s" << endl;
	cout << "Number of samples: " << testedDetector.numSamp << endl;

	cout << "Converting .txt data in .root format..." << endl;

}