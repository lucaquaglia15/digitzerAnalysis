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
#include <tuple>

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
    remove(("Analysis_run_"+to_string(run)+".root").c_str());

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

    tuple<double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double> analysisResults;
    tuple<vector <double>,vector <double>,vector <double>,vector <double>,vector <double>,vector <double>> caenParam;
    vector<double> eff, eEff,muCS,eMuCS,muCM,eMuCM,gammaCS,eGammaCS,gammaCM,eGammaCM,gammaRate,eGammaRate,muCharge,eMuCharge,gammaCharge,eGammaCharge,tot,eTot; 
    vector<double> hvEff,eHVeff,iMon,eImon,hvMon,eHVmon;

    //analyze a single hv point for a test
    for (int i = 0; i < file_count; i++) {
        string hvPoint = scanFolder + "HV" + to_string(i+1) + "_DIGITIZER";
        //eff.push_back(analyzer(times,readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax));
        //analyzer(times,readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax);
        analysisResults = analyzer(times,readoutStrips,hvPoint,i,verbose,testedDetector, muonMin, muonMax, noiseMin, noiseMax);
        eff.push_back(get<0>(analysisResults));
        eEff.push_back(get<1>(analysisResults));
        muCS.push_back(get<2>(analysisResults));
        eMuCS.push_back(get<3>(analysisResults));
        muCM.push_back(get<4>(analysisResults));
        eMuCM.push_back(get<5>(analysisResults));
        gammaCS.push_back(get<6>(analysisResults));
        eGammaCS.push_back(get<7>(analysisResults));
        gammaCM.push_back(get<8>(analysisResults));
        eGammaCM.push_back(get<9>(analysisResults));
        gammaRate.push_back(get<10>(analysisResults));
        eGammaRate.push_back(get<11>(analysisResults));
        muCharge.push_back(get<12>(analysisResults));
        eMuCharge.push_back(get<13>(analysisResults));
        gammaCharge.push_back(get<14>(analysisResults));
        eGammaCharge.push_back(get<15>(analysisResults));    
    }

    //Read CAEN parameters from .root CAEN file
    caenParam = hvReader(scanFolder, run, file_count, verbose);
    hvEff = get<0>(caenParam);
    eHVeff = get<1>(caenParam);
    iMon = get<2>(caenParam);
    eImon = get<3>(caenParam);
    hvMon = get<4>(caenParam);
    eHVmon = get<5>(caenParam);

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
    TGraphErrors *effHv = new TGraphErrors(eff.size(),&hvEff[0],&eff[0],NULL,&eEff[0]);
    effHv->SetMarkerStyle(8);
    effHv->SetMarkerSize(1.4);
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
    TGraphErrors *currHv = new TGraphErrors(iMon.size(),&hvEff[0],&iMon[0],&eHVeff[0],&eImon[0]);
    currHv->SetMarkerStyle(8);
    currHv->SetMarkerSize(1);
    currHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    currHv->GetYaxis()->SetTitle("I [#muA]");
    currHv->GetYaxis()->SetRangeUser(iMon.front() - 0.1*iMon.front(),iMon.back() + 0.1*iMon.back());
    currHv->Draw("AP");

    //Plot muon CS(HV) curve
    TCanvas *cMuCS = new TCanvas();
    cMuCS->cd();
    TGraphErrors *muCsHv = new TGraphErrors(muCS.size(),&hvEff[0],&muCS[0],NULL,&eMuCS[0]);
    muCsHv->SetMarkerStyle(8);
    muCsHv->SetMarkerSize(1);
    muCsHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    muCsHv->GetYaxis()->SetTitle("Muon Cluster Size [strips]");
    muCsHv->GetYaxis()->SetRangeUser(muCS.front() - 0.1*muCS.front(),muCS.back() + 0.1*muCS.back());
    muCsHv->Draw("AP");

    //Plot muon CM(HV) curve
    TCanvas *cMuCM = new TCanvas();
    cMuCM->cd();
    TGraphErrors *muCmHv = new TGraphErrors(muCM.size(),&hvEff[0],&muCM[0],NULL,&eMuCM[0]);
    muCmHv->SetMarkerStyle(8);
    muCmHv->SetMarkerSize(1);
    muCmHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    muCmHv->GetYaxis()->SetTitle("Muon Cluster Multiplicity");
    muCmHv->GetYaxis()->SetRangeUser(muCM.front() - 0.1*muCM.front(),muCM.back() + 0.1*muCM.back());
    muCmHv->Draw("AP");

    //Plot gamma CS(HV) curve
    TCanvas *cGammaCS = new TCanvas();
    cGammaCS->cd();
    TGraphErrors *gammaCsHv = new TGraphErrors(gammaCS.size(),&hvEff[0],&gammaCS[0],NULL,&eGammaCS[0]);
    gammaCsHv->SetMarkerStyle(8);
    gammaCsHv->SetMarkerSize(1);
    gammaCsHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    gammaCsHv->GetYaxis()->SetTitle("Gamma Cluster Size [strips]");
    gammaCsHv->GetYaxis()->SetRangeUser(gammaCS.front() - 0.1*gammaCS.front(),gammaCS.back() + 0.1*gammaCS.back());
    gammaCsHv->Draw("AP");

    //Plot gamma CM(HV) curve
    TCanvas *cGammaCM = new TCanvas();
    cGammaCM->cd();
    TGraphErrors *gammaCmHv = new TGraphErrors(gammaCM.size(),&hvEff[0],&gammaCM[0],NULL,&eGammaCM[0]);
    gammaCmHv->SetMarkerStyle(8);
    gammaCmHv->SetMarkerSize(1);
    gammaCmHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    gammaCmHv->GetYaxis()->SetTitle("Gamma Cluster Multiplicity");
    gammaCmHv->GetYaxis()->SetRangeUser(gammaCM.front() - 0.1*gammaCM.front(),gammaCM.back() + 0.1*gammaCM.back());
    gammaCmHv->Draw("AP");

    //Plot gamma rate(HV) curve
    TCanvas *cGammaRate = new TCanvas();
    cGammaRate->cd();
    TGraphErrors *gammaRateHv = new TGraphErrors(gammaRate.size(),&hvEff[0],&gammaRate[0],NULL,&eGammaRate[0]);
    gammaRateHv->SetMarkerStyle(8);
    gammaRateHv->SetMarkerSize(1);
    gammaRateHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    gammaRateHv->GetYaxis()->SetTitle("Gamma rate [Hz/cm^{2}]");
    gammaRateHv->GetYaxis()->SetRangeUser(gammaRate.front() - 0.1*gammaRate.front(),gammaRate.back() + 0.1*gammaRate.back());
    gammaRateHv->Draw("AP");

    //Plot muon prompt charge(HV) curve
    TCanvas *cMuCharge = new TCanvas();
    cMuCharge->cd();
    TGraphErrors *muChargeHv = new TGraphErrors(muCharge.size(),&hvEff[0],&muCharge[0],NULL,&eMuCharge[0]);
    muChargeHv->SetMarkerStyle(8);
    muChargeHv->SetMarkerSize(1);
    muChargeHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    muChargeHv->GetYaxis()->SetTitle("Muon prompt charge [pC]");
    muChargeHv->GetYaxis()->SetRangeUser(muCharge.front() - 0.1*muCharge.front(),muCharge.back() + 0.1*muCharge.back());
    muChargeHv->Draw("AP");

    //Plot gamma prompt charge(HV) curve
    TCanvas *cGammaCharge = new TCanvas();
    cGammaCharge->cd();
    TGraphErrors *gammaChargeHv = new TGraphErrors(gammaCharge.size(),&hvEff[0],&gammaCharge[0],NULL,&eGammaCharge[0]);
    gammaChargeHv->SetMarkerStyle(8);
    gammaChargeHv->SetMarkerSize(1);
    gammaChargeHv->GetXaxis()->SetTitle("HV_{eff} [V]");
    gammaChargeHv->GetYaxis()->SetTitle("Gamma prompt charge [pC]");
    gammaChargeHv->GetYaxis()->SetRangeUser(gammaCharge.front() - 0.1*gammaCharge.front(),gammaCharge.back() + 0.1*gammaCharge.back());
    gammaChargeHv->Draw("AP");

    TFile *fAnalysisOut = new TFile(("Analysis_run_"+to_string(run)+".root").c_str(),"RECREATE");
    fAnalysisOut->cd();
    cEff->Write(("Efficiency_run_"+to_string(run)).c_str());
    cCurr->Write(("Current_run_"+to_string(run)).c_str());
    cMuCS->Write(("Muon_CS_run_"+to_string(run)).c_str());
    cMuCM->Write(("Muon_CM_run_"+to_string(run)).c_str());
    cGammaCS->Write(("Gamma_CS_run_"+to_string(run)).c_str());
    cGammaCM->Write(("Gamma_CM_run_"+to_string(run)).c_str());
    cGammaRate->Write(("Gamma_rate_run_"+to_string(run)).c_str());
    cMuCharge->Write(("Muon_prompt_charge_run_"+to_string(run)).c_str());
    cGammaCharge->Write(("Gamma_prompt_charge_run_"+to_string(run)).c_str());

    fAnalysisOut->Close();
}
