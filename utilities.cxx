#include "utilities.h"
#include "PeakFinder.h"
#include "SGSmooth.hpp"

//--------------------//
//                    //
//    Produce tree    //
//                    //
//--------------------//

void treeProducer(int numFiles, string folderPath, string trigFile, string waveFile, int hv, bool verbose) {

	string fileName;

	gSystem->cd(folderPath.c_str());
	if (verbose) cout << folderPath << endl;

	TFile *fData = new TFile(("digitizerData_HV_" + to_string(hv+1) + ".root").c_str(),"RECREATE"); //root file to save data
	TTree *dataTree = new TTree("dataTree","Data from digitizer"); //tree

	vector<vector<string>> data;
	vector<string> strip;
	vector<double> counts[8];

	for (int i = 0; i < 8; i++) {
		string branchName = "counts_strip_" + to_string(i);
		dataTree->Branch(branchName.c_str(), &counts[i]);
		if (verbose) cout << &(counts[i]) << endl;
	}
	//string tp = ""; //new position
	for (int i = 0; i < 8; i++) {

		fstream dataFile; //For txt input file -> to stream data
		fileName = waveFile + to_string(i) + ".txt";
		dataFile.open(fileName.c_str());

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			if (verbose) cout << "Opening file: " << fileName << endl;
			string tp = ""; //old position, working

			while(getline(dataFile,tp)) {
				strip.push_back(tp);
			}
			data.push_back(strip);
			strip.clear();
		}

		else { //File opening not succesful, close the program
			cout << "Error in opening wave file" << endl;
			exit(0);
		}
		dataFile.close();
	}

	/*for (unsigned int i = 0; i < data.size(); i++) { //Debug printouts
 		cout << i << endl; //8 = number of readout strips
		cout << data[i].size() << endl; //number of lines in each file
	}*/

	int counter = 0;

	if (verbose) cout << data[0].size() << endl;

	for (unsigned int i = 0; i < data[0].size(); i++) {
		if ((i % 1032 >= 8 && i % 1032 <= 1031)) { //ADC data
			for (int j = 0; j < 8; j++) {
				counts[j].push_back(stod(data[j][i]));
				counter++;
				if (counter == 1024*8) { //event finished
					dataTree->Fill();
					counter = 0;
					for (unsigned k = 0; k < 8; k++) counts[k].clear();
				}
			}
		}
	} //<- WORKING VERSION!!! ->

	/*for (unsigned int i = 0; i < data[0].size(); i++) {
		if ((i % 1032 >= 8 && i % 1032 <= 1031)) { //ADC data
			for (int j = 0; j < 8; j++) {
				counts[j].push_back(stod(data[j][i]));
				counter++;
				if (counter == 1024*8) { //event finished
					dataTree->Fill();
					counter = 0;
					for (unsigned k = 0; k < 8; k++) counts[k].clear();
				}
			}
		}
	}*/

	fData->cd();
	dataTree->Write();
	//dataTree->Write("", TObject::kOverwrite); //In theory this saves the data without creating the same tree twice (i.e. dataTree10, dataTree11)
	fData->Close();
}

//--------------------//
//                    //
//Count # of HV values//
//                    //
//--------------------//

int hvCounter(string folderPath, string ext, bool verbose) {
	const char* entry; 
	const char* filename;
    int file_count = 0;
    TString str;
      
    char* dir = gSystem->ExpandPathName(folderPath.c_str());
    void* dirp = gSystem->OpenDirectory(dir);
    
    while((entry=gSystem->GetDirEntry(dirp))) { //if the file ends with ".root" -> increase the file count by 1
    	str = entry;
      	if(str.EndsWith(ext.c_str())) {
      		file_count++;
    	}
    }
	return file_count;
}

//--------------------//
//                    //
//    Analyze data    //
//                    //
//--------------------//

//double analyzer(int numFiles, string folderPath, int hv, bool verbose, chamber1G rpc, int muonMin, int muonMax, int noiseMin, int noiseMax) { //1 Gs/s
double analyzer(vector<float> times, int numFiles, string folderPath, int hv, bool verbose, chamber1G rpc, int muonMin, int muonMax, int noiseMin, int noiseMax) { //1 Gs/s
//double analyzer(int numFiles, string folderPath, int hv, bool verbose, chamber rpc, int muonMin, int muonMax, int noiseMin, int noiseMax) { // 5 Gs/s

	cout << endl << "Analyzing HV value: " << hv+1 << endl;
	bool graphMakeMuon = true, graphMakeGamma = true;

	if (graphMakeMuon == true) cout << "Graph making is enabled for muons" << endl;
	else cout << "Graph making is disabled for muons" << endl;

	if (graphMakeGamma == true) cout << "Graph making is enabled for gammas" << endl;
	else cout << "Graph making is disabled for gammas" << endl; 

	/////////////
	// General //
	/////////////
	int samples = rpc.numSamp; //number of samples
	int stripArea = rpc.stripArea;
	float freq = rpc.sampFreq; //sampling frequency
	float amplitue = rpc.vpp;
	float convFactor = amplitue/4096; //From ADC to mV
	double threshold = 0.;//This will be used once we have the conversion between FEERIC threshold and actual mV
	double clusTime = 15.; //Clustering time for muons and gammas, can be changed manually

	///////////
	// Muons //
	///////////
	double meanNoise = 0., muonPeak = 0., sq_sum = 0, stdev = 0, rms = 0;
	int eventsWitHit = 0, eventsWtihoutHit = 0, muonCounter = 0;
	bool muonEff[7] = {false,false,false,false,false,false,false}; //strip by strip efficiency
	bool isMuon = false, hasHitMuon = false;
	vector<float> muonPeakValue;
	vector<int> muonPeakPos;
	//TH1I *muonPeakTimeAllStrips = new TH1I("muonPeakTimeAllStrips","muonPeakTimeAllStrips",300,0.5,300.5);
	//TH1I *muonPeakTimePerStrip[7];
	//for (int i = 0; i < 7; i++) muonPeakTimePerStrip[i] = new TH1I(("muonStrips"+to_string(i+1)).c_str(),("muonStrips"+to_string(i+1)).c_str(),300,0.5,300.5);
	int clusSizeMuon = 1, clusMultMuon = 0;
	bool isClusMuon = false;
	vector<int> clustersMuon;
	TH1F *muonProfile = new TH1F("muonProfile","muonProfile",7,0.5,7.5); //muon strip profile
	TH1I *muonClusterMult = new TH1I("muonClusterMult","muonClusterMult",8,0.5,8.5); //muon cluster multiplicity
	TH1F *muonClusterSize = new TH1F("muonClusterSize","muonClusterSize",8,0.5,8.5); //muon cluster size

	////////////
	// Gammas //
	////////////
	double meanNoiseGamma = 0., gammaPeak = 0., sq_sum_gamma = 0, stdevGamma = 0, rmsGamma = 0;
	int eventsWithGammaHit = 0, eventsWtihoutGammaHit = 0, gammaCounter = 0;
	bool gammaEff[7] = {false,false,false,false,false,false,false};
	bool isGamma = false, hasHitGamma = false;
	vector<float> gammaPeakValue;
	vector<int> gammaPeakPos;
	//TH1I *gammaPeakTimeAllStrips = new TH1I("muonPeakTimeAllStrips","muonPeakTimeAllStrips",300,0.5,300.5);
	//TH1I *gammaPeakTimePerStrip[7];
	//for (int i = 0; i < 7; i++) gammaPeakTimePerStrip[i] = new TH1I(("gammaStrips"+to_string(i+1)).c_str(),("gammaStrips"+to_string(i+1)).c_str(),300,0.5,300.5);
	int clusSizeGamma = 1, clusMultGamma = 0;
	bool isClusGamma = false;
	vector<int> clustersGamma;
	TH1F *gammaProfile = new TH1F("gammaProfile","gammaProfile",7,0.5,7.5); //gamma strip profile
	TH1F *gammaRateProfile = new TH1F("gammaRateProfile","gammaRateProfile",7,0.5,7.5); //gamma rate profile
	TH1I *gammaClusterMult = new TH1I("gammaClusterMult","gammaClusterMult",8,0.5,8.5); //gamma cluster multiplicity
	TH1F *gammaClusterSize = new TH1F("gammaClusterSize","gammaClusterSize",8,0.5,8.5); //gamma cluster size
	
	vector<float> tempMinValue, tempMaxValue;

	//cout << samples << "\t" << freq << "\t" << amplitue << endl;

	//Open root file and get the tree inside
	gSystem->cd(folderPath.c_str());

	TFile *fout = new TFile(("HV_" + to_string(hv+1) + "_ANALYSIS.root").c_str(),"RECREATE");
	TDirectory *dirEff = fout->mkdir("efficientEvents"); //Efficient muon events display
	TDirectory *dirNonEff = fout->mkdir("nonEfficientEvents"); //Non efficient muon events display
	TDirectory *dirEffGamma = fout->mkdir("gammaEfficientEvents"); //Efficient gamma events display
	TDirectory *dirNonEffGamma = fout->mkdir("gammaNonEfficientEvents");//Non efficient gamma events display

	TFile *f = new TFile(("digitizerData_HV_" + to_string(hv+1) + ".root").c_str(),"READ");
	TTree *t = (TTree*)f->Get("dataTree");

	int nEntries = t->GetEntries();
	cout << "Number of entries in the tree: " << nEntries << endl;

	if (verbose) cout << rpc.numSamp << endl;

	vector<float> *data[8] = {0,0,0,0,0,0,0,0};
	vector<float> datamV[8];
	vector<float> gammaValues[8];
	for (int i = 0; i < 8; i++) {
		gammaValues[i].reserve(1024*12502);
		datamV[i].reserve(1024);
	}
	
	vector<float> stripmV;

	if (verbose) for (unsigned int i = 0; i < times.size(); i++) cout << times.at(i) << endl;

	//Set Branch addresses for the 8 strips
	for (int i = 0; i < numFiles; i++) {
		string branchName = "counts_strip_" + to_string(i);
		t->SetBranchAddress(branchName.c_str(),&data[i]);
		if (verbose) cout << "Setting branch address for: " << branchName << endl;
	}

	//cout << data[0]->size() << endl;

	TCanvas *cStrip = new TCanvas(); //mV(time) for muons, eight per event
	cStrip->cd();
	cStrip->Divide(1,numFiles);
	TCanvas *cStripGamma = new TCanvas(); //mV(time) for gammas, eight per event
	cStripGamma->cd();
	cStripGamma->Divide(1,numFiles);
	TCanvas *cStripGammaTot = new TCanvas(); //mV(time) for gammas, eight in total, joined data from each event
	cStripGammaTot->cd();
	cStripGammaTot->Divide(1,numFiles);

	TGraph *gStrip[numFiles];
	TGraph *gStripGamma[numFiles];
	TGraph *gStripGammaTot[numFiles];
	TLine *noiseLine[numFiles];
	TMarker *peak[numFiles];

	//Analyze the data
	// 1) Get tree entries
	for (int i = 0; i < nEntries; i++) { //cycle on all the events

		if (verbose) cout << "Entry: " << i << endl;
			
		t->GetEntry(i); //Get i-th entry
		
		if (verbose) {
			cout << data[0]->size() << endl;
			cout << data[0] << endl;
		}

		////////////////
		// Muon event //
		////////////////

		//cout << *min_element(data[0]->begin(), data[0]->end()) << endl;
		//Check if the event is muon or gamma by looking at the trigger (on channel 0)
		if (*min_element(data[0]->begin(), data[0]->end()) < 1000.) { //It's a muon - 1 Gs/s, if 5 Gs/s comment this line
			muonCounter++;
			isMuon = true;

			for (int j = 0; j < numFiles; j++) { //Channel 0 excluded because it's the trigger only
				
				//Conversion to mV 
				//Calculate mean oscillations, subtract to each point and convert to mV
				meanNoise = accumulate(&data[j]->at(noiseMin), &data[j]->at(noiseMax), 0.0) / (noiseMax-noiseMin);
				
				for (unsigned int k = 0; k < data[0]->size(); k++) {
					datamV[j].push_back((data[j]->at(k)-meanNoise)*convFactor); //CONVERSION HERE
				}

				//STD DEV (ADC counts)
				//sq_sum = inner_product(&data[j]->at(noiseMin), &data[j]->at(noiseMax), &data[j]->at(noiseMin),0.0);
				//stdev = std::sqrt(sq_sum/(noiseMax-noiseMin) - meanNoise*meanNoise); //std dev in noise window
			
				//RMS (mV) -> Remember that inner_production function is [....) meaning that last element is not included in the calculation
				sq_sum = inner_product(&datamV[j].at(noiseMin), &datamV[j].at(noiseMax), &datamV[j].at(noiseMin),0.0);
				rms = std::sqrt(sq_sum/(noiseMax - noiseMin));

				if (verbose) cout << "Size " + to_string(i) << "\t" << data[j]->size() << endl;

				if (graphMakeMuon) {
					//gStrip[j] = new TGraph(data[j]->size(),&times[0],&data[j]->at(0)); //ADC counts
					gStrip[j] = new TGraph(samples,&times[0],&datamV[j][0]); //mV
					//gStrip[j]->GetXaxis()->SetTitle("Time [ns]");
					gStrip[j]->GetXaxis()->SetLabelSize(0.1);
					gStrip[j]->GetXaxis()->SetTickSize(0.01);
					gStrip[j]->GetYaxis()->SetNdivisions(510);
					//gStrip[j]->GetYaxis()->SetTitle("Amplitude [a.u.]");
					gStrip[j]->GetYaxis()->SetLabelSize(0.1);
					gStrip[j]->GetYaxis()->SetTickSize(0.01);
					gStrip[j]->GetYaxis()->SetNdivisions(2);
					gStrip[j]->SetMarkerColor(kRed);
					gStrip[j]->SetLineColor(kRed);
					gStrip[j]->SetTitle("");
					cStrip->cd(j+1);
					gStrip[j]->Draw("APL");
					//noiseLine[j] = new TLine(times.front(),meanNoise + 5*stdev,times.back(),meanNoise + 5*stdev); //ADC counts
					noiseLine[j] = new TLine(0,5*rms,1024,5*rms); //mV
					noiseLine[j]->SetLineColor(kBlack);
					noiseLine[j]->Draw("SAME");
				}
				
				if (j == 0) { //Trigger event, simply draw it and continue						
					datamV[j].clear();
					continue;
				}

				if (verbose) {
					//cout << "Event: \t" << i << "\t strip: \t" << j+1 << "\t muon pulse: \t" << muonPeak << "\t noise (with 5 sigma): \t" << meanNoise - 5*stdev << endl; //Negative pulse polarity
					cout << "Event: \t" << i << "\t strip: \t" << j+1 << "\t muon pulse: \t" << muonPeak << "\t noise (with 5 sigma): \t" << meanNoise + 5*stdev << endl; //Positive pulse polarity
				}

				//Value corresponding to muon peak
				//muonPeak = *min_element(&data[j]->at(muonMin), &data[j]->at(muonMax)); //Min in muon window - Negative pulse polarity - ADC counts
				//muonPeak = *max_element(&data[j]->at(muonMin), &data[j]->at(muonMax)); //Max in muon window - Positive pule polarity - ADC counts
				//muonPeak = *min_element(&datamV[j]->at(muonMin), &datamV[j]->at(muonMax)); //Min in muon window - Negative pulse polarity - mV
				muonPeak = *max_element(&datamV[j].at(muonMin), &datamV[j].at(muonMax)); //Max in muon window - Positive pule polarity - mV

				//Was the chamber efficient?
				//if (muonEff[j-1] == false && muonPeak <= (meanNoise - 5*stdev)) { //Negative pulse polarity - ADC counts
				//if (muonEff[j-1] == false && muonPeak >= (meanNoise + 5*stdev)) { //Positive pulse polarity - ADC counts
				//if (muonEff[j-1] == false && muonPeak <= threshold) { //FEERIC equivalent negative threshold - mV
				//if (muonEff[j-1] == false && muonPeak >= threshold) { //FEERIC equivalent positive threshold - mV
				//if (muonEff[j-1] == false && muonPeak <= 5*rms) { //Negative pulse polarity - mV
				if (muonEff[j-1] == false && muonPeak >= 5*rms) { //Positive pulse polarity - mV
					muonEff[j-1] = true;
					muonPeakPos.push_back(find(datamV[j].begin(),datamV[j].end(),muonPeak)-datamV[j].begin());
					muonPeakValue.push_back(muonPeak);
					if (verbose) cout << "Evento " << i << " strip " << j << " peak pos " << find(datamV[j].begin(),datamV[j].end(),muonPeak) - datamV[j].begin()  << endl;
					//muonPeakTimeAllStrips->Fill(find(datamV[j].begin(),datamV[j].end(),muonPeak)-datamV[j].begin());
					//muonPeakTimePerStrip[j-1]->Fill(find(datamV[j].begin(),datamV[j].end(),muonPeak)-datamV[j].begin());
					if (hasHitMuon == false) hasHitMuon = true;
				}

				else { //not efficient
					muonPeakPos.push_back(0); //push 0 if there is no muon, to have always 7 elements in the vector
					muonPeakValue.push_back(0); //push 0 if there is no muon, to have always 7 elements in the vector
				}

				//Min and max values of each event in order to rescale y axis in graph
				//tempMinValue.push_back(*min_element(&data[j]->front(), &data[j]->back()));
				//tempMaxValue.push_back(*max_element(&data[j]->front(), &data[j]->back()));

				if (verbose) {
					cout << "MIN: " << *min_element(&data[j]->front(), &data[j]->back()) << endl;
					cout << "MAX: " << *max_element(&data[j]->front(), &data[j]->back()) << endl;
					cout << "size " << tempMinValue.size() << "\t" << tempMaxValue.size() << endl;
				}
				meanNoise = 0.;
				muonPeak = 0.;
				sq_sum = 0.;
				stdev = 0.;	
				datamV[j].clear(); //clear values in mV for next event, done strip by strip		
			}

			//float minValue = *min_element(&tempMinValue.front(), &tempMinValue.back());
			//float maxValue = *max_element(&tempMaxValue.front(), &tempMaxValue.back());
			
			if (verbose) {
				cout << "REAL min: " << *min_element(&tempMinValue.front(), &tempMinValue.back()) << endl;
				cout << "REAL max: " << *max_element(&tempMaxValue.front(), &tempMaxValue.back()) << endl;
			}

			/*for (int i = 1; i < numFiles; i++) { //Graph cosmetics + memory clearance
				//gStrip[i]->GetYaxis()->SetRangeUser(minValue-0.05*minValue,maxValue+0.05*maxValue);
				gStrip[i]->GetYaxis()->SetRangeUser(minValue,maxValue);
				cStrip->cd(i+1);
				gStrip[i]->Draw("APL");
				noiseLine[i]->Draw("SAME");
			}
			tempMinValue.clear();
			tempMaxValue.clear();*/
		} //End of muon event

		/////////////////
		// Gamma event //
		/////////////////

		else { //It's a rate event
			//continue;
			gammaCounter++;
			isGamma = true;

			for (int k = 0; k < numFiles; k++) {

				meanNoiseGamma = 0;
				sq_sum_gamma = 0;
				stdevGamma = 0;

				meanNoiseGamma = accumulate(&data[k]->at(noiseMin), &data[k]->at(noiseMax), 0.0) / (noiseMax-noiseMin);
				
				for (unsigned int j = 0; j < samples; j++) {
					datamV[k].push_back((data[k]->at(j)-meanNoiseGamma)*convFactor); //CONVERSION HERE
				}

				//STD DEV (ADC counts)
				//sq_sum_gamma = inner_product(&data[k]->at(noiseMin), &data[k]->at(noiseMax), &data[k]->at(noiseMin),0.0);
				//stdevGamma = std::sqrt(sq_sum_gamma/(noiseMax-noiseMin) - meanNoiseGamma*meanNoiseGamma); //std dev in noise window
				
				//RMS (mV) -> Remember that inner_production function is [....) meaning that last element is not included in the calculation
				sq_sum_gamma = inner_product(&datamV[k].at(noiseMin), &datamV[k].at(noiseMax), &datamV[k].at(noiseMin),0.0);
				rmsGamma = std::sqrt(sq_sum_gamma/(noiseMax-noiseMin)); //RMS in noise window;

				if (graphMakeGamma) {
					//gStripGamma[k] = new TGraph(samples,&times[0],&data[k]->at(0)); //ADC counts
					gStripGamma[k] = new TGraph(samples,&times[0],&datamV[k][0]); //mV
					gStripGamma[k]->GetXaxis()->SetLabelSize(0.1);
					gStripGamma[k]->GetXaxis()->SetTickSize(0.01);
					gStripGamma[k]->GetYaxis()->SetNdivisions(510);
					gStripGamma[k]->GetYaxis()->SetLabelSize(0.1);
					gStripGamma[k]->GetYaxis()->SetTickSize(0.01);
					gStripGamma[k]->GetYaxis()->SetNdivisions(2);
					//gStripGamma[k]->GetYaxis()->SetRangeUser(meanNoiseGamma-20*stdevGamma,meanNoiseGamma+20*stdevGamma);
					gStripGamma[k]->SetMarkerColor(kRed);
					gStripGamma[k]->SetLineColor(kRed);
					gStripGamma[k]->SetTitle("");
					cStripGamma->cd(k+1);
					gStripGamma[k]->Draw("APL");
					//noiseLine[k] = new TLine(0,meanNoiseGamma + 5*stdevGamma,1024,meanNoiseGamma + 5*stdevGamma); //ADC counts
					noiseLine[k] = new TLine(0,5*rmsGamma,1024,5*rmsGamma); //mV
					noiseLine[k]->SetLineColor(kBlack);
					noiseLine[k]->Draw("SAME");
				}
				

				if (k == 0) { //Strip 0 is the int+ext coincidence, draw it but nothing will appear because there is no coincidence in a gamma event
					datamV[k].clear();
					continue;
				}

				//Value corresponding to gamma peak
				gammaPeak = *max_element(&datamV[k].at(0), &datamV[k].at(980)); //Max in window - Positive pule polarity - mV
				//gammaPeak = *min_element(&datamV[k].at(0), &datamV[k].at(980)); //Min in window - Negative pulse polarity - mV
				//gammaPeak = *min_element(&data[k]->at(0), &data[k]->at(980)); //Min in window - Negative pulse polarity - ADC counts
				//gammaPeak = *max_element(&data[k]->at(0), &data[k]->at(980)); //Max in window - Positive pule polarity - ADC counts

				//Did the chamber see a gamma?
				//if (gammaEff[k-1] == false && gammaPeak <= (meanNoiseGamma - 5*stdevGamma)) { //negative polarity - ADC counts
				//if (gammaEff[k-1] == false && gammaPeak <= 5*rms) { //Negative pulse polarity - mV
				if (gammaEff[k-1] == false && gammaPeak >= 5*rms) { //Positive pulse polarity - mV
				//if (gammaEff[k-1] == false && gammaPeak <= threshold) { //FEERIC equivalent negative threshold - mV
				//if (gammaEff[k-1] == false && gammaPeak >= threshold) { //FEERIC equivalent positive threshold - mV
				//if (gammaEff[k-1] == false && gammaPeak >= (meanNoiseGamma + 5*stdevGamma)) { //positive polarity - ADC counts
					//cout << "Strip " << k << " max elem " << *max_element(data[k]->begin(), data[k]->end()) << " 5 sigma " << 5*stdevGamma << endl;
					if (hasHitGamma == false) hasHitGamma = true;
					gammaEff[k-1] = true;
					//peak[k] = new TMarker(max_element(&data[k]->at(0), &data[k]->at(980))-&data[k]->at(0),*max_element(&data[k]->at(0), &data[k]->at(980)),8);
					peak[k] = new TMarker(max_element(&datamV[k].at(0), &datamV[k].at(980))-&datamV[k].at(0),*max_element(&datamV[k].at(0), &datamV[k].at(980)),8);
					peak[k]->SetMarkerColor(kBlue);
					peak[k]->SetMarkerSize(1);
					peak[k]->Draw("SAME");
					gammaPeakPos.push_back(find(datamV[k].begin(),datamV[k].end(),gammaPeak)-datamV[k].begin());
					gammaPeakValue.push_back(gammaPeak);
					//gammaPeakTimeAllStrips->Fill(find(datamV[j].begin(),datamV[j].end(),gammaPeak)-datamV[j].begin());
					//gammaPeakTimePerStrip[j-1]->Fill(find(datamV[j].begin(),datamV[j].end(),gammaPeak)-datamV[j].begin());
				}
				
				else { //not efficient
	 				gammaPeakPos.push_back(0); //push 0 if there is no gamma, to have always 7 elements in the vector
					gammaPeakValue.push_back(0); //push 0 if there is no gamma, to have always 7 elements in the vector
				}

				datamV[k].clear();

				//Create long vector to join alla gamma events for visulization purposes
				//cout << "Size gammaValues before: " << gammaValues[k].size() << endl;
				//if (gammaValues[k].size() == 0) gammaValues[k] = *data[k]; //uncomment to create long vector
				//else gammaValues[k].insert(gammaValues[k].end(), data[k]->begin(), data[k]->end()); //uncomment to create long vector
				//cout << "Size gammaValues after: " << gammaValues[k].size() << endl;
			}
			if (verbose) cout << "Number of gamma events: " << gammaCounter << endl;
		} //End of gamma event

		if (isMuon == true) { //it's a muon event
			isMuon = false;

			if (hasHitMuon == true) { //at least one strip saw a muon
				hasHitMuon = false;
				eventsWitHit++;
				if (graphMakeMuon) {
					fout->cd();
					dirEff->cd();
					cStrip->Write(("Event_"+to_string(i)).c_str());
				}
				for (int s = 0; s < 7; s++) {
					if (muonEff[s] == true) {
						muonEff[s] = false; //reset counters for next event
						muonProfile->Fill(s+1,1); //fill muon strip profile strip by strip
					}
				}

				//Muon clustering calculation
				for (int cs = 0; cs < 7; cs++) {
					//cs 0  1  2  3  4  5  6 -> few examples for logic building
					//   -------------------
					//   x  x  x  x  x  x  x
					//   1, 1, 0, 1, 0, 0, 1
					//   0, 0, 1, 0, 0, 1, 1
					//   1, 0, 1, 0, 0, 1, 0
					//   0, 0, 0, 0, 0, 0, 0
					//   1, 0, 0, 0, 0, 0, 0

					if (muonPeakPos[cs] != 0) { //element is different from 0
						if (cs == 6 && isClusMuon == true && clusMultMuon == 0) { //There is a single big cluster (i.e. 0, 0, 0, 1, 1, 1, 1)
							clustersMuon.push_back(clusSizeMuon);
							clusMultMuon++;
							isClusMuon = false;
							break;
						}

						else if (cs == 6 && isClusMuon == true && clusMultMuon > 0) { //There is a single big cluster (i.e. 0, 0, 1, 1, 0, 1, 1)
							clustersMuon.push_back(clusSizeMuon);
							clusMultMuon++;
							isClusMuon = false;
							break;
						}

						else if (cs == 6 && isClusMuon == false) { //For example this case (0, 0, 0, 1, 1, 0, 1).
							//CM does not matter since it's a new cluster, formed only by the last element
							clusMultMuon++;
							clusSizeMuon = 1;
							clustersMuon.push_back(clusSizeMuon);
							break;
						}

						if (muonPeakPos[cs+1] != 0) { //following element is != 0
							//Time difference is < 15 ns?
							if (TMath::Abs(muonPeakPos[cs+1] - muonPeakPos[cs]) <= clusTime) { //yes
								isClusMuon = true;
								clusSizeMuon++;
								continue;
							}

							else { //no
								clusMultMuon++;
								clustersMuon.push_back(clusSizeMuon);
								clusSizeMuon = 1;
								isClusMuon = false;
								continue;
							}

						}

						else if (muonPeakPos[cs+1] == 0) { //following element is = 0 -> cluster is done
							clusMultMuon++;
							clustersMuon.push_back(clusSizeMuon);
							clusSizeMuon = 1;
							isClusMuon = false;
						}
					}

					else if (muonPeakPos[cs] == 0) {
						continue;
					}
				}
				//Calculate average cluster size
 				float avgMuonClusSize = accumulate(clustersMuon.begin(),clustersMuon.end(),0.0)/clustersMuon.size();
				
				//Print out debugs
				/*cout << "Event # " << i << endl;
				cout << "CM " << clusMultMuon << endl;
				cout << "Elements in cs vector " << clusters.size() << endl;
				for (unsigned int pp = 0; pp < 7; pp++) {
					cout << muonPeakPos[pp] << "\t";
				}
				cout << endl;
				for (unsigned int pp = 0; pp < clustersMuon.size(); pp++) {
					cout << clustersMuon[pp] << "\t";
				}
				cout << endl;
				cout << "avg clus size " << avgMuonClusSize << endl;*/

				muonClusterSize->Fill(avgMuonClusSize,clusMultMuon); //Fill cluster size histo (weight given by multiplicity)
				muonClusterMult->Fill(clusMultMuon); //Fill cluster multiplicity histo
				
				//Reset all values for next event
				clusSizeMuon = 1; 
				clusMultMuon = 0;
				isClusMuon = false;
				clustersMuon.clear();
			}

			else if (hasHitMuon == false) { //no strip saw a muon
				eventsWtihoutHit++;
				if (graphMakeMuon) { 
					fout->cd();
					dirNonEff->cd();
					cStrip->Write(("Event_"+to_string(i)).c_str());
				}
			}
		} //end of saving muon plots

		else if (isGamma == true) { //it's a gamma event
			isGamma = false;

			if(hasHitGamma == true) {
				hasHitGamma = false;
				eventsWithGammaHit++;
				if (graphMakeGamma) { 
					fout->cd();
					dirEffGamma->cd();
					cStripGamma->Write(("Event_"+to_string(i)).c_str());
				}
				for (int s = 0; s < 7; s++) {
					if (gammaEff[s] == true) {
						gammaEff[s] = false; //reset counters for next event	
						gammaProfile->Fill(s+1,1);
					}
				}

				//Gamma clustering calculation
				for (int cs = 0; cs < 7; cs++) {
					if (gammaPeakPos[cs] != 0) { //element is different from 0
						if (cs == 6 && isClusGamma == true && clusMultGamma == 0) { //There is a single big cluster (i.e. 0, 0, 0, 1, 1, 1, 1)
							clustersGamma.push_back(clusSizeGamma);
							clusMultGamma++;
							isClusGamma = false;
							break;
						}

						else if (cs == 6 && isClusGamma == true && clusMultGamma > 0) { //There is a single big cluster (i.e. 0, 0, 1, 1, 0, 1, 1)
							clustersGamma.push_back(clusSizeGamma);
							clusMultGamma++;
							isClusGamma = false;
							break;
						}

						else if (cs == 6 && isClusGamma == false) { //For example this case (0, 0, 0, 1, 1, 0, 1).
							//CM does not matter since it's a new cluster, formed only by the last element
							clusMultGamma++;
							clusSizeGamma = 1;
							clustersGamma.push_back(clusSizeGamma);
							break;
						}

						if (gammaPeakPos[cs+1] != 0) { //following element is != 0
							//Time difference is < 15 ns?
							if (TMath::Abs(gammaPeakPos[cs+1] - gammaPeakPos[cs]) <= clusTime) { //yes
								isClusGamma = true;
								clusSizeGamma++;
								continue;
							}

							else { //no
								clusMultGamma++;
								clustersGamma.push_back(clusSizeGamma);
								clusSizeGamma = 1;
								isClusGamma = false;
								continue;
							}

						}

						else if (gammaPeakPos[cs+1] == 0) { //following element is = 0 -> cluster is done
							clusMultGamma++;
							clustersGamma.push_back(clusSizeGamma);
							clusSizeGamma = 1;
							isClusGamma = false;
						}
					}

					else if (gammaPeakPos[cs] == 0) {
						continue;
					}
				}
				//Calculate average cluster size
 				float avgGammaClusSize = accumulate(clustersGamma.begin(),clustersGamma.end(),0.0)/clustersGamma.size();
				
				//Print out debugs
				/*cout << "Event # " << i << endl;
				cout << "CM " << clusMultGamma << endl;
				cout << "Elements in cs vector " << clustersGamma.size() << endl;
				for (unsigned int pp = 0; pp < 7; pp++) {
					cout << gammaPeakPos[pp] << "\t";
				}
				cout << endl;
				for (unsigned int pp = 0; pp < clustersGamma.size(); pp++) {
					cout << clustersGamma[pp] << "\t";
				}
				cout << endl;
				cout << "avg clus size " << avgGammaClusSize << endl;*/

				gammaClusterSize->Fill(avgGammaClusSize,clusMultGamma); //Fill cluster size histo (weight given by multiplicity)
				gammaClusterMult->Fill(clusMultGamma); //Fill cluster multiplicity histo
				
				//Reset all values for next event
				clusSizeGamma = 1; 
				clusMultGamma = 0;
				isClusGamma = false;
				clustersGamma.clear();
			}

			else if (hasHitGamma == false) {
				eventsWtihoutGammaHit++;
				if (graphMakeGamma) {
					fout->cd();
					dirNonEffGamma->cd();
					cStripGamma->Write(("Event_"+to_string(i)).c_str());
				}
			}
		} //end of saving gamma plots

		muonPeakPos.clear();
		muonPeakValue.clear();
		gammaPeakPos.clear();
		gammaPeakValue.clear();
	}

	//Join all gamma events together in a single long vector and count the peaks

	/*vector<float> testValues, smoothData;
	vector <float> avgNoise, avg;
	vector<int> peakIndex, noisePeakIndex, smoothPeakIndex; //peakIndex = indeces of gamma peaks, noisePeakIndex = indeces of common noise peaks
	testValues.reserve(1024*12502);
	for (unsigned int i = 0; i < gammaValues[0].size(); i++) {
		testValues.push_back((int)i);
	}

	if (hv == 9) { //last hv point, just to test rate measurement and smoothing
		
		for (int i = 0; i < numFiles; i++) {
			gStripGammaTot[i] = new TGraph(gammaValues[i].size(),&testValues[0],&gammaValues[i].at(0));
			gStripGammaTot[i]->GetXaxis()->SetLabelSize(0.1);
			gStripGammaTot[i]->GetXaxis()->SetTickSize(0.01);
			gStripGammaTot[i]->GetYaxis()->SetNdivisions(510);
			gStripGammaTot[i]->GetYaxis()->SetLabelSize(0.1);
			gStripGammaTot[i]->GetYaxis()->SetTickSize(0.01);
			gStripGammaTot[i]->GetYaxis()->SetNdivisions(2);
			gStripGammaTot[i]->SetMarkerColor(kBlack);
			gStripGammaTot[i]->SetLineColor(kBlack);
			gStripGammaTot[i]->SetTitle("");
			cStripGammaTot->cd(i+1);
			gStripGammaTot[i]->Draw("APL");
			if (verbose) cout << "size cumulative " << i << "\t" << gammaValues[i].size() << endl;
			
			if (i == 0 || i == 3) {
   				if (i == 0) {
   					cout << "Looking for common noise..." << endl;
   					cout << "I found " << PeakFinder::findPeaks(gammaValues[i],noisePeakIndex,true,5) << " noise peaks" << endl;
   				}
   				else {
   					cout << "I found " << PeakFinder::findPeaks2(gammaValues[i],avgNoise,avg,peakIndex,true,5) << " gamma peaks" << endl; //old version
   					cout << "Noise size " << noisePeakIndex.size() << " peak size " << peakIndex.size() << endl;
   					cout << "Trying to smooth out the data for strip: " << i << "..." << endl;

   					smoothData = sg_smooth(gammaValues[i],5,3);
   					cout << "Original size " << gammaValues[i].size() << " new size " << smoothData.size() << endl;

   					cout << "Smoothing done" << endl;
   					cout << "Finding peaks in smoothed data..." << endl;
   					cout << "I found " << PeakFinder::findPeaks(smoothData,smoothPeakIndex,true,5) << " peaks in smothed data (common noise NOT subracted)" << endl;;
   					//cout << "I found " << PeakFinder::findPeaksSmoothed(smoothData,avgNoise,avg,smoothPeakIndex,true,5) << " peaks in smothed data (common noise NOT subracted)" << endl;;

   					for (unsigned int kk = 0; kk < peakIndex.size(); kk++) {
   						//cout << "sono qui" << endl;
   						for (unsigned int zz = 0; zz < noisePeakIndex.size(); zz++) {
   							//cout << "sono qui 2" << endl;
   							if (noisePeakIndex[zz] != peakIndex[kk]) {
   								continue;
   							}

   							else if (noisePeakIndex[zz] == peakIndex[kk]) {
   								peakIndex.erase(peakIndex.begin() + kk); //eliminate from gamma rate vector
   								//cout << "same" << endl;
   								break;
   							}
   						}
   					}
   					cout << "New gamma size " << peakIndex.size() << endl;
   					//vector<TLine*> testRate;
   					TLine *testRate;
   					for (unsigned int zzz = 0; zzz < peakIndex.size(); zzz++) {
   						if (verbose) cout << peakIndex[zzz] << endl;
   						testRate = new TLine(peakIndex[zzz],2334.8,peakIndex[zzz],2549.3);
   						testRate->SetLineColor(kRed);
   						cStripGammaTot->cd(4);
   						testRate->Draw("SAME");
   					}

   					TCanvas *cTest = new TCanvas();
   					TMultiGraph *smoothNonSmooth = new TMultiGraph();
   					//cTest->cd();
   					//cTest->Divide(1,2);
   					//cTest->cd(1);
   					//gStripGammaTot[i]->Draw("APL");
   					//cTest->cd(2);
   					TGraph *testSmooth = new TGraph(smoothData.size(),&testValues[0],&smoothData.at(0));
   					testSmooth->GetXaxis()->SetLabelSize(0.1);
					testSmooth->GetXaxis()->SetTickSize(0.01);
					testSmooth->GetYaxis()->SetNdivisions(510);
					testSmooth->GetYaxis()->SetLabelSize(0.1);
					testSmooth->GetYaxis()->SetTickSize(0.01);
					testSmooth->GetYaxis()->SetNdivisions(2);
					testSmooth->SetMarkerColor(kRed);
					testSmooth->SetLineColor(kRed);
					testSmooth->SetTitle("Smoothed");
					//testSmooth->Draw("APL");
					smoothNonSmooth->Add(gStripGammaTot[i]);
					smoothNonSmooth->Add(testSmooth);
					smoothNonSmooth->Draw("APL");
					//gStripGammaTot[i]->Draw("SAME");
   				}
			}
		}
	}*/
	//Save big graph with all gamma events in the same plot
	//fout->cd();
	//cStripGammaTot->Write("gammaAllEvents");	
	//fout->Close();

	//Strip profiles

	//Muon strip profile
	TCanvas *cMuonProfile = new TCanvas();
	cMuonProfile->cd();
	muonProfile->GetXaxis()->SetTitle("Strip number");
	muonProfile->GetYaxis()->SetTitle("Muon counts");
	muonProfile->GetYaxis()->SetRangeUser(0,muonProfile->GetMaximum()+(0.1*muonProfile->GetMaximum()));
	muonProfile->Draw("HISTO");

	//Gamma strip profile
	TCanvas *cGammaProfile = new TCanvas();
	cGammaProfile->cd();
	gammaProfile->GetXaxis()->SetTitle("Strip number");
	gammaProfile->GetYaxis()->SetTitle("Gamma counts");
	gammaProfile->GetYaxis()->SetRangeUser(0,gammaProfile->GetMaximum()+(0.1*gammaProfile->GetMaximum()));
	gammaProfile->Draw("HISTO");

	//Calculate gamma strip rate
	for (int i = 0; i < 7; i++) {
		gammaRateProfile->Fill(i+1,gammaProfile->GetBinContent(i+1)/(1E-9*gammaCounter*stripArea*980));
	}
	TCanvas *cGammaRateProfile = new TCanvas();
	cGammaRateProfile->cd();
	gammaRateProfile->GetXaxis()->SetTitle("Strip number");
	gammaRateProfile->GetYaxis()->SetTitle("Gammma rate [Hz/cm^{2}]");
	gammaRateProfile->GetYaxis()->SetRangeUser(0,gammaRateProfile->GetMaximum()+(0.1*gammaRateProfile->GetMaximum()));
	gammaRateProfile->Draw("HISTO");

	//Muon cluster multiplicity
	TCanvas *cMuonClusterMult = new TCanvas();
	cMuonClusterMult->cd();
	muonClusterMult->GetXaxis()->SetTitle("Muon CM");
	muonClusterMult->GetXaxis()->SetTitle("Counts");
	muonClusterMult->Draw("HISTO");

	//Muon cluster multiplicity
	TCanvas *cMuonClusterSize = new TCanvas();
	cMuonClusterSize->cd();
	muonClusterSize->GetXaxis()->SetTitle("Muon CS [strips]");
	muonClusterSize->GetXaxis()->SetTitle("Counts");
	muonClusterSize->Draw("HISTO");

	//Gamma cluster multiplicity
	TCanvas *cGammaClusterMult = new TCanvas();
	cGammaClusterMult->cd();
	gammaClusterMult->GetXaxis()->SetTitle("Gamma CM");
	gammaClusterMult->GetXaxis()->SetTitle("Counts");
	gammaClusterMult->Draw("HISTO");

	//Gamma cluster size
	TCanvas *cGammaClusterSize = new TCanvas();
	cGammaClusterSize->cd();
	gammaClusterSize->GetXaxis()->SetTitle("Gamma CS [strips]");
	gammaClusterSize->GetXaxis()->SetTitle("Counts");
	gammaClusterSize->Draw("HISTO");

	//Muon peak time (all strips)
	/*TCanvas *cMuonPeakTime = new TCanvas();
	muonPeakTimeAllStrips->Draw("HISTO");

	//Muon peak time (per strip)
	TCanvas *cMuonPeakTimePerStrip = new TCanvas();
	cMuonPeakTimePerStrip->Divide(2,4);
	for (int i = 0; i < 7; i++) {
		cMuonPeakTimePerStrip->cd(i+1);
		muonPeakTimePerStrip[i]->Draw("HISTO");
	}*/

	//Gamma peak time (all strips)
	/*TCanvas *cGammaPeakTime = new TCanvas();
	gammaPeakTimeAllStrips->Draw("HISTO");

	//Gamma peak time (per strip)
	TCanvas *cGammaPeakTimePerStrip = new TCanvas();
	cGammaPeakTimePerStrip->Divide(2,4);
	for (int i = 0; i < 7; i++) {
		cGammaPeakTimePerStrip->cd(i+1);
		gammaPeakTimePerStrip[i]->Draw("HISTO");
	}*/

	fout->cd();
	cMuonProfile->Write(("Muon_strip_profile_hv_"+to_string(hv+1)).c_str());
	cMuonClusterMult->Write(("Muon_clus_mult_hv_"+to_string(hv+1)).c_str());
	cMuonClusterSize->Write(("Muon_clus_size_hv_"+to_string(hv+1)).c_str());
	cGammaProfile->Write(("Gamma_strip_profile_hv_"+to_string(hv+1)).c_str());
	cGammaRateProfile->Write(("Gamma_rate_profile_hv_"+to_string(hv+1)).c_str());	
	cGammaClusterMult->Write(("Gamma_clus_mult_hv_"+to_string(hv+1)).c_str());
	cGammaClusterSize->Write(("Gamma_clus_size_hv_"+to_string(hv+1)).c_str());
	fout->Close();

	//Uncomment if full analysis is launched (i.e. more than 1 hv point)
	//delete muonProfile;
	//delete gammaProfile;
	//delete gammaRateProfile;
	//delete muonClusterMult;
	//delete muonClusterSize;
	//delete gammaClusterMult;
	//delete gammaClusterSize;

	cout << "Muon counter: " << muonCounter << " efficient events: " << eventsWitHit << " non efficient :" << eventsWtihoutHit << " efficiency: " << eventsWitHit/(double)muonCounter << endl;

	return (eventsWitHit/(double)muonCounter)*100;
}

//--------------------//
//                    //
//   Get HV values    //
//                    //
//--------------------//

vector <double> hvReader(string scanFolder, int scan, int file_count, bool verbose) {

	vector<double> hvEff;

	gSystem->cd(scanFolder.c_str());

	for (int i = 0; i < file_count; i++) {
		
			TFile *f = new TFile(("Scan000" + to_string(scan) + "_HV" + to_string(i+1) +"_CAEN.root").c_str(),"READ");
			TH1F *HVeff = (TH1F*)f->Get("HVeff_ALICE-2-0-GAP");
			TH1F *HVmon = (TH1F*)f->Get("HVmon_ALICE-2-0-GAP");
			TH1F *Imon = (TH1F*)f->Get("Imon_ALICE-2-0-GAP");
			hvEff.push_back(HVeff->GetMean());
		
		
		//else hvEff.push_back(9500); //Specific for scan 243, hv point 8 .root file was produced
	}
	return hvEff;
}
