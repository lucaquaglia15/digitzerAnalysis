#include "utilities.h"

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
	string tp = ""; //new position
	for (int i = 0; i < 8; i++) {

		fstream dataFile; //For txt input file -> to stream data
		fileName = waveFile + to_string(i) + ".txt";
		dataFile.open(fileName.c_str());

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			if (verbose) cout << "Opening file: " << fileName << endl;
			//string tp = ""; //old position

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
	}*/ //WORKING VERSION!!!

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
	}

	fData->cd();
	dataTree->Write();
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

double analyzer(int numFiles, string folderPath, int hv, bool verbose, chamber1G rpc, int muonMin, int muonMax, int noiseMin, int noiseMax) { //1 Gs/s
//double analyzer(int numFiles, string folderPath, int hv, bool verbose, chamber rpc, int muonMin, int muonMax, int noiseMin, int noiseMax) { // 5 Gs/s

	cout << endl << "Analyzing HV value: " << hv+1 << endl;

	int samples = rpc.numSamp; //number of samples
	float freq = rpc.sampFreq; //sampling frequency
	float amplitue = rpc.vpp;
	float convFactor = amplitue/4096; //From ADC to mV
	float meanNoise = 0., muonPeak = 0., sq_sum, stdev; //
	bool hasHit = false;
	int eventsWitHit = 0, eventsWtihoutHit = 0, muonCounter = 0, rateCounter = 0;

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
	cout << "Number of entries in the tree: " << nEntries << endl; //Debug

	if (verbose) cout << rpc.numSamp << endl;

	//vector<vector<double>*> data = 0.0;
	vector<float> *data[8] = {0,0,0,0,0,0,0,0};
	//vector<float> datamV[8] = {0,0,0,0,0,0,0,0};
	vector<float> datamV[8];
	//vector<float> *data = 0; //data in ADC counts

	//vector<vector<float>> datamV; //data in mV
	vector<float> stripmV;

	//Compute time vector (for x axis)
	float time = 0;
	vector<float> times;
	for (int i = 0; i < samples; i++) {
		if (i == 0) times.push_back(0);
		else {
			time = time + (samples/freq)/samples;
			times.push_back(time);
		}
	}

	if (verbose) for (unsigned int i = 0; i < times.size(); i++) cout << times.at(i) << endl;
	
	//t->SetBranchAddress("counts_strip_0",&data);

	//Set Branch addresses for the 8 strips
	for (int i = 0; i < numFiles; i++) {
		string branchName = "counts_strip_" + to_string(i);
		t->SetBranchAddress(branchName.c_str(),&data[i]);
		if (verbose) cout << "Setting branch address for: " << branchName << endl;
	}

	//cout << data[0]->size() << endl;

	TCanvas *cStrip = new TCanvas();
	cStrip->cd();
	cStrip->Divide(1,numFiles);
	TGraph *gStrip[numFiles];

	//Analyze the data
	// 1) Get tree entries
	for (int i = 0; i < nEntries; i++) { //cycle on all the events
		//if (i == 0) {
		if (verbose) cout << "Entry: " << i << endl;
			
			t->GetEntry(i); //Get i-th entry
			//cout << data[0]->size() << endl;

			//cout << *min_element(data[0]->begin(), data[0]->end()) << endl;

			if (*min_element(data[0]->begin(), data[0]->end()) < 1000.) { //It's a muon - 1 Gs/s, if 5 Gs/s comment this line
				muonCounter++;

				for (int j = 0; j < numFiles; j++) { //Channel 0 excluded because it's the trigger only
					//Conversion to mV 
					/*for (unsigned int k = 0; k < data[0]->size(); k++) {
						//datamV[i].push_back(data[j]->at(k) * ); //CONVERSION HERE
					}*/
					if (verbose) cout << "Size " + to_string(i) << data[j]->size() << endl;

					//Check if the event is rate or efficiency

					gStrip[j] = new TGraph(data[j]->size(),&times[0],&data[j]->at(0));
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

					if (j == 0) continue;

					meanNoise = accumulate(&data[j]->at(noiseMin), &data[j]->at(noiseMax), 0.0) / (noiseMax-noiseMin);

					sq_sum = inner_product(&data[j]->at(noiseMin), &data[j]->at(noiseMax), &data[j]->at(noiseMin),0.0);
					stdev = std::sqrt(sq_sum/(noiseMax-noiseMin) - meanNoise*meanNoise); //std dev in noise window
					
					//muonPeak = *min_element(&data[j]->at(muonMin), &data[j]->at(muonMax)); //Min in muon window - Negative pulse polarity
					muonPeak = *max_element(&data[j]->at(muonMin), &data[j]->at(muonMax)); //Max in muon window - Positive pule polarity

					//cout << "Event: \t" << i << "\t strip: \t" << j+1 << "\t muon pulse: \t" << muonPeak << "\t noise (with 5 sigma): \t" << meanNoise - 5*stdev << endl; //Negative pulse polarity
					//cout << "Event: \t" << i << "\t strip: \t" << j+1 << "\t muon pulse: \t" << muonPeak << "\t noise (with 5 sigma): \t" << meanNoise + 5*stdev << endl; //Positive pulse polarity

					//if (hasHit == false && muonPeak <= (meanNoise - 5*stdev)) { //Negative pulse polarity
					if (hasHit == false && muonPeak >= (meanNoise + 5*stdev)) { //Positive pulse polarity
						hasHit = true;
					}

					if (verbose) {
						cout << "noise: " << stdev << endl;
						cout << "muon: " << muonPeak << endl;
						//for (unsigned t = 0; t < data[j]->size(); t++) cout << data[j]->at(t) << endl;
					}

					meanNoise = 0.;
					muonPeak = 0.;
					sq_sum = 0.;
					stdev = 0.;			
				}
			}

			else { //It's a rate event, for the moment we ignore it
				rateCounter++;
				continue; 
			}
		//}

	//else break; //To test only on the first event

		if (hasHit == true) {
			eventsWitHit++;
			hasHit = false;

			fout->cd();
			dirEff->cd();
			cStrip->Write(("Event_"+to_string(i)).c_str());
		}

		else if (hasHit == false) {
			eventsWtihoutHit++;
			fout->cd();
			dirNonEff->cd();
			cStrip->Write(("Event_"+to_string(i)).c_str());
		}
	}

	fout->Close();

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
