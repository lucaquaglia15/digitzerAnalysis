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

	for (int i = 0; i < 8; i++) {

		fstream dataFile; //For txt input file -> to stream data
		fileName = waveFile + to_string(i) + ".txt";
		dataFile.open(fileName.c_str());

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			if (verbose) cout << "Opening file: " << fileName << endl;
			string tp = ""; 

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

void analyzer(string folderPath, int hv, bool verbose) {
	//Open root file and get the tree inside
	gSystem->cd(folderPath.c_str());
	TFile *f = new TFile(("digitizerData_HV_" + to_string(hv+1) + ".root").c_str(),"READ");
	TTree *t = (TTree*)f->Get("dataTree");

	int nEntries = t->GetEntries();
	if (verbose) cout << "Number of entries in the tree: " << nEntries << endl; //Debug

	//Convert the counts to mV

	//Calculate backgournd and check if any strip is over 5 sigma

	//Plot the eight strips



}