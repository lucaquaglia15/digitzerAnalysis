#include "utilities.h"

//Produce .root file
void treeProducer (int numFiles, string folderPath, string trigFile, string waveFile) {

	gSystem->cd(folderPath.c_str());

	TFile *fData = new TFile("digitizerData.root","RECREATE"); //root file to save data
	TTree *dataTree = new TTree("dataTree","Data from digitizer"); //tree

	int evNum = 0;
	double count[8][1024];
	// = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	dataTree->Branch("evNum",&evNum,"evNum/I");

	string fileName; //Name of the input file

	for (int i = 0; i < numFiles; i++) {
		int j = 0, lineCounter = 1;
		fstream dataFile; //For txt input file -> to stream data
		string branchName = "counts_strip_" + to_string(i);
		//TBranch *branch = new TBranch (dataTree, branchName.c_str(),&count,"count/D",32000);
		TBranch *branch = dataTree->Branch(branchName.c_str(),&count[i][j]);
		//dataTree->Branch(branchName.c_str(),&count[i]);
		//dataTree->Branch(branchName.c_str(),count+i);
		cout << &count[i][j] << endl;
		fileName = waveFile + to_string(i) + ".txt";
		//fileName = "wave_0.txt";
		dataFile.open(fileName.c_str());

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			cout << "Opening file: " << fileName << endl;
			string tp; 
			while(getline(dataFile,tp)) {
				if (lineCounter % 1032 == 4) {
					//cout << tp << endl; //Event number
				}
				else if ((lineCounter % 1032 >= 9 && lineCounter % 1032 <= 1031) || (lineCounter % 1032 == 0)) {
					//cout << tp << endl;
					count[i][j] = stod(tp);
					j++;
					//cout << count + i << endl;
					//dataTree->Fill();
					//branch->Fill();
					if (j == 1023) {
						dataTree->Fill();
						j = 0;
					}
					//cout << count[i] << endl;
					//cout << count << endl;
				}
				lineCounter++; //Increase line counter
				//dataTree->Fill(); //Fill tree branches
				//branch->Fill();
			}
			dataFile.close();
		} 

		else { //File opening not succesful, close the program
			cout << "Error in opening wave file" << endl;
			exit(0);
		}

		//dataFile.close();
		//fData->cd();
		//dataTree->Write();
	}
	
	//Save tree in root file
	fData->cd();
	dataTree->Write();
	fData->Close();
}

//Convert ADC counts to mV
void countsConverter() {
	//Open root file produced earlier

}