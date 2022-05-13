#include "utilities.h"

//Produce .root file
void treeProducer (int numFiles, string folderPath, string trigFile, string waveFile) {

	gSystem->cd(folderPath.c_str());

	TFile *fData = new TFile("digitizerData.root","RECREATE"); //root file to save data
	TTree *dataTree = new TTree("dataTree","Data from digitizer"); //tree

	fstream dataFile; //For txt input file -> to stream data
	string fileName; //Name of the input file

	for (int i = 0; i< numFiles; i++) {
		fileName = waveFile + to_string(i) + ".txt";
		dataFile.open(fileName.c_str());
		int lineCounter = 1;

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			cout << "Opening file: " << fileName << endl;
			string tp; 
			while(getline(dataFile,tp)) {
				if (lineCounter % 1032 == 4) cout << tp << endl; //Event number
				lineCounter++; //Increase line counter
			}
		} 

		else { //File opening not succesful, close the program
			cout << "Error in opening wave file" << endl;
			exit(0);
		}

		dataFile.close();
	
	}
}

//Convert ADC counts to mV
void countsConverter() {
	//Open root file produced earlier

}