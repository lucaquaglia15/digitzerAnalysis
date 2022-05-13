#include "utilities.h"

void treeProducer (int numFiles, string folderPath, string trigFile, string waveFile) {

	gSystem->cd(folderPath.c_str());

	fstream dataFile;
	string fileName;

	for (int i = 0; i< numFiles; i++) {
		fileName = waveFile + to_string(i) + ".txt";
		dataFile.open(fileName.c_str());

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			cout << "Opening file: " << fileName << endl;
		} 

		else { //File opening not succesful, close the program
			cout << "Error in opening wave file" << endl;
			exit(0);
		}

		dataFile.close();

	}
}

void countsConverter() {
	//Open root file produced earlier

}