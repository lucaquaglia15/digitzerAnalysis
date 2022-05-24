#include "utilities.h"

//Produce .root file
/*void treeProducer (int numFiles, string folderPath, string trigFile, string waveFile) {

	gSystem->cd(folderPath.c_str());

	TFile *fData = new TFile("digitizerData.root","RECREATE"); //root file to save data
	TTree *dataTree = new TTree("dataTree","Data from digitizer"); //tree
	dataTree->SetEntries(-1);
	TBranch *branch[8];

	//int evNum = 0;
	//double count[8][1024]; //-> this works
	double count[1024];
	//double count[8][501][1024];
	//long double count[8][513024];
	//double count[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	string fileName; //Name of the input file

	for (int i = 0; i < numFiles; i++) {
		int evNum = 0;
		int j = 0; 
		int lineCounter = 1;
		int address = 0;
		fstream dataFile; //For txt input file -> to stream data
		string branchName = "counts_strip_" + to_string(i);
		//TBranch *branch = new TBranch (dataTree, branchName.c_str(),&count,"count/D",32000);
		//TBranch *branch = dataTree->Branch(branchName.c_str(),&count[i][j]); // this works
		//TBranch *branch = dataTree->Branch(branchName.c_str(),count,"count/D");
		//branch[i] = dataTree->Branch(branchName.c_str(),&count[i][j],"count[1024]/D"); //This kinda works
		//branch[i] = dataTree->Branch(branchName.c_str(),&count[i],"count[1024]/D");
		cout << "Im here before branch" << endl;
		//branch[i] = dataTree->Branch(branchName.c_str(),&count[i][evNum][j],"count[1024]/D");
		branch[i] = dataTree->Branch(branchName.c_str(),count,"count[1024]/D");
		//branch[i] = dataTree->Branch(branchName.c_str(),&count[i][j]);
		//branch[i] = dataTree->Branch(branchName.c_str(),&count[i][j],"count[8][513024]/D");

		cout << "Im here after branch" << endl;
		//cout << "Index: " << i << "\t Branch address: " << &count[i][j] << endl;

		//cout << "Index: " << i << "\t Branch address: " << &count[i] << endl;
		//dataTree->Branch(branchName.c_str(),&count[i][0]);
		//dataTree->Branch(branchName.c_str(),count+i);
		//cout << &count[i][j] << endl;
		fileName = waveFile + to_string(i) + ".txt";
		//fileName = "wave_0.txt";
		dataFile.open(fileName.c_str());
		if (fileName == "wave_0.txt") dataTree->Branch("evNum",&evNum,"evNum/I");

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			cout << "Opening file: " << fileName << endl;
			string tp = ""; 
			while(getline(dataFile,tp)) {

				if (fileName == "wave_0.txt" && lineCounter % 1032 == 4) {
					//cout << tp.substr(tp.find(": ") + 1) << endl; 
					evNum = stoi(tp.substr(tp.find(": ") + 1));
					//cout << endl << endl << "+++++" << evNum << "+++++" << endl << endl;
					//dataTree->Fill();
				//cout << tp << endl; //Event number
				}

				if ((lineCounter % 1032 >= 9 && lineCounter % 1032 <= 1031) || (lineCounter % 1032 == 0)) {
					//cout << tp << endl;
					//count[i][j] = stod(tp); //This kinds works
					//count[i][evNum][j] = stod(tp);
					//count[i][j] = stod(tp);
					//count[i] = stod(tp);
					count[j] = stod(tp);
					//count[i][j] = stod(tp);
					j++;
					//branch[i]->Fill();
					//dataTree->Fill();
					//cout << i << "\t" << j << endl;
					//cout << count[i][j] << endl;
					//count[j] = stod(tp);
					//cout << count[j] << endl;
					//j++;
					//cout << j << endl;
					//cout << count + i << endl;
					//dataTree->Fill();
					//branch->Fill();
					
					if (j == 1024) {
						//cout << "Filled array!" << endl;
						//count[j] = stod(tp);
						//cout << i << "\t" << j << endl;
						//dataTree->Fill();
						branch[i]->Fill(); //this kinda works
						cout << evNum << endl;
						branch[i]->FlushBaskets();
						//dataTree->Fill();
						//cout << j << endl;
						j = 0;
						address++;
					} //uncomment if needed

					//cout << count[i] << endl;
					//dataTree->Fill();
					//cout << count << endl;
				}
				lineCounter++; //Increase line counter
				//dataTree->Fill(); //Fill tree branches
				//branch->Fill();
			}
			//dataFile.close();
			//dataTree->Fill();
		}

		else { //File opening not succesful, close the program
			cout << "Error in opening wave file" << endl;
			exit(0);
		}

		//branch[i]->Fill();
		//branch[i]->FlushBaskets();

		dataFile.close();
		//cout << j << endl;
		//branch[i]->Fill();
		//dataTree->Fill();
		
		//dataFile.close();
		//fData->cd();
		//dataTree->Write();
	}
	
	//Save tree in root file
	//dataTree->Fill();
	fData->cd();
	dataTree->Write();
	fData->Close();
}*/

//struct strip {
	//double counts[513024]; //This works
//	double count;
//};

void treeProducer(int numFiles, string folderPath, string trigFile, string waveFile) {

	string fileName;

	gSystem->cd(folderPath.c_str());

	TFile *fData = new TFile("digitizerData.root","RECREATE"); //root file to save data
	TTree *dataTree = new TTree("dataTree","Data from digitizer"); //tree

	vector<vector<string>> data;
	vector<string> strip;
	vector<double> counts[8];

	for (int i = 0; i < 8; i++) {
		string branchName = "counts_strip_" + to_string(i);
		dataTree->Branch(branchName.c_str(), &counts[i]);
		cout << &(counts[i]) << endl;
	}

	for (int i = 0; i < 8; i++) {

		fstream dataFile; //For txt input file -> to stream data
		fileName = waveFile + to_string(i) + ".txt";
		dataFile.open(fileName.c_str());

		if (dataFile.is_open()) { //Open file and parse useful information to put in .root file
			cout << "Opening file: " << fileName << endl;
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

	cout << data[0].size() << endl;

	for (unsigned int i = 0; i < data[0].size(); i++) {
		if ((i % 1032 >= 8 && i % 1032 <= 1031)) { //ADC data
			for (int j = 0; j < 8; j++) {
				counts[j].push_back(stod(data[j][i]));
				//cout << data[j][i] << endl;
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

//Convert ADC counts to mV
void analyzer() {
	//Open root file and get the tree inside
	TFile *f = new TFile("digitizerData.root","READ");
	TTree *t = (TTree*)f->Get("dataTree");


	//Convert the counts to mV

	//Calculate backgournd and check if any strip is over 5 sigma

	//Plot the eight strips



}