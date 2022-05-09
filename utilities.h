#ifndef UTILITIES_H
#define UTILITIES_H

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"

#include <fstream>

struct chamber {
	std::string name = "ALICE-2-0";
	int stripNum = 8, numSamp = 1024;
	float sampFreq = 2.5; //Gs/s
	int strips[8] = {1,2,3,4,5,6,7,8};
};

void treeProducer();

#endif