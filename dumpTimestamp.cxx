#include <iostream>
#include <fstream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>

void dumpTimeStamp(const char* filename)
{
	std::stringstream ss; ss << filename << "__dumpedTimeStamp.txt";
    ofstream dumpFile;
    dumpFile.open (ss.str().c_str());

	TFile* file = new TFile(filename);
	TTree* tree = (TTree*) file->Get("Event");
	ULong64_t timeStamp;
	tree->SetBranchAddress("TimeStamp",&timeStamp);
	const Long64_t Nentries = tree->GetEntries();

// 	for (Long64_t n = 0; n < Nentries; n++) {
	for (Long64_t n = 0; n < Nentries; n++) {
		tree->GetEntry(n);
		dumpFile << timeStamp << "\n";
		if(n%10000 == 0) cout << n << endl;
	}
}
