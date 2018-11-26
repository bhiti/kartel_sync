// #include "TROOT.h"
// #include "TFile.h"
// #include "TDirectory.h"
// #include "TTree.h"
// #include "TBranch.h"
// #include "TH2D.h"
// 
// #include <string>
// #include <iostream>
// 
// #include "progressBar.h"
// 
// #define MAXHITS 10000
// 
// int telescopeToProteus(const char *inPath, const char *outFilePath){

//   std::cout << " Input file: " << inPath << std::endl;
//   std::cout << "Output file: " << outFilePath << std::endl;
  std::vector<Int_t> planes;
  planes.push_back(1);
  planes.push_back(3);
  planes.push_back(4);

  //#################################
  //load input file and link branches
  //#################################
  TFile * inputFile = f_tel;
//   TFile * inputFile = NULL;
//   inputFile = new TFile(inPath,"READ");
  if(!inputFile){
    std::cout << __PRETTY_FUNCTION__ << " :: Input file could not be opened!" << std::endl;
    return -1;
  }

//   TTree * inputEvents = (TTree *)inputFile->Get("Event");
  TTree * inputEvents = t_tel_ts;
////   if(!inputEvents){
//     std::cout << "Could not open input event tree!" << std::endl;
//     return -1;
//   }
//   ULong64_t frameNumber;
  if(inputEvents->SetBranchAddress("FrameNumber", &_event->tel_frameNumber) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
//   ULong64_t timestamp;
//   if(inputEvents->SetBranchAddress("TimeStamp", &_event->tel_timestamp) != 0){
//     std::cout << "Branch not found, stopping script..." << std::endl;
//     return -1;
//   }

//   ULong64_t triggerTime;
//   if(t_tel_ts->SetBranchAddress("TriggerOffset", &_event->tel_triggerOffset) != 0){
//     std::cout << "Branch not found, stopping script..." << std::endl;
//     return -1;
//   }
//   Int_t triggerInfo;
  if(inputEvents->SetBranchAddress("TriggerInfo", &_event->tel_triggerInfo) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
//   Int_t triggerOffset;
  if(inputEvents->SetBranchAddress("TriggerOffset", &_event->tel_triggerOffset) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
//   UInt_t triggerBCID;
//   if(inputEvents->SetBranchAddress("TriggerBCID", &triggerBCID) != 0){
//     std::cout << "Branch not found, stopping script..." << std::endl;
//     return -1;
//   }
//   UInt_t triggerL1ID;
//   if(inputEvents->SetBranchAddress("TriggerL1ID", &triggerL1ID) != 0){
//     std::cout << "Branch not found, stopping script..." << std::endl;
//     return -1;
//   }
//   Bool_t invalid = 0;
  if(inputEvents->SetBranchAddress("Invalid", &_event->tel_invalid) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }

  //set up plane dirs and hit trees to read in
  std::vector<TDirectory *> inputPlaneDirs;
  std::vector<TTree *> inputHitTrees;
  for(UInt_t i = 0; i < planes.size(); i++){
    TDirectory * planeDir = new TDirectory();
    std::string name = "Plane" + std::to_string(planes[i]);
    f_tel->GetObject(name.c_str(), planeDir);
    if(!planeDir){
      std::cout << __PRETTY_FUNCTION__ << ": Plane" << planes[i] << " directory not found!" << std::endl;
      return -1;
    }
    inputPlaneDirs.push_back(planeDir);

    TTree * hitsTree = new TTree();
    inputPlaneDirs[i]->GetObject("Hits",hitsTree);
    if(!hitsTree){
      std::cout << __PRETTY_FUNCTION__ << ": Hits tree not found!" << std::endl;
      return -1;
    }
    inputHitTrees.push_back(hitsTree);
  }

  Int_t maxHits = 10000;
  //define variables to link to input file hit tree branches
  Int_t numHits;
  Double_t hitValue[maxHits];
  Double_t hitTiming[maxHits];
  Int_t hitPixX[maxHits];
  Int_t hitPixY[maxHits];
  Int_t hitInCluster[maxHits];

  //##################
  //set up output file
  //##################
  TFile * outputFile = f_res;
//   TFile * outputFile = NULL;
//   outputFile = new TFile(outFilePath,"RECREATE");
//   if(!outputFile){
//     std::cout << __PRETTY_FUNCTION__ << " :: Output file could not be opened!" << std::endl;
//     return -1;
//   }
  
  TTree * outputEvents = new TTree("Event","Event Information");
  if(!outputEvents){
    std::cout << "Could not open output event tree!" << std::endl;
    return -1;
  }
  outputEvents->Branch("FrameNumber", &_event->tel_frameNumber, "FrameNumber/l");
  outputEvents->Branch("TimeStamp", &_event->ts_tel, "TimeStamp/l");
//   outputEvents->Branch("TriggerTime", &triggerTime, "TriggerTime/l");
  outputEvents->Branch("TriggerInfo", &_event->tel_triggerInfo, "TriggerInfo/I");
  outputEvents->Branch("TriggerOffset", &_event->tel_triggerOffset, "TriggerOffset/I");
//   outputEvents->Branch("TriggerL1ID", &triggerL1ID, "TriggerL1ID/i");
//   outputEvents->Branch("TriggerBCID", &triggerBCID, "TriggerBCID/i");
  outputEvents->Branch("Invalid", &_event->tel_invalid, "Invalid/O");

  //link hit tree branches to variables
  std::vector<TDirectory *> outputPlaneDirs;
  std::vector<TTree *> outputHitTrees;
  for(UInt_t k = 0; k < planes.size(); k++){
    std::string name = "Plane" + std::to_string(k);
    outputPlaneDirs.push_back(outputFile->mkdir(name.c_str()));
    outputPlaneDirs[k]->cd();
    TTree * hits = new TTree("Hits","Hits");
    if(!hits){
      std::cout << "Could not create hits tree for plane " << name << "!" << std::endl;
      return -1;
    }
    outputHitTrees.push_back(hits);
    hits->Branch("NHits", &numHits, "NHits/I"); //event->size()
    hits->Branch("Value", hitValue, "Value[NHits]/D"); //(*event)[i].ToT
    hits->Branch("Timing", hitTiming, "Timing[NHits]/D"); //(*event)[i].time1_56ns
    hits->Branch("PixX", hitPixX, "PixX[NHits]/I"); //(*event)[i].X
    hits->Branch("PixY", hitPixY, "PixY[NHits]/I"); //(*event)[i].Y
    hits->Branch("HitInCluster", hitInCluster, "HitInCluster[NHits]/I");
  }
  //#########################################################################
  //loop through events of input file and write selected planes to ouput file
  //#########################################################################
  Int_t nEntries = inputEvents->GetEntries();
  for(Int_t i = 0; i < nEntries; i++){
    loadBar(i, nEntries);
    inputEvents->GetEntry(i);
    for(UInt_t j = 0; j < planes.size(); j++){
      inputHitTrees[j]->SetBranchAddress("NHits", &numHits);
      inputHitTrees[j]->SetBranchAddress("PixX", hitPixX);
      inputHitTrees[j]->SetBranchAddress("PixY", hitPixY);
      Int_t tmp_hitTiming[maxHits];
      inputHitTrees[j]->SetBranchAddress("Timing", tmp_hitTiming);
      Int_t tmp_hitValue[maxHits];
      inputHitTrees[j]->SetBranchAddress("Value", tmp_hitValue);
//       inputHitTrees[j]->SetBranchAddress("HitInCluster", hitInCluster);
      inputHitTrees[j]->SetBranchAddress("InCluster", hitInCluster);
      inputHitTrees[j]->GetEntry(i);
      for(Int_t i = 0; i < numHits; i++){
        hitTiming[i] = (Double_t)tmp_hitTiming[i];
        hitValue[i] = (Double_t)tmp_hitValue[i];
      }
      outputHitTrees[j]->Fill();
    }
    outputEvents->Fill();
  }
//   outputFile->Write();
  outputFile->Close();
  inputFile->Close();

  return 0;
}
