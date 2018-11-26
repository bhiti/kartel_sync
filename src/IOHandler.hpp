// #include <TDirectory.h>
// #include <TTree.h>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TString.h"

#include <string>
#include <iostream>

class IOHandler
{
public:
  IOHandler();
  ~IOHandler();

//   void Init(event* Event, char* fname_drs="", char* fname_tel="", char* fname_anchor="", char* fname_results="", char* fname_ana="", char *brname_drs="", vector<TString> lname_drs=0);
  void Init(event* Event, const char* fname_drs="", const char* fname_tel="", const char* fname_anchor="", const char* fname_results="", const char* fname_ana="", const char *brname_drs="");
  
  void GetEvent(int ie, int rel_offset);
  void GetEvent_s(int ie_s, int rel_offset);
  void GetTimeStamp(int ie); 
  Long64_t GetDrsTimestamp(int ie);
  Long64_t GetDrsDelta(int ie);
  Long64_t GetTelTimestamp(int ie);
  Long64_t GetTelDelta(int ie);
  
  Long64_t GetFEI4Timestamp(int ie);
  Long64_t GetFEI4Delta(int ie);
  
  event* _event;
  analyzed_event* _anaevent;
  void FillAnalyzedEvent(int itr[2], float charge[NCH_DRS], float sigw_min[NCH_DRS]);

// private:
  int nch;
  
  TFile *f_drs;
  TTree* t_drs_wf;
  TTree* v_t_drs_wf[NCH_DRS];  
  TTree* t_drs_ts;
  
  TFile *f_tel;
  TTree* t_tel_ts;
  TTree* t_tel_tracks;
  TTree* t_tel_slopeX;
  TTree* t_tel_slopeY;
  
  TFile *f_anchor;
  TTree* t_anchor_hits;
  TTree* t_anchor_ts;
  TTree* t_anchor;
  
  TFile *f_res;   // results
  TFile *f_ana;   // analyzed data
  TTree* t_ana[NCH_DRS];
  TTree* t_ana_tracks;  
  
  // For synchronizing
  int CreateSyncedFile();
  void FillSyncedFile(int ie, int rel_offset);
  void FillSyncedFileInvalid();
  void FillSyncedFileDummy(int ie, int rel_offset);
  int AddTelescopePlane(int n);
  
  std::vector<Int_t> planes;
  std::vector<TDirectory*> inputPlaneDirs;
  std::vector<TTree*> inputHitTrees;
  TTree* dummyEvent;      // to fill invalid events
  TTree* outputEvents;
  std::vector<TDirectory*> outputPlaneDirs;
  std::vector<TTree*> outputHitTrees;
  
//   bool FnameIsEmpty(char* );
};
