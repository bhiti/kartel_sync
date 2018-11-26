// #include "waveform.cpp"
// #include "event.cpp"
// #include "DUTCorrelation.cxx"
#include <iostream>
#include "TH1I.h"
#include "TH2I.h"

class run
{
public:
  run();
  run(const char* fname_drs, const char* fname_tel, const char *fname_anchor, const char *fname_synced);
  run(char* fname_tel, char *fname_anchor);
  
  run(char* fname_drs, char *brname_drs, char* lname_drs, int ch);
  run(char* fname_drs, char* fname_tel, char *brname_drs, char *brname_tel, vector<TString> lname_drs);
  run(char* fname_drs, char* fname_tel, char* fname_anchor, char* fname_results, char *brname_drs, char *brname_tel, vector<TString> lname_drs);
  run(char* fname_drs, char* fname_tel, char* fname_anchor, char* fname_results, char* fname_ana, char *brname_drs, char *brname_tel, vector<TString> lname_drs);
  ~run();
  
  TH1I *hDeltaTimeStamp;
  TH2I *hDeltaVsInvalid;
  TH1I *hRatio;
  
  IOHandler *IOHand;
  Histograms *histos;
  int nch;
  waveform* avwf;
  event* Event;
  analyzed_event AnaEvent;
  
  void FillAnalyzedEvent(int itr[2], float charge[NCH_DRS]);
  int CreateSyncedFile();
  void FillSyncedFile(int ie);
  void FillSyncedFileInvalid();
  void FillSyncedFileDummy(int ie);
  int AddTelescopePlane(int n);
  void WriteAll();

  int N;
  int rel_offset;     // global offset for sync
  int abs_offset;
  float ratio_mean, ratio_tolerance;
  float ratio_low, ratio_high;
  float SetRatio(float mean, float tolerance);
  
  float GetRatio(int i0, int i1, int mode=0);
  int ProcessTimestamp(int ie, int mode, int verbose);
  int ProcessOffsets(int ie);
  
  float sigw_min[NCH_DRS], sigw_max[NCH_DRS];
  float noisew_min[NCH_DRS], noisew_max[NCH_DRS];
  void FindSignalWindow(int ch, int nsamples, int delta_nbins=50); 
  void FindSignalWindow(int ch, int nsamples, int delta_low, int delta_high, int pulse_polarity=+1); 
  void SetSignalWindow(int ch, float min, float max);
  void SetNoiseWindow(int ch, float min, float max);
  float FindThreshold(int ch, int nsamples=1000, float sigma=4., int mode=0);
  float GetEfficiency(int ch, float thr, int mode=0, int pulse_polarity=+1);
  float GetEfficiencyClustered(float *thr, int mode=0, int pulse_polarity=+1);
  float GetNoise(int ch, int mode=0);

  window ROI;
  void SetROI(float xmin, float xmax, float ymin, float ymax, int nbinsx, int nbinsy=-1);
  
  void GetEvent(int ie); 
  void GetEventDRS(int ie); 
  void GetTimeStamp(int ie);
  float zDUT[NCH_DRS];
  float zAnchor;
  
  Long64_t GetDrsTimestamp(int ie);
  Long64_t GetDrsDelta(int ie);
  Long64_t GetTelTimestamp(int ie);
  Long64_t GetTelDelta(int ie);
  float GetTelRatio(int i0, int i1, int mode);
  
  void GetEvent_s(int ie_s);  // for file synchronized with fei4
  float GetResidualX(int it, int jt);
  float GetResidualY(int it, int jt);
  void FillResiduals();
  void FillResiduals(window anchor_ROI);
  window ResCut;
  bool CheckResiduals(int itr[2], double maxchi2=100);
  bool CheckResiduals(int ie, int itr[2], double maxchi2=100);
  
  DUTCorrelation dutcorr;
  
  int bigOffset=0;
  
  
  //   void SetResidualCut(float xmin, float xmax, float ymin, float ymax);
private:
  std::pair<int, int> BestOffset(int ie, int nsteps, int mode, int verbose);
//int BestOffset(int ie, int nsteps, int mode, int verbose);
//   int* BestOffset_abs(int ie, int nsteps, int mode, int verbose);
  bool RatioInRange(float ratio);
  float TryPairRatio(int i0, int i1);  
  std::pair<int, int> offsets;    // local offsets <rel_offset, abs_offset> for sync
  
};

