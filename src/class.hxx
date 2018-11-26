
// #define DRS_N_POINTS 1024
// #define NCH_DRS 8

// struct window
// {
//   window(){};
//   window(float xmin, float xmax, float ymin, float ymax) : x1(xmin), x2(xmax), y1(ymin), y2(ymax) {};
//   bool Contains(float x, float y) {if(x>=x1 && x<=x2 && y>=y1 && y<=y2 ) return true; else return false; };
//   void Set(float xmin, float xmax, float ymin, float ymax){x1=xmin, x2=xmax, y1=ymin, y2=ymax; return;};
//   float x1, x2, y1, y2;  
// };
// 
// class waveform
// {
// public:
//   waveform();
//   ~waveform();
// 
//   Float_t wf[DRS_N_POINTS];
//   Float_t GetMin();
//   Float_t GetMin(int first);
//   Float_t GetMin(int first, int last);
//   Float_t GetMax();
//   Float_t GetMax(int first);
//   Float_t GetMax(int first, int last);
//   Float_t GetMaxFiltered(int bin_start, int bin_end, float thr);
//   Float_t GetMaxFiltered(int skip, float thr, int window);
//   Float_t GetMaxMultiFiltered(int bin_start, int bin_end, float thr);
//   Float_t GetMaxPos();
//   Float_t GetMaxPos(int skip);
//   Float_t GetMaxPos(int first, int last);
//   Float_t GetMinPos();
//   Float_t GetMinPos(int skip);
//   Float_t GetMinPos(int first, int last);
//   Float_t Integral(int bin_start, int bin_end);
//   Float_t MaxStep(float thr, int first=0, int last=DRS_N_POINTS-1);
//   Float_t MaxStepN(int N, int sign, float thr, int dx, int first=0, int last=DRS_N_POINTS-1);
//   Float_t SearchStep(int sign, float thr, int width, int step, int first=0, int last=DRS_N_POINTS-1);
//   Float_t GetDelta(int bin_start, int bin_end, int N=10);
//   TGraph *Draw();
// };

// class event
// {
// public:
//   event();
//   ~event();
//   waveform wf[NCH_DRS];
//   TGraph *gr[NCH_DRS];
//   Float_t _charge[NCH_DRS];
//   
//   ULong64_t ts_drs;
//   ULong64_t ts_tel;
//   ULong64_t ts_anchor;
//   Int_t tel_trgoffset;
//   Long64_t delta_drs;
//   Long64_t delta_tel;
//   float ratio;
//   bool invalid;       // telescope invalid events
//   
//   Double_t x[1000];
//   Double_t y[1000];
//   int n;
//   
//   ULong64_t ts_tel_s;   // For anchor module
//   Int_t x_a[1000];
//   Int_t y_a[1000];
//   int n_a;
//   Double_t slopeX[1000];
//   Double_t slopeY[1000];
//   
//   void Draw(int ch, char* option="");
// //   float GetCharge(int skip, float thr, int window); 
//   float GetCharge(int ch, int bin_start, int bin_end, int mode=0);
//   float GetChargeCSA(int ch, int bin_start, int bin_end, int mode=0);
//   int SearchTracks(window w);
//   int SearchOneTrack(window w);
//   int TrackInRegion(int it, window w);
//   int NTracksInRegion(window w);
// 
//   bool VetoPresent(int ch, float thr, int first=0, int last=NCH_DRS);
// };
  
class run
{
public:
  run();
  run(char* fname_drs, char *brname_drs, char* lname_drs, int ch);
  run(char* fname_drs, char* fname_tel, char *brname_drs, char *brname_tel, vector<TString> lname_drs);
  run(char* fname_drs, char* fname_tel, char* fname_tel_s, char* fname_anchor, char *brname_drs, char *brname_tel, char* lname_drs, int ch);
  ~run();
  TFile *f_drs;
  int nch;
  
  TTree* t_drs_wf;
  TTree* v_t_drs_wf[NCH_DRS];
  
  TTree* t_drs_ts;
  TFile *f_tel;
  TTree* t_tel_ts;
  TTree* t_tel_trgoffset;
  TTree* t_tel_invalid;
  TTree* t_tel_tracks;
  TTree* t_tel_slopeX;
  TTree* t_tel_slopeY;
  TFile *f_anchor;
  TTree* t_anchor_hits;
  TTree* t_anchor_ts;
  
  waveform avwf;
  event Event;

  std::vector<float> v_posmax;
  std::vector<float> v_max;
  
  int N;
  int rel_offset;
  float ratio_mean, ratio_tolerance;
  float ratio_low, ratio_high;
  float sigw_min[NCH_DRS], sigw_max[NCH_DRS];
  float noisew_min[NCH_DRS], noisew_max[NCH_DRS];
  
  void GetTimeStamp(int ie); 
  Long64_t GetDrsTimestamp(int ie);
  Long64_t GetDrsDelta(int ie);
  Long64_t GetTelTimestamp(int ie);
  Long64_t GetTelDelta(int ie);
  float GetRatio(int i0, int i1);
  void ProcessTimestamp(int ie, int verbose);
  bool ProcessInvalid(int ie, int verbose);
  
  void FindSignalWindow(int ch, int nsamples, int delta_nbins=50); 
  void FindSignalWindow(int ch, int nsamples, int delta_low, int delta_high); 
  void SetSignalWindow(int ch, float min, float max);
  void SetNoiseWindow(int ch, float min, float max);
  float FindThreshold(int ch, int nsamples=1000, float sigma=4., int mode=0);
  float GetEfficiency(int ch, float thr, int mode=0);
//   float GetEfficiency(int it, float thr);
  float GetNoise(int ch, int mode=0);
  TH2S *heff_tracks[NCH_DRS];
  TH2S *heff_hits[NCH_DRS];
  TH2F *heff_efficiency[NCH_DRS];
  window ROI;
  
  void GetEvent(int ie); 
  void GetEventDRS(int ie); 
  TH2S *hTracks;
  TH2S *hTracks_ROI;
  
  void GetEvent_s(int ie_s);  // for file synchronized with fei4
  TH1I *hResX;
  TH1I *hResY;
  TH2I *hRes;
  float GetResidualX(int it);
  float GetResidualY(int it);
  void FillResiduals();
  void FillResiduals(window anchor_ROI);
  window ResCut;
  bool CheckResiduals(int itr[2]);
  
//   DUTCorrelation dutcorr;
  
  //   void SetResidualCut(float xmin, float xmax, float ymin, float ymax);
  
private:
  int BestOffset(int ie, int nsteps, int verbose);
  bool RatioInRange(float ratio);
  float TryPairRatio(int i0, int i1);
  bool in_pair;
  
};

