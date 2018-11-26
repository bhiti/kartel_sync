#include "TGraphErrors.h"
#include "TF1.h"
#include "TPad.h"
// #include "waveform.hpp"

class event
{
public:
  event();
  ~event();
  waveform* wf[NCH_DRS];
  Float_t _charge[NCH_DRS];
  
  Float_t xaxis[DRS_N_POINTS];
//   TGraph *gr[NCH_DRS];
  TGraphErrors *gr[NCH_DRS];
  Float_t yerr[DRS_N_POINTS];
  TF1 *erf;
  
  ULong64_t ts_drs;
  ULong64_t ts_tel;
  ULong64_t ts_anchor;
//   UInt_t ts_anchor;
  Int_t tel_trgoffset;        // to be removed
  ULong64_t tel_frameNumber;
  Int_t tel_triggerOffset;
  Int_t tel_triggerInfo;
  Bool_t tel_invalid;
  
  Long64_t delta_drs;
  Long64_t delta_tel;
  UInt_t delta_FEI4;
  float ratio;
  bool invalid;       // telescope invalid events // to be removed
  
  //telescope tracks
  Double_t x[1000];
  Double_t y[1000];
  Double_t errx[1000];
  Double_t erry[1000];
  Double_t slopeX[1000];
  Double_t slopeY[1000];
  Double_t slopeErrX[1000];
  Double_t slopeErrY[1000];
  Double_t chi2[1000];
  int n;
  
  //telescope plane hits
  Int_t hits_NHits[8];
  Double_t hits_Value[8][1000];
  Int_t hits_Tot[8][1000];       // for FEI4
  Double_t hits_Timing[8][1000];
  Int_t hits_PixX[8][1000];
  Int_t hits_PixY[8][1000];
  Int_t hits_HitInCluster[8][1000];
  
  ULong64_t ts_tel_s;   // For anchor module
//   Int_t PixX[1000];     // FE pixel hit
//   Int_t PixY[1000];
  UShort_t ClusterCharge[1000];    
  Float_t PosX[1000];
  Float_t PosY[1000];
  Double_t xa[1000];    // Cluster coordinates from PixX, PixY, ToT
  Double_t ya[1000];
  int n_a;
  UShort_t ClusterSize[1000];
  
  void Draw(int ch, const char* option="");
//   float GetCharge(int skip, float thr, int window); 
  float GetCharge(int ch, int bin_start, int bin_end, int mode=0, int pulse_polarity=1);
  float GetChargeCSA(int ch, int bin_start, int bin_end, int mode=0);
  int SearchTracks(window w);
  int SearchOneTrack(window w);
  bool TrackInRegion(int it, window w);
  int NTracksInRegion(window w);
  int SearchTrack(window w);
  
  bool VetoPresent(int ch, float thr, int mode=0); 
//   bool VetoPresent(int ch, float thr, int first=0, int last=NCH_DRS);
  bool ClusterizeAnchor();
};

struct analyzed_event
{
  float wf[NCH_DRS][NPTS_ANAWF];
  Float_t charge[NCH_DRS];
  
  Double_t xt;
  Double_t yt;
  Double_t chi2;
  Double_t slopeX;
  Double_t slopeY;
  int nt;
  
  Int_t ToT;    
  Double_t xa;    // Cluster coordinates from PixX, PixY, ToT
  Double_t ya;
  UShort_t clustersize;
  int na;
};
