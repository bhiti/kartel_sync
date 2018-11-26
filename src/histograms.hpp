#include "TH1I.h"
#include "TH1S.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH2S.h"

class Histograms
{
public:
  Histograms();
  ~Histograms();
  TH1I* hC[NCH_DRS];  // charge spectrum
  TH1I* hN[NCH_DRS];  // noise spectrum
  
  TH2S *heff_tracks[NCH_DRS];
  TH2S *heff_hits[NCH_DRS];
  TH2F *heff_efficiency[NCH_DRS];
  TH2F *heff_charge[NCH_DRS];
  TH1S* hpulse_height[NCH_DRS];  // Pulse height spectrum

  TH2S *hclustered_tracks;
  TH2F *hclustered_charge;
  TH2S *hclustered_hits;
  TH2F *hclustered_efficiency;
  
  TH2S *hTracks;
  TH2S *hTracks_ROI;
  TH2S *hClusterSize;
  
  TH1I *hResX;
  TH1I *hResY;
  TH2I *hRes;
  
  void Init(int nch);
  void SetHistosROI(window ROI, int nch, int nbinsx, int nbinsy=-1);
};
