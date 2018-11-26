#include "TH2S.h"

class DUTCorrelation
{
public:
  DUTCorrelation();
  DUTCorrelation(window ROI, int nbinsx, int nbinsy, int nch=NCH_DRS);
  ~DUTCorrelation();
  
  void ProcessEvent(event &Event, float *thr);
  TH2S *AddWeight(TH2S* h);
  float GetCorrelation(int ch_probe, int ch_ref);
  TH2S *_hHits[NCH_DRS][NCH_DRS];
  TH2S *_hWeights;
  
private:
  window _ROI;
  int _nch;
};