#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"

#include <iostream>

struct window
{
  window(){};
  window(float xmin, float xmax, float ymin, float ymax) : x1(xmin), x2(xmax), y1(ymin), y2(ymax) {};
//   window(const window&){};
  bool Contains(float x, float y) {if(x>=x1 && x<=x2 && y>=y1 && y<=y2 ) return true; else return false; };
  void Set(float xmin, float xmax, float ymin, float ymax){x1=xmin, x2=xmax, y1=ymin, y2=ymax; return;};
  float x1, x2, y1, y2;  
};

class waveform
{
public:
  waveform();
  ~waveform();

  Float_t wf[DRS_N_POINTS];
  Float_t GetMin();
  Float_t GetMin(int first);
  Float_t GetMin(int first, int last);
  Float_t GetMax();
  Float_t GetMax(int first);
  Float_t GetMax(int first, int last);
  Float_t GetMaxFiltered(int bin_start, int bin_end, float thr);
  Float_t GetMaxFiltered(int skip, float thr, int window);
  Float_t GetMaxMultiFiltered(int bin_start, int bin_end, float thr);
  Float_t GetMaxPos();
  Float_t GetMaxPos(int skip);
  Float_t GetMaxPos(int first, int last);
  Float_t GetMinPos();
  Float_t GetMinPos(int skip);
  Float_t GetMinPos(int first, int last);
  Float_t Integral(int bin_start, int bin_end);
  Float_t MaxStep(float thr, int first=0, int last=DRS_N_POINTS-1);
  Float_t MaxStepN(int N, int sign, float thr, int dx, int first=0, int last=DRS_N_POINTS-1);
  Float_t SearchStep(int sign, float thr, int width, int step, int first=0, int last=DRS_N_POINTS-1);
  Float_t GetDelta(int bin_start, int bin_end, int N=10);
  TGraph *Draw();
};