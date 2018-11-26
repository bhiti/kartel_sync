// #include "src/histograms.hpp"

Histograms::Histograms(){};
Histograms::~Histograms()
{
  delete hTracks;
  delete hResX;
  delete hResY;
  delete hRes;

};
void Histograms::Init(int nch)
{
  hTracks = new TH2S("tracks", "Selected tracks ; x (#mum); y (#mum); tracks", 100,-10e3,10e3, 100, -6e3,6e3); 
  
  hResX = new TH1I("resX", "Residuals X ; residual X (mm); entries", 200, -20, 5);
  hResY = new TH1I("resY", "Residuals Y ; residual Y (mm); entries", 200, 0, 15);
  hRes = new TH2I("resXY", "Residuals XY ; residual X (mm); residual Y (mm) entries", 100, -10,-9, 100, 6.5,7.5);
    
  for (int i=0; i<nch; i++){
    hpulse_height[i] = new TH1S(Form("hpulse_height_%d", i+1), Form("Pulse height distribution ch. %d ; V (mV) ; entries",i+1),100,-50, 300);
    hpulse_height[i]->SetStats(0);
//     hC[i]=new TH1I(Form("hC%d",i+1), Form(" ; U_{out} (V) ; N_{entries}", i+1), 100,-0.01,0.09);
//     hN[i]=new TH1I(Form("hN%d",i+1), Form(" ; U_{out} (V) ; N_{entries}", i+1), 100,-0.01,0.01);
    hC[i]=new TH1I(Form("hC%d",i+1), Form(" ; U_{out} (V) ; N_{entries}"), 100,-0.05,0.3);
    hN[i]=new TH1I(Form("hN%d",i+1), Form(" ; U_{out} (V) ; N_{entries}"), 100,-0.03,0.03);
  }
  return;
};
  
void Histograms::SetHistosROI(window ROI, int nch, int nbinsx, int nbinsy)
{
  if (nbinsy<0) nbinsy=nbinsx;
  for (int i=0; i<nch; i++)
  {
    heff_tracks[i] = new TH2S(Form("heff_tracks_%d", i+1), Form("tracks through sample ch. %d ; x (#mum) ; y (#mum); nTracks",i+1),nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
    heff_hits[i] = new TH2S(Form("heff_hits_%d", i+1), Form("hits through sample ch. %d ; x (#mum) ; y (#mum); hits",i+1),nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
    heff_efficiency[i] = new TH2F(Form("heff_efficiency_%d", i+1), Form("Sample efficiency ch. %d ; x (#mum) ; y (#mum); efficiency",i+1),nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
    heff_charge[i] = new TH2F(Form("heff_charge_%d", i+1), Form("Pulse height ch. %d ; x (#mum) ; y (#mum); voltage (mV)",i+1),nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
    
    heff_hits[i]->SetStats(0);
    heff_tracks[i]->SetStats(0);
    heff_efficiency[i]->SetStats(0);
    heff_charge[i]->SetStats(0);
  }
  hTracks_ROI = new TH2S(Form("tracks ROI %s","000"), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); tracks", nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
  hClusterSize = new TH2S("cluster_size", "Cluster size in FEI4 ; x (#mum); y (#mum); average cluster size", nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2); 
  
  hclustered_tracks = new TH2S("hclustered_tracks", "hclustered_tracks ; x (#mum) ; y (#mum); nTracks", nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
  hclustered_charge = new TH2F("hclustered_charge", "hclustered_charge ; x (#mum) ; y (#mum); nTracks", nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
  hclustered_hits = new TH2S("hclustered_hits", "hclustered_hits ; x (#mum) ; y (#mum); hits", nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
  hclustered_efficiency = new TH2F("hclustered_efficiency", "hclustered_efficiency ; x (#mum) ; y (#mum); efficiency", nbinsx, ROI.x1, ROI.x2, nbinsy, ROI.y1, ROI.y2);
    
};