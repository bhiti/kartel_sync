// #include "event.hpp"

event::event()
{
  for(int i=0;i<NCH_DRS;i++) wf[i] = new waveform(); 
  erf = new TF1("erf", "[3] + [0]*TMath::Erf((x-[1])/[2])", 0, 1000);
  for (int i=0;i<DRS_N_POINTS; i++) xaxis[i]=1.0*i;
};
event::~event()
{
  for(int i=0;i<2;i++) 
  {
    delete wf[i]; 
    delete gr[i];
    delete erf;
  }
};
void event::Draw(int ch, const char* option)
{
//   Float_t x[DRS_N_POINTS];
//   for (int i=0; i<DRS_N_POINTS; i++) x[i] = 1.0*i;
  
  float labelSize = 0.05;
  if(gr[ch]) delete gr[ch];
//   gr[ch] = new TGraph(DRS_N_POINTS, xaxis, wf[ch]->wf);
  gr[ch] = new TGraphErrors(DRS_N_POINTS, xaxis, wf[ch]->wf,0,0);
  
  gr[ch]->SetLineColor(1+ch);
  gr[ch]->SetLineWidth(2);
  gr[ch]->SetMarkerColor(1+ch);
  gr[ch]->GetYaxis()->SetRangeUser(-0.05,0.2);
  gr[ch]->SetTitle(" ;t (ns); U (V)");
  
  gr[ch]->GetXaxis()->SetTitleSize(labelSize);
  gr[ch]->GetXaxis()->SetTitleOffset(1.);
  gr[ch]->GetXaxis()->SetLabelSize(labelSize);
  gr[ch]->GetYaxis()->SetTitleSize(labelSize);
  gr[ch]->GetYaxis()->SetLabelSize(labelSize);
  gr[ch]->GetYaxis()->SetTitleOffset(1.);
  
  if (strlen(option)) gr[ch]->Draw(option);
  else gr[ch]->Draw("APL");
  return; //*/
};

float event::GetCharge(int ch, int bin_start, int bin_end, int mode, int pulse_polarity)
{
  switch(mode) {
    case 1:
      // charge as peak to peak of the waveform segment between [bin_start, bin_end]
      _charge[ch] = pulse_polarity*wf[ch]->GetDelta(bin_start, bin_end, 10);
      break;
    case 2:
      erf->SetRange((double) bin_start, (double) bin_end);
      erf->SetParameters(pulse_polarity*0.005, 0.5*(bin_start+bin_end), 10, wf[ch]->wf[bin_start]+0.01);   // 0...scale, 1...center, 2...width, 3...offset
      if(gr[ch]) {
        cout << " ----------------------------- " << endl;
        delete gr[ch];
      }
      gr[ch] = new TGraphErrors(bin_end - bin_start +1, xaxis + bin_start, wf[ch]->wf + bin_start, 0, yerr);
//       erf->Draw();
      gr[ch]->Fit("erf", "R");
      gr[ch]->Draw("APL");
      gr[ch]->GetYaxis()->SetRangeUser(0,0.1);
      erf->Draw("SAME");
      gPad->Update();
      
      char c;
      cin >> c;
      
      _charge[ch] = pulse_polarity*wf[ch]->GetDelta(bin_start, bin_end, 10);
      break;
    default:
      // charge as waveform integral between [bin_start, bin_end]
      _charge[ch] = pulse_polarity*wf[ch]->Integral(bin_start, bin_end) / (bin_end - bin_start +1);
      break;
    
  }
//   return (wf[ch]->wf[bin_end] - wf[ch]->wf[bin_start]);
  return _charge[ch];
  //   return wf[ch]->GetMax(bin_start, bin_end);

  //   return wf[ch]->GetMaxMultiFiltered(bin_start, bin_end, thr);
  //return wf[ch]->GetMaxFiltered(skip, thr, window);
};

int event::SearchTrack(window w)
{
  // Returns first track within the selected region
//   int nt=0;
  for(int j=0; j<n; j++)
  {
    if(x[j] >= w.x1 && x[j] <= w.x2)
      if(y[j] >= w.y1 && y[j] <= w.y2)
        return j;
  }
  return -1;
}

int event::SearchOneTrack(window w)
{
  // Returns index of a track within the selected region. If more than one track passed, then return error
//   int nt=0;
  int i=-2;
  for(int j=0; j<n; j++)
  {
    if(x[j] >= w.x1 && x[j] <= w.x2)
      if(y[j] >= w.y1 && y[j] <= w.y2)
      {
        if(i==-2) i=j;    // first track
        else i=-1;        // any additionaly tracks
      }
  }
  return i;
}

bool event::TrackInRegion(int it, window w)
{
  // Determines whether a track passes through the selected region w

  if(x[it] >= w.x1 && x[it] <= w.x2)
    if(y[it] >= w.y1 && y[it] <= w.y2)
      return true;
  
  return false;
}

int event::NTracksInRegion(window w)
{
  // Returns number of tracks through the selected region
  int count=0;
  for (int it=0; it<n; it++)
  {
    if(x[it] >= w.x1 && x[it] <= w.x2)
      if(y[it] >= w.y1 && y[it] <= w.y2)
        count++;
  }
  return count;
}

bool event::VetoPresent(int ch, float thr, int mode)
{
  // Seaches for a veto signal in one of waveforms 
  // Veto signal has a box form so it can be discovered simply by finding the minimum and maximum of the wf
  // mode : 
//   0 ... get minimum
//   1 ... get difference between min and max
  float delta;
  switch (mode){
    case 1:
      delta = (wf[ch]->GetMax() - wf[ch]->GetMin());
      if ( delta > thr) return true;
      break;
    case 2:
      delta = wf[ch]->GetMax();
      if ( delta > thr) return true;
      break;
    default:
      delta = wf[ch]->GetMin();
      if ( delta < thr) return true;
      break;
  }  
   
//   cout << delta << endl;
  return false;
}

/*
bool event::VetoPresent(int ch, float thr, int first, int last)
{
  // Seaches for a veto signal in one of waveforms 
  // Veto signal has a box form so it can be discovered simply by finding the minimum and maximum of the wf
//   float delta = wf[ch]->GetMax(first) - wf[ch]->GetMin(first);
  
  float delta = (wf[ch]->GetMax(first) - wf[ch]->GetMin(first));
  cout << delta << endl;
  if ( delta > thr) return true;
//   if(delta > thr){
//     cout << delta << " " << thr << " " << wf[ch]->GetMax(first) << " " << wf[ch]->GetMin(first) << " " << wf[ch]->GetMax(first) - wf[ch]->GetMin(first) << endl ;
//     return true;
//   }
  return false;
}
*/