// #include "DUTCorrelation.cxx"
#include "class.hxx"

// waveform::waveform(){ for (int i=0; i<DRS_N_POINTS; i++) wf[i]=0; };
// waveform::~waveform(){};
// Float_t waveform::GetMin(){ return TMath::MinElement(DRS_N_POINTS, wf);};
// Float_t waveform::GetMin(int first){ return TMath::MinElement(DRS_N_POINTS-2*first, wf+first);};
// Float_t waveform::GetMin(int first, int last){ return TMath::MinElement(last-first, wf+first);};
// Float_t waveform::GetMax(){ return TMath::MaxElement(DRS_N_POINTS, wf);};
// Float_t waveform::GetMax(int first){ return TMath::MaxElement(DRS_N_POINTS-2*first, wf+first);};
// Float_t waveform::GetMax(int first, int last){ return TMath::MaxElement(last-first, wf+first);};
// 
// 
// Float_t waveform::GetMaxFiltered(int bin_start, int bin_end, float thr)
// {
//   // Filter out spikes
//   // maximum searched in interval [bin_start ~50, bin_end ~ 300] maximal bin_end = 1024
//   // thr ... amplitude above which to search for spikes at all - approx 0.003 (=3 mV)
//   
//   //if(bin_start > bin_end) bin_start=bin_end;
//   
//   int search_window=20;		// number of points around maximum to search for spike
//   int loc_max = (int)(TMath::LocMax(bin_end-bin_start, wf+bin_start));
//    
//   Float_t max = wf[ TMath::Min( loc_max+bin_start, (int)(bin_end-1) )];
//   if (max<thr) return max;
//   
//  // cout << i << " ....... " << loc_max << " ... " << max << " ... " << wf[loc_max + skip + search_window] << endl;
//   if (wf[TMath::Min((int)(loc_max + bin_start + search_window), 1023)] < 0.7*max)
//   {    
//     // If spike detected enter here
//     if (loc_max-search_window < 1) 
//     {
//       // take edges into account
//       max = TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window);
//     }
//     else 
//     {
//       max = TMath::Max(
//         TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+bin_start),
//         TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window)
//       );
//     }
//   }
//   
//   return max;
//   
// };
// Float_t waveform::GetMaxFiltered(int skip, float thr, int window)
// {
//   // Filter out spikes
//   // this method fails if there are two spikes in a single wf
//   // thr ... amplitude above which to search for spikes at all - approx 0.003 (=3 mV)
//   // skip ... number of initial steps skipped (some waveforms are messy in the beginning) - approx. 20
// 
//   /*
//   int search_window=20;		// number of points around maximum to search for spike
//   int loc_max = (int)(TMath::LocMax(DRS_N_POINTS-2*skip, wf+skip));
//   Float_t max = wf[ TMath::Min( loc_max+skip, (int)(DRS_N_POINTS-1) )];
//   if (max<thr) return max;
//   
//   //cout << i << " ....... " << loc_max << endl;
//   if (wf[loc_max + skip + search_window] < 0.7*max)
//   {
//     max = TMath::Max(
//       TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+skip),
//       TMath::MaxElement( TMath::Max(DRS_N_POINTS-2*skip-loc_max-search_window, (long)(0)) , wf+skip+loc_max+search_window)
//     );
//   }
//   
//   return max;
// */
//   
//   int search_window=20;		// number of points around maximum to search for spike
//   int loc_max = (int)(TMath::LocMax(window-2*skip, wf+skip));
//     
//   Float_t max = wf[ TMath::Min( loc_max+skip, (int)(window-1) )];
//   if (max<thr) return max;
//   
//  // cout << i << " ....... " << loc_max << " ... " << max << " ... " << wf[loc_max + skip + search_window] << endl;
//   if (wf[loc_max + skip + search_window] < 0.7*max)
//   {
//     if (loc_max-search_window < 1) max = TMath::MaxElement( TMath::Max(window-2*skip-loc_max-search_window, 0) , wf+skip+loc_max+search_window);
//     else 
//     {
//       max = TMath::Max(
// 	TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+skip),
// 	TMath::MaxElement( TMath::Max(window-2*skip-loc_max-search_window, 0) , wf+skip+loc_max+search_window)
//       );
//     }
//   }
//   
//   return max;
//   
// };
// Float_t waveform::GetMaxMultiFiltered(int bin_start, int bin_end, float thr)
// {
//   // Filter out spikes
//   // maximum searched in interval [bin_start ~50, bin_end ~ 300] maximal bin_end = 1024
//   // thr ... amplitude above which to search for spikes at all - approx 0.003 (=3 mV)
//   
//   //if(bin_start > bin_end) bin_start=bin_end;
//   
//   int search_window=20;		// number of points around maximum to search for spike
//   int loc_max = (int)(TMath::LocMax(bin_end-bin_start, wf+bin_start));
//    
// //   cout << bin_end << endl;
//   Float_t max = wf[ TMath::Min( loc_max+bin_start, (int)(bin_end) )];
//   if (max<thr) return max;
//   
//  // cout << i << " ....... " << loc_max << " ... " << max << " ... " << wf[loc_max + skip + search_window] << endl;
//   if (wf[TMath::Min((int)(loc_max + bin_start + search_window), 1023)] < 0.7*max)
//   {    
//     // If spike detected enter here
//     if (loc_max-search_window < 1) 
//     {
//       // take edges into account
// //       cout << "case 1" << endl;
// //       cout << bin_start+loc_max+search_window << " " << bin_end << endl;
// //       max = TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window);
//       max = GetMaxFiltered(
// 	TMath::Min(bin_start+loc_max+search_window, bin_end), 
// 	bin_end, 
// 	thr);
//       
// //       cout << "ended case 1" << endl;
// 
//     }
//     else 
//     {
// //       max = TMath::Max(
// // 	TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+bin_start),
// // 	TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window)
// //       );
//             
// //       cout << "case 2" << endl;
// 
//       max = TMath::Max(
// 	GetMaxFiltered(bin_start, bin_start+loc_max-search_window, thr),
// 	GetMaxFiltered(bin_start+loc_max+search_window, bin_end, thr)	
// 		   ); 
//            
// //       cout << "ended case 2" << endl;
// 
//     }
//   }
//   
//   return max;
//   
// };
// Float_t waveform::GetMinPos(){ return (float)(TMath::LocMin(DRS_N_POINTS, wf));};
// Float_t waveform::GetMinPos(int skip){ return (float)(TMath::LocMin(DRS_N_POINTS-2*skip, wf+skip));};
// Float_t waveform::GetMinPos(int first, int last){ return (float)(TMath::LocMin(first-last, wf+first));};
// Float_t waveform::GetMaxPos(){ return (float)(TMath::LocMax(DRS_N_POINTS, wf));};
// Float_t waveform::GetMaxPos(int skip){ return (float)(TMath::LocMax(DRS_N_POINTS-2*skip, wf+skip));};
// Float_t waveform::GetMaxPos(int first, int last){ return (float)(TMath::LocMax(first-last, wf+first));};
// Float_t waveform::Integral(int bin_start, int bin_end)
// {
//   Float_t integral=0;
//   for (int i=bin_start; i<=bin_end; i++) integral+=wf[i];
//   return integral;
// }
// Float_t waveform::MaxStep(float thr, int first, int last)
// {
//   // Get difference between max and min of the waveform
//   // thr:   threshold for step detection (minimal delta V)
//   // first: first sampled bin
//   // last:  last sampled bin
//   float dV=0;   // maximal delta U
//   dV = GetMax(first, last) - GetMin(first, last);
//   if (dV > thr) return dV;
//   return 0.;
// }
// Float_t waveform::MaxStepN(int N, int sign, float thr, int dx, int first, int last)
// {
//   // Get difference between max and min of the waveform
//   // Max and min are obtained by integrating over N consequtive points
//   // N:     number of integrated points
//   // sign:  step polarity (+1/-1)
//   // thr:   threshold for step detection (minimal delta V)
//   // dx:    increment step
//   // first: first sampled bin
//   // last:  last sampled bin
//   float sum=0;   // maximal delta U
//   float min=22222;
//   float max=-22222;
//   for (int i=TMath::Min(first, (int)(DRS_N_POINTS-N-10)); i<last-N; i+=dx)
//   {
//     sum=0;
//     for (int j=0; j<N; j++) sum+=wf[i+j];
//     sum/=N;
//     if (i==first) {
//       min=sum; max=sum;
//     }
//     if (sum>max) max=sum;
//     else if (sum<min) min=sum;
//   }
// //   cout << first << " " << last << " " << max << " " << min << endl; 
//   float dV=max-min;
// //   cout << first << " " << dV << endl;
//   if (dV > thr) return dV;
//   return 0.;
// }
// Float_t waveform::SearchStep(int sign, float thr, int width, int step, int first, int last)
// {
//   // search waveform for a step between two straight lines
//   // sign:  step polarity (+1/-1)
//   // thr:   threshold for step detection (minimal delta V)
//   // width: maximal width of the step (high-pass threshold) in bins
//   // step:  step (in bins) at which the waveform is sampled
//   // first: first sampled bin
//   // last:  last sampled bin
//   float dV=0;   // maximal delta U
//   dV = GetMax(first, last) - GetMin(first, last);
//   if (dV > thr) return dV;
//   return 0.;
// }
// Float_t waveform::GetDelta(int first, int last, int N)
// {
//   // Get the difference in the voltage between bins first and last
//   // N bins in the interval [first, first+N-1] and [last, last+N-1] are averaged
//   float v1=0, v2=0;
//   for (int j=0; j<N; j++) {
//     v1+=wf[first+j];
//     v2+=wf[last+j];
//   }
//   float delta=(v2-v1)/N;
//   return delta;
// }
// 
// TGraph* waveform::Draw()
// {
//   float t[DRS_N_POINTS];  
//   for(int ipt=0; ipt<DRS_N_POINTS; ipt++) 
//   {
// //     wf[ipt] /= r->N;
//     t[ipt] = ipt;
//   }  
//   TGraph* g = new TGraph (DRS_N_POINTS, t, wf);  
//   g->GetXaxis()->SetTitle("t (ns)");
//   g->GetYaxis()->SetTitle("U (mV)");
//   g->SetTitle("Average waveform");
//   DrawTH1(g->GetHistogram(), "");  // Cosmetics for plotting
//   return g;
// }

// event::event(){for(int i=0;i<NCH_DRS;i++) wf[i] = new waveform(); };
// event::~event()
// {
//   for(int i=0;i<2;i++) 
//   {
//     delete wf[i]; 
//     delete gr[i];
//   }
// };
// void event::Draw(int ch, char* option)
// {
//   Float_t x[DRS_N_POINTS];
//   for (int i=0; i<DRS_N_POINTS; i++) x[i] = 1.0*i;
//   
//   if(gr[ch]) delete gr[ch];
//   gr[ch] = new TGraph(DRS_N_POINTS, x, wf[ch].wf);
//   
//   gr[ch]->SetLineColor(1+ch);
//   gr[ch]->SetMarkerColor(1+ch);
//   gr[ch]->GetYaxis()->SetRangeUser(-0.05,0.2);
//   gr[ch]->SetTitle(" ;t (ns); U (V)");
//   
//   if (strlen(option)) gr[ch]->Draw(option);
//   else gr[ch]->Draw("APL");
//   return; //*/
// };
// 
// float event::GetCharge(int ch, int bin_start, int bin_end, int mode)
// {
//   switch(mode) {
//     default:
// //       return wf[ch].Integral(bin_start, bin_end) / (bin_end - bin_start +1);
//       _charge[ch] = wf[ch].Integral(bin_start, bin_end) / (bin_end - bin_start +1);
//       break;
//     case 1:
// //       return wf[ch].MaxStepN(10, 1, 0., 5, bin_start, bin_end);
// //       return wf[ch].GetDelta(bin_start, bin_end, 10);
//       _charge[ch] = wf[ch].GetDelta(bin_start, bin_end, 10);
//       break;
//   }
//   return _charge[ch];
//   //   return wf[ch].GetMax(bin_start, bin_end);
// 
//   //   return wf[ch].GetMaxMultiFiltered(bin_start, bin_end, thr);
//   //return wf[ch].GetMaxFiltered(skip, thr, window);
// };
// 
// int event::SearchTrack(window w)
// {
//   // Returns first track within the selected region
// //   int nt=0;
//   for(int j=0; j<n; j++)
//   {
//     if(x[j] >= w.x1 && x[j] <= w.x2)
//       if(y[j] >= w.y1 && y[j] <= w.y2)
//         return j;
//   }
//   return -1;
// }
// 
// int event::SearchOneTrack(window w)
// {
//   // Returns index of a track within the selected region. If more than one track passed, then return error
// //   int nt=0;
//   int i=-2;
//   for(int j=0; j<n; j++)
//   {
//     if(x[j] >= w.x1 && x[j] <= w.x2)
//       if(y[j] >= w.y1 && y[j] <= w.y2)
//       {
//         if(i==-2) i=j;    // first track
//         else i=-1;        // any additionaly tracks
//       }
//   }
//   return i;
// }
// 
// bool event::TrackInRegion(int it, window w)
// {
//   // Determines whether a track passes through the selected region w
// 
//   if(x[it] >= w.x1 && x[it] <= w.x2)
//     if(y[it] >= w.y1 && y[it] <= w.y2)
//       return true;
//   
//   return false;
// }
// 
// int event::NTracksInRegion(window w)
// {
//   // Returns number of tracks through the selected region
//   int count=0;
//   for (int it=0; it<n; it++)
//   {
//     if(x[it] >= w.x1 && x[it] <= w.x2)
//       if(y[it] >= w.y1 && y[it] <= w.y2)
//         count++;
//   }
//   return count;
// }
// 
// bool event::VetoPresent(int ch, float thr, int first, int last)
// {
//   // Seaches for a veto signal in one of waveforms 
//   // Veto signal has a box form so it can be discovered simply by finding the minimum and maximum of the wf
// //   float delta = wf[ch].GetMax(first) - wf[ch].GetMin(first);
//   
//   if ( (wf[ch].GetMax(first) - wf[ch].GetMin(first)) > thr) return true;
// //   if(delta > thr){
// //     cout << delta << " " << thr << " " << wf[ch].GetMax(first) << " " << wf[ch].GetMin(first) << " " << wf[ch].GetMax(first) - wf[ch].GetMin(first) << endl ;
// //     return true;
// //   }
//   return false;
// }

run::run(){Event = new event(); avwf = new waveform();};
run::run(char* fname_drs, char *brname_drs, char* lname_drs, int ch)
{
  Event = new event();
  avwf = new waveform();
  
  f_drs = new TFile(fname_drs);
  t_drs_wf = new TTree(fname_drs, "t_drs_wf");
  t_drs_wf = (TTree*) f_drs->Get(brname_drs);
//   t_drs_wf->SetBranchAddress(Form("%s%d",lname_drs,ch), Event.wf[0].wf);
  t_drs_wf->SetBranchAddress(Form("%s",lname_drs), Event.wf[ch].wf);
  N = t_drs_wf->GetEntries();
  
  t_drs_wf_ch2 = new TTree(fname_drs, "t_drs_wf");    // hack
  t_drs_wf_ch2 = (TTree*) f_drs->Get(brname_drs);
  t_drs_wf_ch2->SetBranchAddress("trigger2", Event.wf[1].wf);
    
};

run::run(char* fname_drs, char* fname_tel, char *brname_drs, char *brname_tel, vector<TString> lname_drs)
{
  Event = new event();
  avwf = new waveform();
  
  f_drs = new TFile(fname_drs);
  f_tel = new TFile(fname_tel);	
  
  nch = lname_drs.size();
  
  //t_drs_wf = new TTree(fname_drs, "t_drs_wf");
//   t_drs_wf = (TTree*) f_drs->Get(brname_drs);
//   t_drs_wf->SetBranchAddress(Form("%s%d",lname_drs,ch), Event.wf[0].wf);
//   t_drs_wf->SetBranchAddress(Form("%s",lname_drs[ch].Data()), Event.wf[ch].wf);
  
  for (int i=0; i<nch; i++)
  {
    v_t_drs_wf[i]= (TTree*) f_drs->Get(brname_drs);
    v_t_drs_wf[i]->SetBranchAddress(Form("%s",lname_drs[i].Data()), Event.wf[i].wf);
  }
  
  t_drs_ts = (TTree*) f_drs->Get("Event");
  t_drs_ts->SetBranchStatus("TimeStamp",1);
  t_drs_ts->SetBranchAddress("TimeStamp", &Event.ts_drs);
  
  t_tel_ts = (TTree*) f_tel->Get("Event");
  t_tel_ts->SetBranchStatus("TimeStamp",1);
  t_tel_ts->SetBranchAddress("TimeStamp", &Event.ts_tel);
  
  t_tel_invalid = (TTree*) f_tel->Get("Event");
  t_tel_invalid->SetBranchStatus("Invalid",1);
  t_tel_invalid->SetBranchAddress("Invalid", &Event.invalid);
  
  t_tel_tracks = (TTree*) f_tel->Get("Tracks");
  t_tel_tracks->SetBranchStatus("Origin*",1);
  t_tel_tracks->SetBranchAddress("OriginX", &Event.x);
  t_tel_tracks->SetBranchAddress("OriginY", &Event.y);
  t_tel_tracks->SetBranchAddress("NTracks", &Event.n);
  
  //t_tel_ts->AddFriend(t_tel_tracks);
  
  N = TMath::Min(t_drs_ts->GetEntries(), t_tel_ts->GetEntries());
  rel_offset = 0;
  
  ratio_mean = 1.25;		// expected delta_drs / delta_tel
  ratio_tolerance = 1.5;		// allowed deviation from expected ratio (factor)
  ratio_low = ratio_mean / ratio_tolerance;
  ratio_high = ratio_mean*ratio_tolerance;
  
  for (int i=0; i<NCH_DRS; i++) 
  {
    sigw_min[i]=0; sigw_max[i]=1;
    noisew_min[i]=0; noisew_max[i]=1;
  }
  //in_pair=false;
  
  hTracks = new TH2S(Form("tracks %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); tracks", 100,-10e3,10e3, 100, -6e3,6e3); 
//   hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); U_{CSA} (V)", 50,-6500,-6000, 50,-2200,-1800); 
//     hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); U_{CSA} (V)", 50,-6500,-6000, 50,-2200,-1750); 

};

run::run(char* fname_drs, char* fname_tel, char* fname_anchor, char *brname_drs, char *brname_tel, vector<TString> lname_drs)
{
  Event = new event();
  avwf = new waveform();
  
  f_drs = new TFile(fname_drs);
  f_tel = new TFile(fname_tel);	
  f_anchor = new TFile(fname_anchor);	

  nch = lname_drs.size();
  
  //t_drs_wf = new TTree(fname_drs, "t_drs_wf");
//   t_drs_wf = (TTree*) f_drs->Get(brname_drs);
//   t_drs_wf->SetBranchAddress(Form("%s%d",lname_drs,ch), Event.wf[0].wf);
//   t_drs_wf->SetBranchAddress(Form("%s",lname_drs[ch].Data()), Event.wf[ch].wf);
  
  for (int i=0; i<nch; i++)
  {
    v_t_drs_wf[i]= (TTree*) f_drs->Get(brname_drs);
    v_t_drs_wf[i]->SetBranchAddress(Form("%s",lname_drs[i].Data()), Event.wf[i].wf);
  }
  
  t_drs_ts = (TTree*) f_drs->Get("Event");
  t_drs_ts->SetBranchStatus("TimeStamp",1);
  t_drs_ts->SetBranchAddress("TimeStamp", &Event.ts_drs);
  
  t_tel_ts = (TTree*) f_tel->Get("Event");
  t_tel_ts->SetBranchStatus("TimeStamp",1);
  t_tel_ts->SetBranchAddress("TimeStamp", &Event.ts_tel);
  
  t_tel_trgoffset = (TTree*) f_tel->Get("Event");
  t_tel_trgoffset->SetBranchStatus("TriggerOffset",1);
  t_tel_trgoffset->SetBranchAddress("TriggerOffset", &Event.tel_trgoffset);
  
  t_tel_invalid = (TTree*) f_tel->Get("Event");
  t_tel_invalid->SetBranchStatus("Invalid",1);
  t_tel_invalid->SetBranchAddress("Invalid", &Event.invalid);
  
  t_tel_tracks = (TTree*) f_tel->Get("Tracks");
  t_tel_tracks->SetBranchStatus("Origin*",1);
  t_tel_tracks->SetBranchAddress("OriginX", &Event.x);
  t_tel_tracks->SetBranchAddress("OriginY", &Event.y);
  t_tel_tracks->SetBranchAddress("NTracks", &Event.n);
  
  t_anchor_ts = (TTree*) f_anchor->Get("Event");
  t_anchor_ts->SetBranchStatus("TimeStamp",1);
  t_anchor_ts->SetBranchAddress("TimeStamp", &Event.ts_anchor);
  
  //t_tel_ts->AddFriend(t_tel_tracks);
  
  t_anchor_hits = (TTree*) f_anchor->Get("Plane0/Hits");
  t_anchor_hits->SetBranchAddress("PixX", &Event.x_a);
  t_anchor_hits->SetBranchAddress("PixY", &Event.y_a);
  t_anchor_hits->SetBranchAddress("NHits", &Event.n_a);
  
  N = TMath::Min(t_drs_ts->GetEntries(), t_tel_ts->GetEntries());
  rel_offset = 0;
  
  ratio_mean = 1.25;		// expected delta_drs / delta_tel
  ratio_tolerance = 1.2;		// allowed deviation from expected ratio (factor)
  ratio_low = ratio_mean / ratio_tolerance;
  ratio_high = ratio_mean*ratio_tolerance;
  
  for (int i=0; i<NCH_DRS; i++) 
  {
    sigw_min[i]=0; sigw_max[i]=1;
    noisew_min[i]=0; noisew_max[i]=1;
  }
  //in_pair=false;
  
  hTracks = new TH2S(Form("tracks %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); tracks", 100,-10e3,10e3, 100, -6e3,6e3); 
//   hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); U_{CSA} (V)", 50,-6500,-6000, 50,-2200,-1800); 
//     hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); U_{CSA} (V)", 50,-6500,-6000, 50,-2200,-1750); 

  hResX = new TH1I("resX", "Residuals X ; residual X (mm); entries", 200, -20, 5);
  hResY = new TH1I("resY", "Residuals Y ; residual Y (mm); entries", 200, 0, 15);
  hRes = new TH2I("resXY", "Residuals XY ; residual X (mm); residual Y (mm) entries", 100, -10,-9, 100, 6.5,7.5);
};

run::~run()
{
  f_drs->Close();
  delete f_drs; 
  delete t_drs_wf; 
  delete t_drs_ts; 
  f_tel->Close();
  delete f_tel; 
  delete t_tel_ts; 
  f_anchor->Close();
  delete f_anchor; 
  delete t_anchor_hits; 
  delete Event;
  delete avwf;
};  
  
void run::GetEvent(int ie)
{
//   t_drs_wf->GetEntry(ie - rel_offset);
  for(int i=0; i<nch;i++) v_t_drs_wf[i]->GetEntry(ie - rel_offset);
  t_tel_tracks->GetEntry(ie);
  return;
};

void run::GetEvent_s(int ie_s)
{
  // Get data from files synchronized with the anchor module
//   t_anchor_ts->GetEntry(ie_s-rel_offset);
  t_anchor_hits->GetEntry(ie_s - rel_offset);
  return;
};

void run::GetEventDRS(int ie)
{

  t_drs_wf->GetEntry(ie);
  t_drs_wf_ch2->GetEntry(ie); // hack
  return;
};

void run::GetTimeStamp(int ie)
{
  //t_drs_ts->GetEntry(ie);
  GetDrsTimestamp(ie);
  //t_tel_ts->GetEntry(ie);
  GetTelTimestamp(ie);
  return;
};
  
bool run::ProcessInvalid(int ie, int verbose)
{
  t_tel_invalid->GetEntry(ie);
  if (Event.invalid) 
  {
    if (verbose) cout << ie << " invalid" << endl;
    return false;
  }
  
  for(int i=0; i<nch;i++) v_t_drs_wf[i]->GetEntry(ie);
  t_tel_tracks->GetEntry(ie);
  
  return true;
};

void run::ProcessTimestamp(int ie, int verbose)
{
  GetRatio(ie, ie+rel_offset);			// Calculate ts ratio
  /*
  if (!RatioInRange(Event.ratio)) 
    TryPairRatio(ie, ie+rel_offset);		// If ratios dont match first try ratio of pairs
  */
  if (!RatioInRange(Event.ratio)) 
    rel_offset += BestOffset(ie,3, verbose);		// If pair ratios dont match try offseting
  
  return;
};

Long64_t run::GetDrsTimestamp(int ie) 
{
  if (ie<0 || ie > t_drs_ts->GetEntries()) return 0;
  
  t_drs_ts->GetEntry(ie);
  Long64_t _ts = Event.ts_drs;
  Long64_t _tsold = 0;
    
  // Check for 1e4 bug (sometimes - rarely - the 5th last digit of the timestamp has to be incremented by 1
  if (ie>0)
  {
    t_drs_ts->GetEntry(ie-1);
    _tsold = Event.ts_drs;
    if(_ts < _tsold) _ts += 1e4;		// bug fix : drs timestamps should be strictly increasing, but seems like a bug, which sometimes does not increment 5th least significant digit in time. This line seems to fix this
    Event.ts_drs = _ts;
  }
  return _ts;
};

Long64_t run::GetDrsDelta(int ie) 
{
  Long64_t _delta =  - (GetDrsTimestamp(ie-1) - GetDrsTimestamp(ie));
  Event.delta_drs = _delta;
  return _delta;
};

Long64_t run::GetTelTimestamp(int ie)
{
  if (ie<0 || ie > t_tel_ts->GetEntries()) 
  {
    //cout << "ts 0" << endl;
    return 0;
  }
  t_tel_ts->GetEntry(ie);
//   t_tel_trgoffset->GetEntry(ie);
  return Event.ts_tel;
//   return Event.ts_tel - Event.tel_trgoffset;
};

Long64_t run::GetTelDelta(int ie) 
{
  Long64_t _tsold = GetTelTimestamp(ie-1);
  Long64_t _ts = GetTelTimestamp(ie);
  
  if (_ts < _tsold) _ts += (ULong64_t)(pow(2,32));	// correct counter overflow
  
  Long64_t _delta = _ts - _tsold;
  Event.delta_tel = _delta;
  return _delta;
};

float run::GetRatio(int i0, int i1)
{
  GetDrsDelta(i0);
  GetTelDelta(i1);

  if (Event.delta_tel == 0) Event.delta_tel = 1;
  
  float ratio = 1e4 * Event.delta_drs / Event.delta_tel;	// 1e4 to get the ratio in the region of ~1
  Event.ratio = ratio;
  return ratio;
};

bool run::RatioInRange(float ratio)
{
  if(ratio > ratio_low && ratio < ratio_high) return true;
  return false;
};

float run::TryPairRatio(int i0, int i1)
{
  // Sometimes the ratios of two consecutive events are not within limits, 
  // but if the order of the timestamps is flipped in one of the two sets the ratios are ok.
  // This function shall be called when ratio fails. 
  // It calculates the ratio (ts_drs1+ts_drs2) / (ts_tel1+ts_tel2), which might be within the range
  // Added also a check with the neighbour before
  
//   cout << "pair " << i0 << " " << i1 << endl;
  
  Long64_t _delta_drs1, _delta_drs0, _delta_drs2, _delta_tel1, _delta_tel0, _delta_tel2;
  _delta_drs2 = GetDrsDelta(i0+1);
  _delta_drs0 = GetDrsDelta(i0-1);
  _delta_drs1 = GetDrsDelta(i0);
  _delta_tel2 = GetTelDelta(i1+1);
  _delta_tel0 = GetTelDelta(i1-1);
  _delta_tel1 = GetTelDelta(i1);
      
  float ratio_after = 1e4 * (_delta_drs1+_delta_drs2) / TMath::Max(_delta_tel1+_delta_tel2, (Long64_t)(1));	// 1e4 to get the ratio in the region of ~1
  float ratio_before = 1e4 * (_delta_drs0+_delta_drs1) / TMath::Max(_delta_tel0+_delta_tel1, (Long64_t)(1));	// 1e4 to get the ratio in the region of ~1
  
  //if(!RatioInRange(ratio)) cout << "pair(NO) " << i0 << " " << i1 << endl;
  if(RatioInRange(ratio_after) )
  {
    Event.ratio = ratio_after;
  }
  else Event.ratio = ratio_before;
  
  
  //in_pair = true;	// the first member of a pair 
  
  return Event.ratio;
};

int run::BestOffset(int ie, int nsteps, int verbose)
{
  // This function evaluates the ratio delta_drs / delta_tel
  // If the ratio is out of expected range (typical interval [0.65,2.6]) it tries to shift one of the timestamp columns (either drs or tel) to achieve better agreement
  // The ratio is calculated for <w> consecutive entires and the number of times it is within the acceptable range is counted by <n>
  // The relative offsets <o1> and <o2> are taken from the interval [0,nsteps) and the best set of offsets is one which yields highest <n>
  
  int w=10;		// number of consecutive events for inspection
  int n_best = -1;	// number of timestamps in agreement
  int n;		// used in loops
    
  int o0=0, o1=0, o0_best, o1_best;		
  
  if (verbose) cout << ie << "\t" << Event.ts_drs << "\t" << Event.delta_drs*1e4 << "\t" << Event.delta_tel << "\t" << Event.delta_drs*1.0e4/TMath::Max(Event.delta_tel, (Long64_t)(1));
    
  // in the first loop try shifting both sets for a same offset, in the second loop also include relative offsets
  for (o0=0; o0<nsteps; o0++)
    for (o1=0; o1<nsteps; o1++)
    {
      if (o0 != o1) continue;
      n=0;
      for (int i=0; i<w; i++)
      {
	GetRatio(ie+o0+i, ie+o1+i+rel_offset);
	if(RatioInRange(Event.ratio)) n++;
      }
      if (n>n_best)
      {
	n_best = n;
	o0_best=o0;
	o1_best=o1;
      }
    }
  
  for (o0=0; o0<nsteps; o0++)
    for (o1=0; o1<nsteps; o1++)
    {
      if (o0 == o1) continue;
      n=0;
      for (int i=0; i<w; i++)
      {
	GetRatio(ie+o0+i, ie+o1+i+rel_offset);
	if(RatioInRange(Event.ratio)) n++;
      }
      if (n>n_best)
      {
	n_best = n;
	o0_best=o0;
	o1_best=o1;
      }
    }
    
  if (verbose)
  {
    cout << "  \t" << n_best << "\t" << o0_best << " " << o1_best << " " << rel_offset;
    if( o1_best - o0_best) 
    {
      cout << "\t incr. " <<  o1_best - o0_best;
      if(n_best <7) cout << " *";
    }
    cout << endl;
  }
  if(n_best <7) return 0;
  return o1_best - o0_best;
};

void run::FindSignalWindow(int nsamples, int delta_nbins, int ch)
{
  // Search for the position of signal peak in first nsamples waveforms
  // delta_nbins = number of bins around peak in avgwf to sample for maximum in individual events
  // In subsequent analysis the pulse height will be sampled only in the interval [peak-delta, peak+delta]
  // Only channel "ch" is sampled. If there are large delays between separate channels, then the funciton needs to be modified
  waveform w;
  int delta_noise=2*delta_nbins + 30;     // time before the signal sampling window in which to sample for noise
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf.wf[ipt] = 0;    // reset average waveform;

  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break;
    GetEvent(ie);
    
    for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf.wf[ipt] += Event.wf[ch].wf[ipt];
//     for(int ipt=0; ipt<DRS_N_POINTS; ipt++) w.wf[ipt] += Event.wf[ch].wf[ipt];
  }
//   for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf.wf[ipt] /= (ie+1);    // rescale in the end;
//   SetSignalWindow(w.GetMaxPos()-delta_nbins, w.GetMaxPos()+delta_nbins);
  SetSignalWindow(ch, avwf.GetMaxPos()-delta_nbins, avwf.GetMaxPos()+delta_nbins);
  SetNoiseWindow(ch, avwf.GetMaxPos()-delta_noise-delta_nbins, avwf.GetMaxPos()-delta_noise+delta_nbins);
  
  return;
}

void run::FindSignalWindow(int ch, int nsamples, int delta_low, int delta_high)
{
  // Search for the position of signal peak in first nsamples waveforms
  // delta_low/high = number of bins around peak in avgwf to sample for maximum in individual events
  // In subsequent analysis the pulse height will be sampled only in the interval [peak-delta_low, peak+delta_high]
  
//   waveform w;
  int delta_noise=(delta_low+delta_high) + 30;     // time before the signal sampling window in which to sample for noise
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf.wf[ipt] = 0;    // reset average waveform;

  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break;
    GetEvent(ie);
    
    for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf.wf[ipt] += Event.wf[ch].wf[ipt];
//     for(int ipt=0; ipt<DRS_N_POINTS; ipt++) w.wf[ipt] += Event.wf[ch].wf[ipt];
  }
//   for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf.wf[ipt] /= (ie+1);    // rescale in the end;
//   SetSignalWindow(w.GetMaxPos()-delta_nbins, w.GetMaxPos()+delta_nbins);
  SetSignalWindow(ch, avwf.GetMaxPos()-delta_low, avwf.GetMaxPos()+delta_high);
  SetNoiseWindow(ch, avwf.GetMaxPos()-delta_noise-delta_low, avwf.GetMaxPos()-delta_noise+delta_high);
  
  return;
}

void run::SetSignalWindow(int ch, float min, float max)
{
  // Set time window in which to search for signal peak in waveforms
  sigw_min[ch]=TMath::Max(min,0.0);
  sigw_max[ch]=TMath::Min(max, DRS_N_POINTS-1.0);
//   sigw_min[ch]=0.0;
//   sigw_max[ch]=DRS_N_POINTS-1.0;
  return;
}

void run::SetNoiseWindow(int ch, float min, float max)
{
  // Set time window in which to sample noise
  noisew_min[ch]=TMath::Max(min,0.0);
  noisew_max[ch]=TMath::Min(max, DRS_N_POINTS-1.0);
  return;
}

float run::GetNoise(int ch, int mode)
{  
//   return Event.GetCharge(ch,noisew_min[ch],noisew_max[ch],mode) / (-noisew_min[ch]+noisew_max[ch]+1);
  return Event.GetCharge(ch,noisew_min[ch],noisew_max[ch],mode);
//   float noise = Event.GetCharge(ch,noisew_min[ch],noisew_max[ch],mode);
//   switch(mode)
//   {
//     case 0:
//       return noise / (-noisew_min[ch]+noisew_max[ch]+1);
//       break;
//     case 1:
//       return noise;
//       break;
//   }
};

float run::FindThreshold(int ch, int nsamples, float sigma, int mode)
{
  // Finds discrimination threshold for channel #ch
  // The function samples noise level in first [nsamples] samples and fills it in a histogram
  // The obtained distribution is fitted with a gaussian
  // The threshold is set [sigma] standard deviations above the noise level
  
  TH1I hNoise("hnoise","hnoise",100, -0.04, 0.04);  // histogram for saving noise levels
  TF1 f("fitf","gaus", hNoise.GetXaxis()->GetXmin(), hNoise.GetXaxis()->GetXmax());
  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break; // if requested event exceeds the number of events in the file
    GetEvent(ie);
    float charge = GetNoise(ch, mode);    // Get noise level using standard function -- pulse integral over noisewindow
//     hNoise.Fill(charge / (noisew_max[ch]-noisew_min[ch]+1));  // fill histogram with [integral] / [integration range]
    hNoise.Fill(charge);  // fill histogram with [integral] / [integration range]
  }
  hNoise.Fit("fitf","RQ");
  float thr = f.GetParameter(1) + sigma * f.GetParameter(2);  // Get gaussian width and multiply by desired sigma
  cout << "Threshold ch. " << ch << ": " << thr*1000 << " mV" << endl;
  
//   hNoise->Draw();
//   gPad->Update();
//   gPad->Print("ASD.png");
  return thr;
};

float run::GetEfficiency(int ch, float thr, int mode)
{ 
  int j=Event.SearchTrack(ROI); // index of the first track going through the roi
//   int j=Event.SearchOneTrack(ROI); // index of the track going through the roi
  float charge=-1111;
  if (j>=0)
  {
    float x = Event.x[j];
    float y = Event.y[j];
    heff_tracks[ch]->Fill(x,y);
//     charge = Event.GetCharge(ch,250,320,0);
//     charge = Event.GetCharge(ch,sigw_min[ch],sigw_max[ch],mode);
//     charge = Event.GetCharge(ch,sigw_min[ch],sigw_max[ch],mode) / (-sigw_min[ch]+sigw_max[ch]+1);
    charge = Event.GetCharge(ch,sigw_min[ch],sigw_max[ch],mode);
    if (charge > thr)
      heff_hits[ch]->Fill(x,y);
  }
  return charge;
};

void run::SetROI(float xmin, float xmax, float ymin, float ymax, int nbins)
{
  r->ROI = window(xmin,xmax,ymin,ymax);
  for (int i=0; i<NCH_DRS; i++)
  {
    heff_tracks[i] = new TH2S(Form("heff_tracks_%d", i+1), Form("tracks through sample ch. %d ; x (#mum) ; y (#mum); nTracks",i+1),nbins,xmin,xmax,nbins,ymin,ymax);
    heff_hits[i] = new TH2S(Form("heff_hits_%d", i+1), Form("hits through sample ch. %d ; x (#mum) ; y (#mum); hits",i+1),nbins,xmin,xmax,nbins,ymin,ymax);
    heff_efficiency[i] = new TH2F(Form("heff_efficiency_%d", i+1), Form("Sample efficiency ch. %d ; x (#mum) ; y (#mum); efficiency",i+1),nbins,xmin,xmax,nbins,ymin,ymax);
    heff_hits[i]->SetStats(0);
    heff_tracks[i]->SetStats(0);
    heff_efficiency[i]->SetStats(0);
  }
  hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); tracks", 50,xmin,xmax, 50,ymin,ymax);
};

float run::GetResidualX(int it, int jt)
{
//   float res = Event.x[it] + (-300000)*Event.slopeX[it] - 250.0*Event.x_a[jt];
  float res = Event.x[it] - 250.0*Event.x_a[jt];
//   return  0.1*Event.x_a[jt];
  return 1e-3*res;
};
float run::GetResidualY(int it, int jt)
{
//   float res = Event.y[it] + (-300000)*Event.slopeY[it] + 50.0*Event.y_a[jt];    // + sign because ROI module is mounted upside down
  float res = Event.y[it] + 50.0*Event.y_a[jt];
  return 1e-3*res;
};
void run::FillResiduals()
{
  for (int it=0; it<Event.n; it++)
    for (int jt=0; jt<Event.n_a; jt++)
    {
//       if(Event.n > 1 || Event.n_a>1) continue;
      hResX->Fill( GetResidualX(it,jt) );
      hResY->Fill( GetResidualY(it,jt) );
      hRes->Fill( GetResidualX(it,jt), GetResidualY(it,jt) );
    }
  return;
};

void run::FillResiduals(window anchor_ROI)
{
  for (int jt=0; jt<Event.n_a; jt++)
  {
    if ( !anchor_ROI.Contains(Event.x_a[jt], Event.y_a[jt])) continue;
    for (int it=0; it<Event.n; it++)
    {
//       if(Event.n > 1 || Event.n_a>1) continue;
      hResX->Fill( GetResidualX(it,jt) );
      hResY->Fill( GetResidualY(it,jt) );
      hRes->Fill( GetResidualX(it,jt), GetResidualY(it,jt) );
    }
  }
  return;
};

bool run::CheckResiduals(int itr[2])
{
  // Check if x,y residuals (between telescope track and FEI4 hit) are within the allowed window <ROI>. 
  // Only events with exactly one telescope track in the ROI are considered, to minimize the errors due to non-reconstructed tracks.
  // The function returns true if exactly one track with right residuals is found and false otherwise. If return true, the indices of the correcposnding tracks are saved in itr[0] (tele) and itr[1] (anchor).
  
  if (Event.NTracksInRegion(ROI) != 1) return false;  // take only events with a single telescope track in ROI
  int npassed=0;
  for (int it=0; it<Event.n; it++)
  {
    if( ! Event.TrackInRegion(it, ROI)) continue;  // Check if the track <it> goes through ROI
    for (int jt=0; jt<Event.n_a; jt++)
    {
      if (! ResCut.Contains(GetResidualX(it, jt), GetResidualY(it, jt))) continue;  // Check if the track passes the residual test
      if (npassed) return false;  // take only event with exactly 1 hit in fei4
      
      itr[0] = it;
      itr[1] = jt;
      npassed++;
    }
  }
  if (npassed == 1) return true;
  return false;
};

/*
int run::BestOffset(int ie, int nsteps)
{
  // This function evaluates the ratio delta_drs / delta_tel
  // If the ratio is out of expected range (typical interval [0.65,2.6]) it tries to shift one of the timestamp columns (either drs or tel) to achieve better agreement
  // The ratio is calculated for <w> consecutive entires and the number of times it is within the acceptable range is counted by <n>
  // The relative offsets <o1> and <o2> are taken from the interval [0,nsteps) and the best set of offsets is one which yields highest <n>
  
  int w=10;		// number of consecutive events for inspection
  int n_best = -1;	// number of timestamps in agreement
  int n;		// used in loops
    
  int o0=0, o1=0, o0_best, o1_best;		
  
  cout << ie << "\t" << Event.ts_drs << "\t" << Event.delta_drs*1e4 << "\t" << Event.delta_tel << "\t" << Event.delta_drs*1.0e4/TMath::Max(Event.delta_tel, (Long64_t)(1));
    
  for (o0=0; o0<nsteps; o0++)
    for (o1=0; o1<nsteps; o1++)
    {
      n=0;
      for (int i=0; i<w; i++)
      {
	GetRatio(ie+o0+i, ie+o1+i+rel_offset);
	if(RatioInRange(Event.ratio)) n++;
      }
      if (n>n_best)
      {
	n_best = n;
	o0_best=o0;
	o1_best=o1;
      }
    }
    
  
  cout << "  \t" << n_best << "\t" << o0_best << " " << o1_best << " " << rel_offset;
  if( o1_best - o0_best) 
  {
    cout << "\t incr. " <<  o1_best - o0_best;
    if(n_best <8) cout << " *";
  }
  cout << endl;
  
  if(n_best <8) return 0;
  return o1_best - o0_best;
};*/