waveform::waveform(){ for (int i=0; i<DRS_N_POINTS; i++) wf[i]=0; };
waveform::~waveform(){};
Float_t waveform::GetMin(){ return TMath::MinElement(DRS_N_POINTS, wf);};
Float_t waveform::GetMin(int first){ return TMath::MinElement(DRS_N_POINTS-2*first, wf+first);};
Float_t waveform::GetMin(int first, int last){ return TMath::MinElement(last-first, wf+first);};
Float_t waveform::GetMax(){ return TMath::MaxElement(DRS_N_POINTS, wf);};
Float_t waveform::GetMax(int first){ return TMath::MaxElement(DRS_N_POINTS-2*first, wf+first);};
Float_t waveform::GetMax(int first, int last){ return TMath::MaxElement(last-first, wf+first);};


Float_t waveform::GetMaxFiltered(int bin_start, int bin_end, float thr)
{
  // Filter out spikes
  // maximum searched in interval [bin_start ~50, bin_end ~ 300] maximal bin_end = 1024
  // thr ... amplitude above which to search for spikes at all - approx 0.003 (=3 mV)
  
  //if(bin_start > bin_end) bin_start=bin_end;
  
  int search_window=20;		// number of points around maximum to search for spike
  int loc_max = (int)(TMath::LocMax(bin_end-bin_start, wf+bin_start));
   
  Float_t max = wf[ TMath::Min( loc_max+bin_start, (int)(bin_end-1) )];
  if (max<thr) return max;
  
 // cout << i << " ....... " << loc_max << " ... " << max << " ... " << wf[loc_max + skip + search_window] << endl;
  if (wf[TMath::Min((int)(loc_max + bin_start + search_window), 1023)] < 0.7*max)
  {    
    // If spike detected enter here
    if (loc_max-search_window < 1) 
    {
      // take edges into account
      max = TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window);
    }
    else 
    {
      max = TMath::Max(
        TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+bin_start),
        TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window)
      );
    }
  }
  
  return max;
  
};
Float_t waveform::GetMaxFiltered(int skip, float thr, int window)
{
  // Filter out spikes
  // this method fails if there are two spikes in a single wf
  // thr ... amplitude above which to search for spikes at all - approx 0.003 (=3 mV)
  // skip ... number of initial steps skipped (some waveforms are messy in the beginning) - approx. 20

  /*
  int search_window=20;		// number of points around maximum to search for spike
  int loc_max = (int)(TMath::LocMax(DRS_N_POINTS-2*skip, wf+skip));
  Float_t max = wf[ TMath::Min( loc_max+skip, (int)(DRS_N_POINTS-1) )];
  if (max<thr) return max;
  
  //cout << i << " ....... " << loc_max << endl;
  if (wf[loc_max + skip + search_window] < 0.7*max)
  {
    max = TMath::Max(
      TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+skip),
      TMath::MaxElement( TMath::Max(DRS_N_POINTS-2*skip-loc_max-search_window, (long)(0)) , wf+skip+loc_max+search_window)
    );
  }
  
  return max;
*/
  
  int search_window=20;		// number of points around maximum to search for spike
  int loc_max = (int)(TMath::LocMax(window-2*skip, wf+skip));
    
  Float_t max = wf[ TMath::Min( loc_max+skip, (int)(window-1) )];
  if (max<thr) return max;
  
 // cout << i << " ....... " << loc_max << " ... " << max << " ... " << wf[loc_max + skip + search_window] << endl;
  if (wf[loc_max + skip + search_window] < 0.7*max)
  {
    if (loc_max-search_window < 1) max = TMath::MaxElement( TMath::Max(window-2*skip-loc_max-search_window, 0) , wf+skip+loc_max+search_window);
    else 
    {
      max = TMath::Max(
	TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+skip),
	TMath::MaxElement( TMath::Max(window-2*skip-loc_max-search_window, 0) , wf+skip+loc_max+search_window)
      );
    }
  }
  
  return max;
  
};
Float_t waveform::GetMaxMultiFiltered(int bin_start, int bin_end, float thr)
{
  // Filter out spikes
  // maximum searched in interval [bin_start ~50, bin_end ~ 300] maximal bin_end = 1024
  // thr ... amplitude above which to search for spikes at all - approx 0.003 (=3 mV)
  
  //if(bin_start > bin_end) bin_start=bin_end;
  
  int search_window=20;		// number of points around maximum to search for spike
  int loc_max = (int)(TMath::LocMax(bin_end-bin_start, wf+bin_start));
   
//   cout << bin_end << endl;
  Float_t max = wf[ TMath::Min( loc_max+bin_start, (int)(bin_end) )];
  if (max<thr) return max;
  
 // cout << i << " ....... " << loc_max << " ... " << max << " ... " << wf[loc_max + skip + search_window] << endl;
  if (wf[TMath::Min((int)(loc_max + bin_start + search_window), 1023)] < 0.7*max)
  {    
    // If spike detected enter here
    if (loc_max-search_window < 1) 
    {
      // take edges into account
//       cout << "case 1" << endl;
//       cout << bin_start+loc_max+search_window << " " << bin_end << endl;
//       max = TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window);
      max = GetMaxFiltered(
	TMath::Min(bin_start+loc_max+search_window, bin_end), 
	bin_end, 
	thr);
      
//       cout << "ended case 1" << endl;

    }
    else 
    {
//       max = TMath::Max(
// 	TMath::MaxElement( TMath::Max(loc_max-search_window, 0) , wf+bin_start),
// 	TMath::MaxElement( TMath::Max(bin_end-bin_start-2*bin_start-loc_max-search_window, 0) , wf+bin_start+loc_max+search_window)
//       );
            
//       cout << "case 2" << endl;

      max = TMath::Max(
	GetMaxFiltered(bin_start, bin_start+loc_max-search_window, thr),
	GetMaxFiltered(bin_start+loc_max+search_window, bin_end, thr)	
		   ); 
           
//       cout << "ended case 2" << endl;

    }
  }
  
  return max;
  
};
Float_t waveform::GetMinPos(){ return (float)(TMath::LocMin(DRS_N_POINTS, wf));};
Float_t waveform::GetMinPos(int skip){ return (float)(TMath::LocMin(DRS_N_POINTS-2*skip, wf+skip));};
Float_t waveform::GetMinPos(int first, int last){ return (float)(TMath::LocMin(last-first, wf+first));};
Float_t waveform::GetMaxPos(){ return (float)(TMath::LocMax(DRS_N_POINTS, wf));};
Float_t waveform::GetMaxPos(int skip){ return (float)(TMath::LocMax(DRS_N_POINTS-2*skip, wf+skip));};
Float_t waveform::GetMaxPos(int first, int last){ return (float)(TMath::LocMax(last-first, wf+first));};
Float_t waveform::Integral(int bin_start, int bin_end)
{
  Float_t integral=0;
  for (int i=bin_start; i<=bin_end; i++) integral+=wf[i];
  return integral;
}
Float_t waveform::MaxStep(float thr, int first, int last)
{
  // Get difference between max and min of the waveform
  // thr:   threshold for step detection (minimal delta V)
  // first: first sampled bin
  // last:  last sampled bin
  float dV=0;   // maximal delta U
  dV = GetMax(first, last) - GetMin(first, last);
  if (dV > thr) return dV;
  return 0.;
}
Float_t waveform::MaxStepN(int N, int sign, float thr, int dx, int first, int last)
{
  // Get difference between max and min of the waveform
  // Max and min are obtained by integrating over N consequtive points
  // N:     number of integrated points
  // sign:  step polarity (+1/-1)
  // thr:   threshold for step detection (minimal delta V)
  // dx:    increment step
  // first: first sampled bin
  // last:  last sampled bin
  float sum=0;   // maximal delta U
  float min=22222;
  float max=-22222;
  for (int i=TMath::Min(first, (int)(DRS_N_POINTS-N-10)); i<last-N; i+=dx)
  {
    sum=0;
    for (int j=0; j<N; j++) sum+=wf[i+j];
    sum/=N;
    if (i==first) {
      min=sum; max=sum;
    }
    if (sum>max) max=sum;
    else if (sum<min) min=sum;
  }
//   cout << first << " " << last << " " << max << " " << min << endl; 
  float dV=max-min;
//   cout << first << " " << dV << endl;
  if (dV > thr) return dV;
  return 0.;
}
Float_t waveform::SearchStep(int sign, float thr, int width, int step, int first, int last)
{
  // search waveform for a step between two straight lines
  // sign:  step polarity (+1/-1)
  // thr:   threshold for step detection (minimal delta V)
  // width: maximal width of the step (high-pass threshold) in bins
  // step:  step (in bins) at which the waveform is sampled
  // first: first sampled bin
  // last:  last sampled bin
  float dV=0;   // maximal delta U
  dV = GetMax(first, last) - GetMin(first, last);
  if (dV > thr) return dV;
  return 0.;
}
Float_t waveform::GetDelta(int first, int last, int N)
{
  // Get the difference in the voltage between bins first and last
  // N bins in the interval [first, first+N-1] and [last, last+N-1] are averaged
  float v1=0, v2=0;
  for (int j=0; j<N; j++) {
    v1+=wf[first+j];
    v2+=wf[last+j];
  }
  float delta=(v2-v1)/N;
//   cout << first << " " << last << " " << delta << " " << wf[last] - wf[first] << endl;;
  return delta;
}

TGraph* waveform::Draw()
{
  float t[DRS_N_POINTS];  
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) 
  {
//     wf[ipt] /= r->N;
    t[ipt] = ipt;
  }  
  TGraph* g = new TGraph (DRS_N_POINTS, t, wf);  
  g->GetXaxis()->SetTitle("t (ns)");
  g->GetYaxis()->SetTitle("U (mV)");
  g->SetTitle("Average waveform");
//   DrawTH1(g->GetHistogram(), "");  // Cosmetics for plotting
  return g;
}