// #include "DUTCorrelation.cxx"
// #include "run.hpp"

run::run(){
  Event = new event(); avwf = new waveform();
};

run::run(const char* fname_drs, const char* fname_tel, const char *fname_anchor, const char *fname_results)
{
  Event = new event();
  
  IOHand = new IOHandler();
  IOHand->Init(Event, fname_drs, fname_tel, fname_anchor, fname_results);
  
//  N = IOHand->t_anchor_ts->GetEntries();
  N = TMath::Min(IOHand->t_anchor_ts->GetEntries(), IOHand->t_tel_ts->GetEntries());
  
  rel_offset = 0;
  abs_offset = 0;
  ratio_mean = 20;		// expected delta_drs / delta_tel
  ratio_tolerance = 1.2;	// Bojan 1.2	// allowed deviation from expected ratio (factor)
  ratio_low = ratio_mean / ratio_tolerance;
  ratio_high = ratio_mean*ratio_tolerance;
  
  hRatio = new TH1I("hRatio", "hRatio", 500, 0, 0.001);
  
//   hDeltaTimeStamp = new TH1I("hDelta", "Delta Time Stamp Telescope", 100, 0, 10000);
//   hDeltaVsInvalid = new TH2I("hDelta2D", "Delta Time Stamp Telescope", 1000, 0, 500000, 2, -0.5, 1.5);
  
};

/*
run::run(char* fname_drs, char *brname_drs, char* lname_drs, int ch)
{
  Event = new event();
  avwf = new waveform();
  
  f_drs = new TFile(fname_drs);
  t_drs_wf = new TTree(fname_drs, "t_drs_wf");
  t_drs_wf = (TTree*) f_drs->Get(brname_drs);
//   t_drs_wf->SetBranchAddress(Form("%s%d",lname_drs,ch), Event->wf[0].wf);
  t_drs_wf->SetBranchAddress(Form("%s",lname_drs), Event->wf[ch]->wf);
  N = t_drs_wf->GetEntries();
  
  t_drs_wf_ch2 = new TTree(fname_drs, "t_drs_wf");    // hack
  t_drs_wf_ch2 = (TTree*) f_drs->Get(brname_drs);
  t_drs_wf_ch2->SetBranchAddress("trigger2", Event->wf[1].wf);
    
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
//   t_drs_wf->SetBranchAddress(Form("%s%d",lname_drs,ch), Event->wf[0].wf);
//   t_drs_wf->SetBranchAddress(Form("%s",lname_drs[ch].Data()), Event->wf[ch]->wf);
  
  for (int i=0; i<nch; i++)
  {
    v_t_drs_wf[i]= (TTree*) f_drs->Get(brname_drs);
    v_t_drs_wf[i]->SetBranchAddress(Form("%s",lname_drs[i].Data()), Event->wf[i].wf);
  }
  
  t_drs_ts = (TTree*) f_drs->Get("Event");
  t_drs_ts->SetBranchStatus("TimeStamp",1);
  t_drs_ts->SetBranchAddress("TimeStamp", &Event->ts_drs);
  
  t_tel_ts = (TTree*) f_tel->Get("Event");
  t_tel_ts->SetBranchStatus("TimeStamp",1);
  t_tel_ts->SetBranchAddress("TimeStamp", &Event->ts_tel);
  
  t_tel_invalid = (TTree*) f_tel->Get("Event");
  t_tel_invalid->SetBranchStatus("Invalid",1);
  t_tel_invalid->SetBranchAddress("Invalid", &Event->invalid);
  
  t_tel_tracks = (TTree*) f_tel->Get("Tracks");
  t_tel_tracks->SetBranchStatus("Origin*",1);
  t_tel_tracks->SetBranchAddress("OriginX", &Event->x);
  t_tel_tracks->SetBranchAddress("OriginY", &Event->y);
  t_tel_tracks->SetBranchAddress("NTracks", &Event->n);
  
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
  
//   hTracks = new TH2S(Form("tracks %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); tracks", 100,-10e3,10e3, 100, -6e3,6e3); 
// //   hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); U_{CSA} (V)", 50,-6500,-6000, 50,-2200,-1800); 
// //     hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); U_{CSA} (V)", 50,-6500,-6000, 50,-2200,-1750); 
//   for (int i=0; i<NCH_DRS; i++){
//     hpulse_height[i] = new TH1S(Form("hpulse_height_%d", i+1), Form("Pulse height distribution ch. %d ; V (mV) ; entries",i+1),100,-50, 300);
//     hpulse_height[i]->SetStats(0);
//   }
};

run::run(char* fname_drs, char* fname_tel, char* fname_anchor, char* fname_results, char *brname_drs, char *brname_tel, vector<TString> lname_drs) : rel_offset(0)
{
  Event = new event();
  avwf = new waveform();
  
  f_drs = new TFile(fname_drs);
  f_tel = new TFile(fname_tel);	
  f_anchor = new TFile(fname_anchor);	
  f_res= new TFile(fname_results, "RECREATE");
  
  nch = lname_drs.size();
  
  //t_drs_wf = new TTree(fname_drs, "t_drs_wf");
//   t_drs_wf = (TTree*) f_drs->Get(brname_drs);
//   t_drs_wf->SetBranchAddress(Form("%s%d",lname_drs,ch), Event->wf[0].wf);
//   t_drs_wf->SetBranchAddress(Form("%s",lname_drs[ch].Data()), Event->wf[ch]->wf);
  
  for (int i=0; i<nch; i++)
  {
    v_t_drs_wf[i]= (TTree*) f_drs->Get(brname_drs);
    v_t_drs_wf[i]->SetBranchAddress(Form("%s",lname_drs[i].Data()), Event->wf[i].wf);
  }
  
  t_drs_ts = (TTree*) f_drs->Get("Event");
  t_drs_ts->SetBranchStatus("TimeStamp",1);
  t_drs_ts->SetBranchAddress("TimeStamp", &Event->ts_drs);
  
  t_tel_ts = (TTree*) f_tel->Get("Event");
  t_tel_ts->SetBranchStatus("TimeStamp",1);
  t_tel_ts->SetBranchAddress("TimeStamp", &Event->ts_tel);
  
//   t_tel_trgoffset = (TTree*) f_tel->Get("Event");
//   t_tel_trgoffset->SetBranchStatus("TriggerOffset",1);
//   t_tel_trgoffset->SetBranchAddress("TriggerOffset", &Event->tel_trgoffset);
//   
//   t_tel_invalid = (TTree*) f_tel->Get("Event");
//   t_tel_invalid->SetBranchStatus("Invalid",1);
//   t_tel_invalid->SetBranchAddress("Invalid", &Event->invalid);
//   
  t_tel_tracks = (TTree*) f_tel->Get("Tracks");
  t_tel_tracks->SetBranchStatus("Origin*",1);
  t_tel_tracks->SetBranchAddress("OriginX", &Event->x);
  t_tel_tracks->SetBranchAddress("OriginY", &Event->y);
  t_tel_tracks->SetBranchAddress("NTracks", &Event->n);
  
//   t_anchor_ts = (TTree*) f_anchor->Get("Event");
//   t_anchor_ts->SetBranchStatus("TimeStamp",1);
//   t_anchor_ts->SetBranchAddress("TimeStamp", &Event->ts_anchor);
  
//   t_anchor = (TTree*) f_anchor->Get("Event");
// //   t_anchor->SetBranchStatus("TimeStamp",1);
//   t_anchor->SetBranchAddress("NClusters", &Event->NClusters_a);
//   t_anchor->SetBranchAddress("ID", &Event->ClusterID_a);
//   t_anchor->SetBranchAddress("PosX", &Event->);
  
  //t_tel_ts->AddFriend(t_tel_tracks);
  
  t_anchor_hits = (TTree*) f_anchor->Get("Plane0/Hits");
//   t_anchor_hits->SetBranchAddress("PixX", &Event->PixX);
//   t_anchor_hits->SetBranchAddress("PixY", &Event->PixY);
  t_anchor_hits->SetBranchAddress("PosX", &Event->PosX);
  t_anchor_hits->SetBranchAddress("PosY", &Event->PosY);
  t_anchor_hits->SetBranchAddress("NClusters", &Event->n_a);
//   t_anchor_hits->SetBranchAddress("Value", &Event->ToT);
  
  N = TMath::Min(t_drs_ts->GetEntries(), t_tel_ts->GetEntries());
  rel_offset = 0;
  
  ratio_mean = 1.25;		// expected delta_drs / delta_tel
  ratio_tolerance = 1.2;		// allowed deviation from expected ratio (factor)
  ratio_low = ratio_mean / ratio_tolerance;
  ratio_high = ratio_mean*ratio_tolerance;
  
  for (int i=0; i<nch; i++) 
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
  
  for (int i=0; i<nch; i++){
    hpulse_height[i] = new TH1S(Form("hpulse_height_%d", i+1), Form("Pulse height distribution ch. %d ; V (mV) ; entries",i+1),100,-50, 300);
    hpulse_height[i]->SetStats(0);
  }
};*/

/*
run::run(char* fname_drs, char* fname_tel, char* fname_anchor, char* fname_results, char* fname_ana, char *brname_drs, char *brname_tel, vector<TString> lname_drs) : rel_offset(0)
{
  nch = lname_drs.size();
  
  Event = new event();
  avwf = new waveform();
  
  IOHand = new IOHandler();
  IOHand->Init(&Event, fname_drs, fname_tel, fname_anchor, fname_results, fname_ana, brname_drs, brname_tel, lname_drs);
  
  histos = new Histograms();
  histos->Init(nch);
  
  N = TMath::Min(IOHand->t_drs_ts->GetEntries(), IOHand->t_tel_ts->GetEntries());
  rel_offset = 0;
  
  ratio_mean = 1.25;		// expected delta_drs / delta_tel
  ratio_tolerance = 1.2;		// allowed deviation from expected ratio (factor)
  ratio_low = ratio_mean / ratio_tolerance;
  ratio_high = ratio_mean*ratio_tolerance;
  
  for (int i=0; i<nch; i++) 
  {
    sigw_min[i]=0; sigw_max[i]=1;
    noisew_min[i]=0; noisew_max[i]=1;
  }
};*/

run::~run()
{
  delete IOHand;
  delete avwf;
  delete Event;
  delete histos;
};  
  
void run::GetEvent(int ie)
{
//   for(int i=0; i<nch;i++) v_t_drs_wf[i]->GetEntry(ie - rel_offset);
//   t_tel_tracks->GetEntry(ie);
  IOHand->GetEvent(ie, rel_offset);
  
  // Apply slope correction to telescope track
//   for(int i=0; i<Event->n; i++){
//     Event->x[i] = Event->x[i] + zDUT * Event->slopeX[i];
//     Event->y[i] = Event->y[i] + zDUT * Event->slopeY[i];
//   }
  return;
};

void run::GetEvent_s(int ie_s)
{
  // Get data from files synchronized with the anchor module
//   t_anchor_ts->GetEntry(ie_s-rel_offset);
//   t_anchor_hits->GetEntry(ie_s - rel_offset);
  IOHand->GetEvent_s(ie_s, rel_offset);
  return;
};

void run::GetTimeStamp(int ie)
{
//   //t_drs_ts->GetEntry(ie);
//   GetDrsTimestamp(ie);
//   //t_tel_ts->GetEntry(ie);
//   GetTelTimestamp(ie);
  IOHand->GetTimeStamp(ie);
  return;
};
  
int run::ProcessTimestamp(int ie, int mode, int verbose)
{
//   cout << "*********************" << endl;
//   cout << ie << endl;
  float ratio = GetRatio(ie+rel_offset, ie, mode);			// Calculate ts ratio  
  hRatio->Fill(ratio);
  
//   cout << ratio << endl;
//   char c;
//   cin >> c;
  
//   hDeltaTimeStamp->Fill(Event->delta_FEI4);
//   hDeltaVsInvalid->Fill(Event->delta_FEI4, Event->tel_invalid ? 1 : 0);
//   if (Event->delta_FEI4 == 0) {
//     cout << ie << " " << Event->delta_FEI4 << " " << Event->delta_tel << endl;
//     bigOffset++;
//   }
//   
  if (!RatioInRange(Event->ratio)){
    offsets = BestOffset(ie,3, mode, verbose);  
    return 1;
  }
  
  return 0;
};

int run::ProcessOffsets(int ie)
{
  // local sync offsets
  int orel = offsets.first;
  int oabs = offsets.second;
    
//   cout << orel << " " << oabs << endl;
  //modify global offsets
  rel_offset += orel;   
  abs_offset += oabs;
  int iskipped=0;   // how many events to skip locally in the main loop
  //
  // first skip the number of events in both files
  //
  for (int i=0; i<oabs; i++){
//     FillSyncedFileDummy(ie); //Bojan
    iskipped++;   
  }
  
  //
  // then skip the number of events in a single file
  //
  for (int i=0; i<orel; i++){
//     FillSyncedFileDummy(ie); //Bojan
    iskipped++;
  }
  // if orel<0 do nothing
  
//   if (iskipped<0) iskipped=0;   // needs to be here in case only skipping in DUT file
  return iskipped;
}

Long64_t run::GetDrsTimestamp(int ie) 
{
  return IOHand->GetDrsTimestamp(ie);
};

Long64_t run::GetDrsDelta(int ie) 
{
  return IOHand->GetDrsDelta(ie);
};

Long64_t run::GetTelTimestamp(int ie)
{
//   if (ie<0 || ie > t_tel_ts->GetEntries()) 
//   {
//     //cout << "ts 0" << endl;
//     return 0;
//   }
//   t_tel_ts->GetEntry(ie);
// //   t_tel_trgoffset->GetEntry(ie);
//   return Event->ts_tel;
// //   return Event->ts_tel - Event->tel_trgoffset;

  return IOHand->GetTelTimestamp(ie);
};

Long64_t run::GetTelDelta(int ie) 
{
//   Long64_t _tsold = GetTelTimestamp(ie-1);
//   Long64_t _ts = GetTelTimestamp(ie);
//   
//   if (_ts < _tsold) _ts += (ULong64_t)(pow(2,32));	// correct counter overflow
//   
//   Long64_t _delta = _ts - _tsold;
//   Event->delta_tel = _delta;
//   return _delta;

  return GetTelDelta(ie);
};

float run::SetRatio(float mean, float tolerance)
{
  ratio_mean = mean;
  ratio_tolerance = tolerance;
  
  ratio_low = ratio_mean / ratio_tolerance;
  ratio_high = ratio_mean*ratio_tolerance;
  
  return ratio_mean;
}

float run::GetRatio(int i0, int i1, int mode)
{
  float d1, d2;
  switch(mode) {
    case 0:
      d1 = IOHand->GetDrsDelta(i0);
      d2 = IOHand->GetTelDelta(i1);
      break;
    case 1:
      d1 = IOHand->GetFEI4Delta(i0);
      d2 = IOHand->GetTelDelta(i1);
      break;
    case 2:
      d1 = IOHand->GetDrsDelta(i0);
      d2 = IOHand->GetFEI4Delta(i1);
      break;
    default:
      d1 = IOHand->GetDrsDelta(i0);
      d2 = IOHand->GetTelDelta(i1);
      break;
  }

//   cout << d1 << " " << d2 << " " << 1.0*d1/d2 << endl;
  if (d2 == 0) d2 = 1;
  
  float ratio = 1.0 * d1 / d2;	// 1e4 to get the ratio in the region of ~1
  Event->ratio = ratio;
//   cout << d1 << " " << d2 << " " << ratio << endl;
  return ratio;
  
};

bool run::RatioInRange(float ratio)
{
//   cout << "..." << ratio << " " << ratio_low << " " << ratio_high << endl;
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
    Event->ratio = ratio_after;
  }
  else Event->ratio = ratio_before;
  
  
  //in_pair = true;	// the first member of a pair 
  
  return Event->ratio;
};

std::pair<int, int> run::BestOffset(int ie, int nsteps, int mode, int verbose)
{
  // This function evaluates the ratio delta_drs / delta_tel
  // If the ratio is out of expected range (typical interval [0.65,2.6]) it tries to shift one of the timestamp columns (either drs or tel) to achieve better agreement
  // The ratio is calculated for <w> consecutive entires and the number of times it is within the acceptable range is counted by <n>
  // The relative offsets <o1> and <o2> are taken from the interval [0,nsteps) and the best set of offsets is one which yields highest <n>
  
  int w=10;		// number of consecutive events for inspection
  int n_best = -1;	// number of timestamps in agreement
  int n;		// used in loops
    
  int o0=0, o1=0, o0_best=0, o1_best=0;		

  if (verbose){
    switch (mode){
      case 1:
        cout  << ie << "\t" 
              << Event->ts_anchor << "\t" 
              << Event->delta_FEI4 << "\t" 
              << Event->delta_tel << "\t" 
              << 1.0*Event->delta_FEI4/TMath::Max(Event->delta_tel, (Long64_t)(1)) << "\t";
        break;
      default:
        cout  << ie << "\t" 
              << Event->ts_drs << "\t" 
              << Event->delta_drs*1e4 << "\t" 
              << Event->delta_tel << "\t" 
              << Event->delta_drs*1.0e4/TMath::Max(Event->delta_tel, (Long64_t)(1)) << "\t";
        break;
    }
  }
  // in the first loop try shifting both sets for a same offset, in the second loop also include relative offsets
  for (o0=0; o0<nsteps; o0++)
    for (o1=0; o1<nsteps; o1++){
      if (o0 != o1) continue;
      n=0;
      for (int i=0; i<w; i++){
        GetRatio(ie+o0+i+rel_offset, ie+o1+i, mode);
        if(RatioInRange(Event->ratio)) n++;
      }
      if (n>n_best){
        n_best = n;
        o0_best=o0;
        o1_best=o1;
      }
    }
//   cout << "A" << endl;
  for (o0=0; o0<nsteps; o0++)
    for (o1=0; o1<nsteps; o1++){
      if (o0 == o1) continue;
      n=0;
      for (int i=0; i<w; i++){
        GetRatio(ie+o0+i+rel_offset, ie+o1+i, mode);
        if(RatioInRange(Event->ratio)) n++;
      }
      if (n>n_best){
        n_best = n;
        o0_best=o0;
        o1_best=o1;
      }
    }
    
  if (verbose){
    cout  << n_best << "\t" 
          << o0_best << " " 
          << o1_best << " " 
          << rel_offset;
    if( o1_best - o0_best) 
    {
      cout << "\t incr. " <<  o0_best - o1_best;
      if(n_best < w-4) cout << " *";
    }
    cout << endl;
  }
        
  int orel, oabs;
  if(n_best < w-4){
    orel = 0;
    oabs = 3; // arbitrarily set to 3
  }
  else{
    orel = o0_best-o1_best;
    oabs = TMath::Min(o0_best, o1_best);
  }
//   cout << n_best << " " << oabs << " " << orel;
//   std::pair<int, int> offsets;
  
  return std::make_pair(orel, oabs);
};

void run::FindSignalWindow(int nsamples, int delta_nbins, int ch)
{
  // Search for the position of signal peak in first nsamples waveforms
  // delta_nbins = number of bins around peak in avgwf to sample for maximum in individual events
  // In subsequent analysis the pulse height will be sampled only in the interval [peak-delta, peak+delta]
  // Only channel "ch" is sampled. If there are large delays between separate channels, then the funciton needs to be modified
  waveform w;
  int delta_noise=2*delta_nbins + 30;     // time before the signal sampling window in which to sample for noise
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf->wf[ipt] = 0;    // reset average waveform;

  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break;
    GetEvent(ie);
    
    for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf->wf[ipt] += Event->wf[ch]->wf[ipt];
//     for(int ipt=0; ipt<DRS_N_POINTS; ipt++) w.wf[ipt] += Event->wf[ch]->wf[ipt];
  }
//   for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf->wf[ipt] /= (ie+1);    // rescale in the end;
//   SetSignalWindow(w.GetMaxPos(100)-delta_nbins, w.GetMaxPos(100)+delta_nbins);
  SetSignalWindow(ch, avwf->GetMaxPos()-delta_nbins, avwf->GetMaxPos()+delta_nbins);
  SetNoiseWindow(ch, avwf->GetMaxPos()-delta_noise-delta_nbins, avwf->GetMaxPos()-delta_noise+delta_nbins);
  
  return;
}

void run::FindSignalWindow(int ch, int nsamples, int delta_low, int delta_high, int pulse_polarity)
{
  // Search for the position of signal peak in first nsamples waveforms
  // delta_low/high = number of bins around peak in avgwf to sample for maximum in individual events
  // In subsequent analysis the pulse height will be sampled only in the interval [peak-delta_low, peak+delta_high]
  
//   waveform w;
  int delta_noise=(delta_low+delta_high) + 30;     // time before the signal sampling window in which to sample for noise
  for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf->wf[ipt] = 0;    // reset average waveform;

  for (int ie=0; ie<nsamples; ie++)
  {
    if (N <= ie) break;
    GetEvent(ie);
    for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf->wf[ipt] += Event->wf[ch]->wf[ipt];
//     for(int ipt=0; ipt<DRS_N_POINTS; ipt++) w.wf[ipt] += Event->wf[ch]->wf[ipt];
  }
//   for(int ipt=0; ipt<DRS_N_POINTS; ipt++) avwf->wf[ipt] /= (ie+1);    // rescale in the end;
//   SetSignalWindow(w.GetMaxPos()-delta_nbins, w.GetMaxPos()+delta_nbins);
  if(pulse_polarity > 0){
    SetSignalWindow(ch, avwf->GetMaxPos()-delta_low, avwf->GetMaxPos()+delta_high);
    SetNoiseWindow(ch, avwf->GetMaxPos()-delta_noise-delta_low, avwf->GetMaxPos()-delta_noise+delta_high);
  }
  else{
    SetSignalWindow(ch, avwf->GetMinPos()-delta_low, avwf->GetMinPos()+delta_high);
    SetNoiseWindow(ch, avwf->GetMinPos()-delta_noise-delta_low, avwf->GetMinPos()-delta_noise+delta_high);
  }
  return;
}

void run::SetSignalWindow(int ch, float min, float max)
{
  // Set time window in which to search for signal peak in waveforms
  sigw_min[ch]=TMath::Max(min,(float) 0.0);
  sigw_max[ch]=TMath::Min(max, (float)(DRS_N_POINTS-1.0));
//   sigw_min[ch]=0.0;
//   sigw_max[ch]=DRS_N_POINTS-1.0;
  return;
}

void run::SetNoiseWindow(int ch, float min, float max)
{
  // Set time window in which to sample noise
  noisew_min[ch]=TMath::Max(min, (float) 0.0);
  noisew_max[ch]=TMath::Min(max, (float)(DRS_N_POINTS-1.0));
  return;
}

float run::GetNoise(int ch, int mode)
{  
//   return Event->GetCharge(ch,noisew_min[ch],noisew_max[ch],mode) / (-noisew_min[ch]+noisew_max[ch]+1);

  //   return Event->GetCharge(ch,noisew_min[ch],noisew_max[ch],mode);
  float noise = Event->GetCharge(ch,noisew_min[ch],noisew_max[ch],mode);
//   float noise = Event->GetCharge(ch,noisew_min[ch],noisew_max[ch],1);
  for (int i=0; i<DRS_N_POINTS; i++) Event->yerr[i] = noise;
  
  return noise;
};

float run::FindThreshold(int ch, int nsamples, float sigma, int mode)
{
  // Finds discrimination threshold for channel #ch
  // The function samples noise level in first [nsamples] samples and fills it in a histogram
  // The obtained distribution is fitted with a gaussian
  // The threshold is set [sigma] standard deviations above the noise level
  
  TH1I hNoise("hnoise","hnoise",400, -0.04, 0.04);  // histogram for saving noise levels
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
//   cout << ch << " " << f.GetParameter(1) << " " << f.GetParameter(2) << " " << f.GetParameter(2)*sigma << endl;
  cout << "Threshold ch. " << ch << ": " << thr*1000 << " mV" << endl;
  return thr;
};

float run::GetEfficiency(int ch, float thr, int mode, int pulse_polarity)
{ 
  int j=Event->SearchTrack(ROI); // index of the first track going through the roi
//   int j=Event->SearchOneTrack(ROI); // index of the track going through the roi
  float charge=-1111;
  if (j>=0)
  {
    // additional track cuts
    double maxErrXY = 18;           // maximum allowable track origin error
    double maxSlopeErrXY = 1.6e-4;  // maximum allowable track slope error
    if ( (Event->errx[j] > maxErrXY) || (Event->erry[j] > maxErrXY)) return charge;
    if ( (Event->slopeErrX[j] > maxSlopeErrXY) || (Event->slopeErrY[j] > maxSlopeErrXY)) return charge;
    
//     float x = Event->x[j];
//     float y = Event->y[j];
    float x = Event->x[j] + zDUT[ch] * Event->slopeX[j];
    float y = Event->y[j] + zDUT[ch] * Event->slopeY[j];
    histos->heff_tracks[ch]->Fill(x,y);
//     charge = Event->GetCharge(ch,250,320,0);
//     charge = Event->GetCharge(ch,sigw_min[ch],sigw_max[ch],mode);
//     charge = Event->GetCharge(ch,sigw_min[ch],sigw_max[ch],mode) / (-sigw_min[ch]+sigw_max[ch]+1);
    charge = Event->GetCharge(ch,sigw_min[ch],sigw_max[ch],mode, pulse_polarity);
//     int binx = heff_charge->FindBin
//     heff_charge[ch]->SetBinContent(x,y, heff_charge[ch]->GetBinContent(x,y) + charge );
    histos->heff_charge[ch]->Fill(x,y,charge*1000);
    if (charge > thr)
      histos->heff_hits[ch]->Fill(x,y);
  }
  return charge;
};

float run::GetEfficiencyClustered(float *thr, int mode, int pulse_polarity)
{
  int j=Event->SearchTrack(ROI); // index of the first track going through the roi
  float charge=-1111;
  bool incluster = false;       // to fill at most one hit per event
  
  if (j>=0)
  {
    // additional track cuts
    double maxErrXY = 18;           // maximum allowable track origin error
    double maxSlopeErrXY = 1.6e-4;  // maximum allowable track slope error
    if ( (Event->errx[j] > maxErrXY) || (Event->erry[j] > maxErrXY)) return charge;
    if ( (Event->slopeErrX[j] > maxSlopeErrXY) || (Event->slopeErrY[j] > maxSlopeErrXY)) return charge;
//     
// //     float x = Event->x[j];
// //     float y = Event->y[j];
    float x = Event->x[j] + zDUT[0] * Event->slopeX[j];
    float y = Event->y[j] + zDUT[0] * Event->slopeY[j];
    histos->hclustered_tracks->Fill(x,y);
//     
    for (int ch=0; ch<nch; ch++){
      charge = Event->GetCharge(ch,sigw_min[ch],sigw_max[ch],mode, pulse_polarity);
      histos->hclustered_charge->Fill(x,y,charge*1000);
      if (charge > thr[ch])
        histos->hclustered_hits->Fill(x,y);
    }
  }
  return charge;
}

void run::SetROI(float xmin, float xmax, float ymin, float ymax, int nbinsx, int nbinsy)
{
  if (nbinsy<0) nbinsy=nbinsx;
  ROI = window(xmin,xmax,ymin,ymax);
  histos->SetHistosROI(ROI, nch, nbinsx, nbinsy);
//   for (int i=0; i<nch; i++)
//   {
// //     heff_tracks[i] = new TH2S(Form("heff_tracks_%d", i+1), Form("tracks through sample ch. %d ; x (#mum) ; y (#mum); nTracks",i+1),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
//     heff_tracks[i] = new TH2S(Form("heff_tracks_%d", i+1), Form("tracks through sample ch. %d ; x (#mum) ; y (#mum); nTracks",i+1),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
//     heff_hits[i] = new TH2S(Form("heff_hits_%d", i+1), Form("hits through sample ch. %d ; x (#mum) ; y (#mum); hits",i+1),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
//     heff_efficiency[i] = new TH2F(Form("heff_efficiency_%d", i+1), Form("Sample efficiency ch. %d ; x (#mum) ; y (#mum); efficiency",i+1),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
//     heff_charge[i] = new TH2F(Form("heff_charge_%d", i+1), Form("Pulse height ch. %d ; x (#mum) ; y (#mum); voltage (mV)",i+1),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
//     
//     heff_hits[i]->SetStats(0);
//     heff_tracks[i]->SetStats(0);
//     heff_efficiency[i]->SetStats(0);
//     heff_charge[i]->SetStats(0);
//   }
//   hTracks_ROI = new TH2S(Form("tracks ROI %s",fname_drs), "Tracks yielding pulses above threshold ; x (#mum); y (#mum); tracks", 50,xmin,xmax, 50,ymin,ymax);
//   
  dutcorr = DUTCorrelation (ROI, nbinsx, nbinsy, nch);
};

float run::GetResidualX(int it, int jt)
{
//   float res = Event->x[it] + (-300000)*Event->slopeX[it] - 250.0*Event->PixX[jt];
  float res = Event->x[it] - 250.0*Event->PosX[jt] + zAnchor*Event->slopeX[it];
//   return  0.1*Event->PixX[jt];
  return 1e-3*res;
};

float run::GetResidualY(int it, int jt)
{
//   float res = Event->y[it] + (-300000)*Event->slopeY[it] + 50.0*Event->PixY[jt];    // + sign because ROI module is mounted upside down
  float res = Event->y[it] + 50.0*Event->PosY[jt] + zAnchor*Event->slopeY[it];
  return 1e-3*res;
};

void run::FillResiduals()
{
  for (int it=0; it<Event->n; it++)
    for (int jt=0; jt<Event->n_a; jt++)
    {
//       if(Event->n > 1 || Event->n_a>1) continue;
      histos->hResX->Fill( GetResidualX(it,jt) );
      histos->hResY->Fill( GetResidualY(it,jt) );
      histos->hRes->Fill( GetResidualX(it,jt), GetResidualY(it,jt) );
    }
  return;
};

void run::FillResiduals(window anchor_ROI)
{
  for (int jt=0; jt<Event->n_a; jt++)
  {
    if ( !anchor_ROI.Contains(Event->PosX[jt], Event->PosY[jt])) continue;
    for (int it=0; it<Event->n; it++)
    {
//       if(Event->n > 1 || Event->n_a>1) continue;
      histos->hResX->Fill( GetResidualX(it,jt) );
      histos->hResY->Fill( GetResidualY(it,jt) );
      histos->hRes->Fill( GetResidualX(it,jt), GetResidualY(it,jt) );
    }
  }
  return;
};

bool run::CheckResiduals(int ie, int itr[2], double maxchi2)
{
  // Check if x,y residuals (between telescope track and FEI4 hit) are within the allowed window <ROI>. 
  // Only events with exactly one telescope track in the ROI are considered, to minimize the errors due to non-reconstructed tracks.
  // The function returns true if exactly one track with right residuals is found and false otherwise. If return true, the indices of the correcposnding tracks are saved in itr[0] (tele) and itr[1] (anchor).
  
  if (Event->NTracksInRegion(ROI) != 1) {
    return false;  // take only events with a single telescope track in ROI
  }
  int npassed=0;
  for (int it=0; it<Event->n; it++)
  {
    if( ! Event->TrackInRegion(it, ROI)) continue;  // Check if the track <it> goes through ROI
    for (int jt=0; jt<Event->n_a; jt++)
    {
      if (! ResCut.Contains(GetResidualX(it, jt), GetResidualY(it, jt))) continue;  // Check if the track passes the residual test
      if (npassed) {
//         cout << "FE" << npassed << "/" << Event->n_a << endl;
        return false;  // skip events with more than 1 matching hit in fei4
      }
      itr[0] = it;
      itr[1] = jt;
      npassed++;
      if (Event->chi2[it] > maxchi2) continue;
    }
  }
  if (npassed == 1) return false;
//   cout << "FE" << npassed << "/" << Event->n_a << endl;
  
  // HACK
//   GetEvent_s(ie - rel_offset - rel_offset_anchor -1);
//   int npassed=0;
//   for (int it=0; it<Event->n; it++)
//   {
//     if( ! Event->TrackInRegion(it, ROI)) continue;  // Check if the track <it> goes through ROI
//     for (int jt=0; jt<Event->n_a; jt++)
//     {
//       if (! ResCut.Contains(GetResidualX(it, jt), GetResidualY(it, jt))) continue;  // Check if the track passes the residual test
//       if (npassed) {
// //         cout << "FE" << npassed << "/" << Event->n_a << endl;
//         return false;  // skip events with more than 1 matching hit in fei4
//       }
//       itr[0] = it;
//       itr[1] = jt;
//       npassed++;
//     }
//   }
//   if (npassed == 1){
//     rel_offset_anchor--;
//     return true;
//   }
  // END OF HACK
  
  return false;     // Skip events with no matching hits in FEI4
  
};

bool run::CheckResiduals(int itr[2], double maxchi2)
{
  // Check if x,y residuals (between telescope track and FEI4 hit) are within the allowed window <ROI>. 
  // Only events with exactly one telescope track in the ROI are considered, to minimize the errors due to non-reconstructed tracks.
  // The function returns true if exactly one track with right residuals is found and false otherwise. If return true, the indices of the correcposnding tracks are saved in itr[0] (tele) and itr[1] (anchor).
  
  if (Event->NTracksInRegion(ROI) != 1) return false;  // take only events with a single telescope track in ROI
  int npassed=0;
  for (int it=0; it<Event->n; it++)
  {
    if( ! Event->TrackInRegion(it, ROI)) continue;  // Check if the track <it> goes through ROI
    for (int jt=0; jt<Event->n_a; jt++)
    {
      if (! ResCut.Contains(GetResidualX(it, jt), GetResidualY(it, jt))) continue;  // Check if the track passes the residual test
      if (npassed) return false;  // take only event with exactly 1 hit in fei4
      
      itr[0] = it;
      itr[1] = jt;
      npassed++;
      if (Event->chi2[it] > maxchi2) return false;;

    }
  }
  if (npassed == 1) return true;
  return false;
};

void run::FillAnalyzedEvent(int itr[2], float charge[NCH_DRS])
{
  IOHand->FillAnalyzedEvent(itr, charge, sigw_min);
  return;
};

int run::CreateSyncedFile()
{
  return IOHand->CreateSyncedFile();
}

int run::AddTelescopePlane(int n)
{
  return IOHand->AddTelescopePlane(n);
}

void run::FillSyncedFile(int ie)
{
  IOHand->FillSyncedFile(ie, rel_offset);
}

void run::FillSyncedFileInvalid()
{
  IOHand->FillSyncedFileInvalid();
}

void run::FillSyncedFileDummy(int ie)
{
  IOHand->FillSyncedFileDummy(ie, rel_offset);
}

void run::WriteAll()
{
  IOHand->f_res->Write();
}
