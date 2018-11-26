IOHandler::IOHandler(){};

IOHandler::~IOHandler()
{
  f_drs->Close();
  delete f_drs; 
  delete t_drs_wf; 
  delete t_drs_ts; 
  f_tel->Close();
  delete f_tel; 
  delete t_tel_ts; 
  f_anchor->Close();
  f_res->Close();
  delete f_anchor; 
  delete t_anchor_hits; 
  f_ana->Close();
  delete f_ana;
  f_res->Close();
  delete f_res;
  
//   delete outputEvents;
//   delete outputHitTrees;
//   delete inputHitTrees;
};  

void IOHandler::Init(event* Event, const char* fname_drs, const char* fname_tel, const char* fname_anchor, const char* fname_results, const char* fname_ana, const char *brname_drs)
{
  _event = Event;
  _anaevent = new analyzed_event();
//   nch = lname_drs.size();
  if (strlen(fname_drs) != 0){
    f_drs = new TFile(fname_drs);
//     for (int i=0; i<nch; i++) {
//       v_t_drs_wf[i]= (TTree*) f_drs->Get(brname_drs);
//       v_t_drs_wf[i]->SetBranchAddress(Form("%s",lname_drs[i].Data()), _event->wf[i].wf);
//     }
    t_drs_ts = (TTree*) f_drs->Get("Event");
    t_drs_ts->SetBranchStatus("TimeStamp",1);
    t_drs_ts->SetBranchAddress("TimeStamp", &_event->ts_drs);
  }
  else cout << "No DRS file specified" << endl;
  
  if (strlen(fname_tel) != 0){
    cout << "Loading telescope file " << fname_tel << endl;
    f_tel = new TFile(fname_tel);	
    t_tel_ts = (TTree*) f_tel->Get("Event");
    t_tel_ts->SetBranchStatus("TimeStamp",1);
    t_tel_ts->SetBranchAddress("TimeStamp", &_event->ts_tel);
    
//     t_tel_tracks = (TTree*) f_tel->Get("Tracks");
//     t_tel_tracks->SetBranchStatus("Origin*",1);
//     t_tel_tracks->SetBranchAddress("OriginX", &_event->x);
//     t_tel_tracks->SetBranchAddress("OriginY", &_event->y);
//     t_tel_tracks->SetBranchAddress("OriginErrX", &_event->errx);
//     t_tel_tracks->SetBranchAddress("OriginErrY", &_event->erry);
//     t_tel_tracks->SetBranchAddress("NTracks", &_event->n);
//     t_tel_tracks->SetBranchAddress("SlopeX", &_event->slopeX);
//     t_tel_tracks->SetBranchAddress("SlopeY", &_event->slopeY);
//     t_tel_tracks->SetBranchAddress("SlopeErrX", &_event->slopeErrX);
//     t_tel_tracks->SetBranchAddress("SlopeErrY", &_event->slopeErrY);
//     t_tel_tracks->SetBranchAddress("Chi2", &_event->chi2);
  }
  else cout << "No telescope file specified" << endl;
  
  if (strlen(fname_anchor) != 0){
    cout << "Loading FEI4 file      " << fname_anchor << endl;
    f_anchor = new TFile(fname_anchor);	
    t_anchor_ts = (TTree*) f_anchor->Get("Event");
    t_anchor_ts->SetBranchStatus("TimeStamp",1);
    t_anchor_ts->SetBranchAddress("TimeStamp", &_event->ts_anchor);
    
    t_anchor_hits = (TTree*) f_anchor->Get("Plane0/Hits");
//   //   t_anchor_hits->SetBranchAddress("PixX", &_event->PixX);
//   //   t_anchor_hits->SetBranchAddress("PixY", &_event->PixY);
//     t_anchor_hits->SetBranchAddress("PosX", &_event->PosX);
//     t_anchor_hits->SetBranchAddress("PosY", &_event->PosY);
//     t_anchor_hits->SetBranchAddress("NClusters", &_event->n_a);
//     t_anchor_hits->SetBranchAddress("ClusterSize", &_event->ClusterSize);
//     t_anchor_hits->SetBranchAddress("Value", &_event->ClusterCharge);

//     t_anchor_hits->SetBranchAddress("PixX", &_event->PixX);
//     t_anchor_hits->SetBranchAddress("PixY", &_event->PixY);
//     t_anchor_hits->SetBranchAddress("NHits", &_event->n_a);
//     t_anchor_hits->SetBranchAddress("Timing", &_event->ClusterSize);
//     t_anchor_hits->SetBranchAddress("InCluster", &_event->ClusterSize);
//     t_anchor_hits->SetBranchAddress("Value", &_event->ClusterCharge);
  }
  else cout << "No FEI4 anchor file specified" << endl;
      
  if (strlen(fname_results) != 0){
    f_res= new TFile(fname_results, "RECREATE");
  }
  else cout << "No plots output file specified" << endl;
      
  if (strlen(fname_ana) != 0){
    cout << "Creating output file " << fname_ana << endl;
    f_ana= new TFile(fname_ana, "RECREATE");
  
    for (int i=0; i<nch; i++){
      t_ana[i] = new TTree(Form("ch%d", i+1), Form("ch%d", i+1));
      t_ana[i]->Branch("waveform", _anaevent->wf[i], TString::Format("wf[%i]/F", NPTS_ANAWF));
      t_ana[i]->Branch("charge10ns", &_anaevent->charge[i], "charge/F");
    }
    t_ana_tracks = new TTree("Tracks","Tracks");
    t_ana_tracks->Branch("xt", &_anaevent->xt, "xt/D");
    t_ana_tracks->Branch("yt", &_anaevent->yt, "yt/D");
    t_ana_tracks->Branch("chi2", &_anaevent->chi2, "chi2/D");
    t_ana_tracks->Branch("NTracks", &_anaevent->nt, "NTracks/I");
    t_ana_tracks->Branch("slopeX", &_anaevent->slopeX, "slopeX/D");
    t_ana_tracks->Branch("slopeY", &_anaevent->slopeY, "slopeY/D");
    t_ana_tracks->Branch("xa", &_anaevent->xa, "xa/D");
    t_ana_tracks->Branch("ya", &_anaevent->ya, "ya/D");    
    t_ana_tracks->Branch("NClusters_a", &_anaevent->na, "nclusters_a/I");    
    t_ana_tracks->Branch("ClusterSize_a", &_anaevent->clustersize, "clustersize_a/s");    
  }
  else cout << "No analyzed output file specified" << endl;
  
  if (f_res != 0) gDirectory->cd(Form("%s:", fname_results));  
  
};

int IOHandler::CreateSyncedFile()
{
  cout << "Creating synced file " << f_res->GetName() << endl;

  //#################################
  //load input file and link branches
  //#################################
  TFile * inputFile = f_tel;
  if(!inputFile){
    std::cout << __PRETTY_FUNCTION__ << " :: Input file could not be opened!" << std::endl;
    return -1;
  }

  TTree * inputEvents = t_tel_ts;
//   if(inputEvents->SetBranchAddress("TimeStamp", &_event->tel_ts) != 0){
//     std::cout << "Branch not found, stopping script..." << std::endl;
//     return -1;
//   }
  if(inputEvents->SetBranchAddress("FrameNumber", &_event->tel_frameNumber) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(inputEvents->SetBranchAddress("TriggerInfo", &_event->tel_triggerInfo) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(inputEvents->SetBranchAddress("TriggerOffset", &_event->tel_triggerOffset) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(inputEvents->SetBranchAddress("Invalid", &_event->tel_invalid) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }

  //set up plane dirs and hit trees to read in
  // telescope planes
  for(UInt_t i = 0; i < planes.size(); i++){
    TDirectory * planeDir = new TDirectory();
    std::string name = "Plane" + std::to_string(planes[i]);
    f_tel->GetObject(name.c_str(), planeDir);
    if(!planeDir){
      std::cout << __PRETTY_FUNCTION__ << ": Plane" << planes[i] << " directory not found!" << std::endl;
      return -1;
    }
    inputPlaneDirs.push_back(planeDir);

    TTree * hitsTree = new TTree();
    inputPlaneDirs[i]->GetObject("Hits",hitsTree);
    if(!hitsTree){
      std::cout << __PRETTY_FUNCTION__ << ": Hits tree not found!" << std::endl;
      return -1;
    }
    
    if(hitsTree->SetBranchAddress("NHits", &_event->hits_NHits[i]) != 0){
      std::cout << "Branch not found, stopping script..." << std::endl;
      return -1;
    }
    if(hitsTree->SetBranchAddress("Value", &_event->hits_Value[i]) != 0){
      std::cout << "Branch not found, stopping script..." << std::endl;
      return -1;
    }
    if(hitsTree->SetBranchAddress("Timing", &_event->hits_Timing[i]) != 0){
      std::cout << "Branch not found, stopping script..." << std::endl;
      return -1;
    }
    if(hitsTree->SetBranchAddress("PixX", &_event->hits_PixX[i]) != 0){
      std::cout << "Branch not found, stopping script..." << std::endl;
      return -1;
    }
    if(hitsTree->SetBranchAddress("PixY", &_event->hits_PixY[i]) != 0){
      std::cout << "Branch not found, stopping script..." << std::endl;
      return -1;
    }
    if(hitsTree->SetBranchAddress("InCluster", &_event->hits_HitInCluster[i]) != 0){
      std::cout << "Branch not found, stopping script..." << std::endl;
      return -1;
    }
    
    inputHitTrees.push_back(hitsTree);
  }
  
  // FEI4 plane
  planes.push_back(0);
  int k = planes.size()-1;
  std::string name = "Plane0";
  TDirectory * planeDir = new TDirectory();
  f_anchor->GetObject(name.c_str(), planeDir);
  if(!planeDir){
    std::cout << __PRETTY_FUNCTION__ << ": FEI4 Plane0 directory not found!" << std::endl;
    return -1;
  }
  inputPlaneDirs.push_back(planeDir);  
  
  TTree * hitsTree = new TTree();
  inputPlaneDirs.back()->GetObject("Hits",hitsTree);
  if(!hitsTree){
    std::cout << __PRETTY_FUNCTION__ << ": Hits tree not found!" << std::endl;
    return -1;
  }
  if(hitsTree->SetBranchAddress("NHits", &_event->hits_NHits[k]) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(hitsTree->SetBranchAddress("Value", &_event->hits_Tot[k]) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(hitsTree->SetBranchAddress("Timing", &_event->hits_Timing[k]) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(hitsTree->SetBranchAddress("PixX", &_event->hits_PixX[k]) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(hitsTree->SetBranchAddress("PixY", &_event->hits_PixY[k]) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  if(hitsTree->SetBranchAddress("InCluster", &_event->hits_HitInCluster[k]) != 0){
    std::cout << "Branch not found, stopping script..." << std::endl;
    return -1;
  }
  
  inputHitTrees.push_back(hitsTree);

  //##################
  //set up output file
  //##################
  TFile * outputFile = f_res;
  
//   TTree * outputEvents = new TTree("Event","Event Information");
  outputEvents = new TTree("Event","Event Information");
  if(!outputEvents){
    std::cout << "Could not open output event tree!" << std::endl;
    return -1;
  }
  outputEvents->Branch("FrameNumber", &_event->tel_frameNumber, "FrameNumber/l");
  outputEvents->Branch("TimeStamp", &_event->ts_tel, "TimeStamp/l");
  outputEvents->Branch("TriggerInfo", &_event->tel_triggerInfo, "TriggerInfo/I");
  outputEvents->Branch("TriggerOffset", &_event->tel_triggerOffset, "TriggerOffset/I");
  outputEvents->Branch("Invalid", &_event->tel_invalid, "Invalid/O");

  //link hit tree branches to variables
  for(UInt_t k = 0; k < planes.size(); k++){
    std::string name = "Plane" + std::to_string(k);
    outputPlaneDirs.push_back(outputFile->mkdir(name.c_str()));
    outputPlaneDirs[k]->cd();
    TTree * hits = new TTree("Hits","Hits");
    if(!hits){
      std::cout << "Could not create hits tree for plane " << name << "!" << std::endl;
      return -1;
    }
    outputHitTrees.push_back(hits);
    
    hits->Branch("NHits", &_event->hits_NHits[k], "NHits/I"); //event->size()
    hits->Branch("Value", &_event->hits_Value[k], "Value[NHits]/D"); //(*event)[i].ToT
    hits->Branch("Timing", &_event->hits_Timing[k], "Timing[NHits]/D"); //(*event)[i].time1_56ns
    hits->Branch("PixX", &_event->hits_PixX[k], "PixX[NHits]/I"); //(*event)[i].X
    hits->Branch("PixY", &_event->hits_PixY[k], "PixY[NHits]/I"); //(*event)[i].Y
    hits->Branch("HitInCluster", &_event->hits_HitInCluster[k], "HitInCluster[NHits]/I");
  }
    
  return 0;
}

void IOHandler::GetEvent(int ie, int rel_offset){
  for(int i=0; i<nch;i++) v_t_drs_wf[i]->GetEntry(ie - rel_offset);
  t_tel_tracks->GetEntry(ie);
  return; 
}

void IOHandler::GetEvent_s(int ie_s, int rel_offset)
{
  // Get data from files synchronized with the anchor module
  t_anchor_ts   ->GetEntry(ie_s - rel_offset);
  t_anchor_hits ->GetEntry(ie_s - rel_offset);
  return;
};

void IOHandler::GetTimeStamp(int ie)
{
  //t_drs_ts->GetEntry(ie);
  GetDrsTimestamp(ie);
  //t_tel_ts->GetEntry(ie);
  GetTelTimestamp(ie);
  return;
};

Long64_t IOHandler::GetDrsTimestamp(int ie) 
{
  if (ie<0 || ie > t_drs_ts->GetEntries()) return 0;
  
  t_drs_ts->GetEntry(ie);
  Long64_t _ts = _event->ts_drs;
  Long64_t _tsold = 0;
    
  // Check for 1e4 bug (sometimes - rarely - the 5th last digit of the timestamp has to be incremented by 1
  if (ie>0)
  {
    t_drs_ts->GetEntry(ie-1);
    _tsold = _event->ts_drs;
    if(_ts < _tsold) _ts += 1e4;		// bug fix : drs timestamps should be strictly increasing, but seems like a bug, which sometimes does not increment 5th least significant digit in time. This line seems to fix this
    _event->ts_drs = _ts;
  }
  return _ts;
};

Long64_t IOHandler::GetDrsDelta(int ie) 
{
  Long64_t _delta =  - (GetDrsTimestamp(ie-1) - GetDrsTimestamp(ie));
  _event->delta_drs = _delta;
  return _delta;
};

Long64_t IOHandler::GetTelTimestamp(int ie)
{
  if (ie<0 || ie > t_tel_ts->GetEntries()) 
  {
    //cout << "ts 0" << endl;
    return 0;
  }
  t_tel_ts->GetEntry(ie);
//   t_tel_trgoffset->GetEntry(ie);
  return _event->ts_tel;
//   return _event->ts_tel - _event->tel_trgoffset;
};

Long64_t IOHandler::GetTelDelta(int ie) 
{
  Long64_t _tsold = GetTelTimestamp(ie-1);
  Long64_t _ts = GetTelTimestamp(ie);
 
  if (_ts < _tsold) _ts += (ULong64_t)(pow(2,32));	// correct counter overflow
  
  Long64_t _delta = _ts - _tsold;
  _event->delta_tel = _delta;
  return _delta;
};

Long64_t IOHandler::GetFEI4Timestamp(int ie)
{
  if (ie<0 || ie > t_anchor_ts->GetEntries()) {
    return 0;
  }
  t_anchor_ts->GetEntry(ie);
//  cout << _event->ts_anchor << endl;
  return (Long64_t)(_event->ts_anchor);
};

Long64_t IOHandler::GetFEI4Delta(int ie) 
{
  // FEI4 overflows at pow(2,31) approx. every 50 seconds (40 MHz timestamp clock)
  Long64_t _tsold = GetFEI4Timestamp(ie-1);
  Long64_t _ts = GetFEI4Timestamp(ie);
  
  if (_ts < _tsold) _ts += (Long64_t)(pow(2,32));	// correct counter overflow
  
  Long64_t _delta = _ts - _tsold;
  _event->delta_FEI4 = _delta;
  return _delta;
};

void IOHandler::FillAnalyzedEvent(int itr[2], float charge[NCH_DRS], float sigw_min[NCH_DRS])
{
  for (int ch=0; ch<nch; ch++){
    int pt_start = sigw_min[ch]-40;   // estimate for the start of the waveform
    for (int i=0; i<NPTS_ANAWF; i++) _anaevent->wf[ch][i] = _event->wf[ch]->wf[pt_start+i];
    _anaevent->charge[ch] = charge[ch];
    t_ana[ch]->Fill();
  }
  _anaevent->xt = _event->x[itr[0]];
  _anaevent->yt = _event->y[itr[0]];
  _anaevent->slopeX = _event->slopeX[itr[0]];
  _anaevent->slopeY =_event->slopeY[itr[0]];
  _anaevent->chi2 =_event->chi2[itr[0]];
  _anaevent->nt =_event->n;
  _anaevent->xa =_event->PosX[itr[1]];
  _anaevent->ya =_event->PosY[itr[1]];
  _anaevent->clustersize =_event->ClusterSize[itr[1]];
  _anaevent->na =_event->n_a;
  t_ana_tracks->Fill();
  return;
};

int IOHandler::AddTelescopePlane(int n)
{
  planes.push_back(n);
  return planes.size();
}

void IOHandler::FillSyncedFile(int ie, int rel_offset)
{
  t_tel_ts->GetEntry(ie - rel_offset);
  
  for (UInt_t k=0; k<planes.size(); k++){
    if (k<planes.size()-1) inputHitTrees[k]->GetEntry(ie - rel_offset);
    else{
      inputHitTrees[k]->GetEntry(ie);
      for (int j=0; j<_event->hits_NHits[k]; j++) _event->hits_Value[k][j] = (Double_t)(_event->hits_Tot[k][j]);
    }
    outputHitTrees[k]->Fill();
//     if (k==6 && _event->hits_NHits[k]>0) cout << _event->hits_Value[k][0]<< endl;
  }
  outputEvents->Fill();
}

void IOHandler::FillSyncedFileInvalid()
{  
  for (UInt_t k=0; k<planes.size(); k++){
//     inputHitTrees[k]->GetEntry(ie);
    outputHitTrees[k]->Fill();
  }
  _event->tel_invalid = true;
  outputEvents->Fill();
}

void IOHandler::FillSyncedFileDummy(int ie, int rel_offset)
{
  for (UInt_t k=0; k<planes.size(); k++){
    if (k<planes.size()-1) inputHitTrees[k]->GetEntry(ie - rel_offset);
    else inputHitTrees[k]->GetEntry(ie);
    _event->hits_NHits[k]=0;
    outputHitTrees[k]->Fill();
  }
  t_tel_ts->GetEntry(ie - rel_offset);
  outputEvents->Fill();
}

