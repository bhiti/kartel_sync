{
  gErrorIgnoreLevel = kWarning;
  #define NCH_DRS 8
  #define DRS_N_POINTS 1024
  #define NPTS_ANAWF 150
  
  gStyle->SetPalette(55);
  TH1* dummy;   // root bs, script crashes if this is not present *here*

  #include "src/waveform.cpp"
  #include "src/event.cpp"
  #include "src/histograms.cpp"
  #include "src/DUTCorrelation.cxx"
  #include "src/IOHandler.cpp"
  #include "src/run.cpp"
  #include "PSTCT-macro.C"
    
  TH2I *hClusterSizeVsEff[NCH_DRS];
  for (int i=0; i<NCH_DRS;i++) hClusterSizeVsEff = new TH2I(Form("CSvsEff%d", i+1), Form("CSvsEff%d; cluster size; efficiency; Nentries", i+1), 11,-0.5,10.5, 100,0,1);
    
  char fname_drs[64], fname_tel[64], fname_anchor[64], fname_res[64], fname_ana[64];
  char brname_drs[32], lname_drs[32];
  char* run_number="087";
  TString plotsDir = Form("plots/run%s-10ns", run_number);
  gSystem->MakeDirectory(plotsDir.Data());
  int ch[NCH_DRS]={0,1,-2,-3,-4,-1,-1,-1};   // ch/drs_input: 0/1, 1/2, 2/3, 3/4
  int charge_mode=0;    // 0 ... current sensitive amp->waveforms will be integrated, 1 ... charge sensitive amplifier -> find steps in waveforms   
  int pulse_polarity = 1;   // EDIT HERE: positive or negative signals
  int Nevents = 3000000;     // EDIT HERE: number of events to process
  double maxchi2 = 2;
//   float zDut = 120000;      // DUT position along z-axis in micrometers (from designated 0)
    
  int verbose=0;

  bool using_anchor = true;
//   bool using_anchor = false;
  
  sprintf(fname_drs, Form("data/conv/run%s.root", run_number));         // DUT data
  sprintf(fname_tel, Form("data/ana/ref%s-p.root", run_number));        // Telescope file for synchronization with DUT
  sprintf(fname_anchor, Form("data/conv/anchor%s.root", run_number));   // Anchor module data (FEI4)
  sprintf(fname_res, Form("plots/results%s.root", run_number));         // File to save plots
  sprintf(fname_ana, Form("data/ana/ana%s.root", run_number));          // analized, matched and merged file
  
  sprintf(brname_drs, "Plane0/Waveforms");
//   sprintf(lname_drs, "%s%d", ch[ich]<2 ? "waveform" : "trigger", 1+ch[ich]%2); 
    
  vector<TString> v_lname_drs;
  for (int i=0; i<NCH_DRS; i++) if(ch[i]>=0) v_lname_drs.push_back(Form("%s%d", ch[i]<2 ? "waveform" : "trigger", 1+ch[i]%2));
//   for (int i=0; i<NCH_DRS; i++) if(ch[i]>=0) v_lname_drs.push_back(Form("waveform%d", ch[i]));
      
  ///////////////////////////////////
//   if (! using_anchor) run *r = new run(fname_drs, fname_tel, brname_drs, brname_drs, v_lname_drs);
//   else run *r = new run(fname_drs, fname_tel, fname_anchor, fname_res, brname_drs, brname_drs, v_lname_drs);
  if (! using_anchor) run *r = new run(fname_drs, fname_tel, brname_drs, brname_drs, v_lname_drs);
  else run *r = new run(fname_drs, fname_tel, fname_anchor, fname_res, fname_ana, brname_drs, brname_drs, v_lname_drs);

//   r->zDUT = zDut;
  
  /////////////////////////////////////
  
  int nbinsROI = 50;
//   r->SetROI(-700 ,1500,-800,1500,nbinsROI);    // 2x2 mm2   run 052
//   r->SetROI(-1500 ,-400,2100,2420,nbinsROI);    // edge   run 058
//   r->SetROI(-1200 ,50,2060,2500,50 );    // edge   run 075
//   r->SetROI(-2300 ,2000,-800,3500,nbinsROI);    // 4x4 mm2   run 069
//   r->SetROI(-1700 ,-300, 2000, 2500, nbinsROI);    // edge   run 058
    r->SetROI(-600,350, 1100, 2050, 50);    // 1x1   run 087
//   r->SetROI(-1500 ,3000,-1500,2800,nbinsROI);    // 4x4 mm2   run 090
//   r->SetROI(-500 ,2000,-500,2000,nbinsROI);    // 2x2 mm2   run 090
//   r->SetROI(-160 ,-40,1820,1950,nbinsROI);    // run 087 ch2 HE region
//   r->SetROI(-300 ,-150,1600,1700,nbinsROI);    // run 087 ch2 LE region
//   r->SetROI(-1450,-1300, 2450, 2600,nbinsROI);    // CMOS   run 097

  window fiducial_region[NCH_DRS];    // used in the plotting script
//   fiducial_region[0].Set(-1900,1750,-400,3100); fiducial_region[1].Set(-1000,800,400,2200); // run 084
//   fiducial_region[0].Set(-1100,2500,-1100,2400); fiducial_region[1].Set(-100,1600,-100,1600); // run 090
  fiducial_region[0].Set(-600,350, 1100, 2050); fiducial_region[1].Set(-600,350, 1100, 2050); // run 087
//   fiducial_region[0].Set(-1375,-1325, 2505, 2555); fiducial_region[1].Set(-1425,-1375, 2505, 2555); // run 097

  if (using_anchor) r->ResCut = window(-9.86, -9.61, 7.0, 7.05);    //window(-9.75, -9.35, 6.9, 7.05);
  window anchor_ROI(15,50,50,250);      // mask out noisy pixels in anchor module. Pixels outside of this window will be ignored
  
  // regions with low charge / high charge
  window wlc(600,700,950,1150);
  window whc(700,800,1300,1400);
  waveform wflc;
  waveform wfhc;
  int nlc=0, nhc=0;
  
  float thr[NCH_DRS];     // Signal thresholds for efficiency measurement
  int istart = 2;
//   int istart = 290000;
//   r->rel_offset = -336;
  
  int passed=0;
    
  // Get the correct time window for pulse height measurement
//   int delta_low=5;                        // number of bins before signal peak to sample - should be >0
//   int delta_high=10;                      // number of bins after signal peak to sample - should be >0
  int delta_low=4;                        // number of bins before signal peak to sample - should be >0
  int delta_high=6;                      // number of bins after signal peak to sample - should be >0
  float delta=delta_high+delta_low+1.;    // total width of integration window
  float sigma_thr=4.0;
  for (int i=0; i<r->nch; i++) 
  {
    if (charge_mode==0)
      r->FindSignalWindow(i, 1000, delta_low,delta_high, pulse_polarity);  // (ch, nsamples, delta_low, delta_high)
    else {
//       r->SetSignalWindow(i,350, 450);
//       r->SetNoiseWindow(i,100., 200.);
      r->SetSignalWindow(i,370, 450);
      r->SetNoiseWindow(i,100., 180.);
    }
  }
  for (int i=0; i<r->nch; i++) 
  {
    thr[i] = r->FindThreshold(i,1000, sigma_thr, charge_mode);   // FindThreshold(int ch, int nsamples=1000, float sigma=4.). After signal and noise sampling windows are determined, determine noise level. The window for noise sampling is defined in FindSignalWindow
//     thr[i] = 0.004;
  }
  
//   TCanvas *c1 = new TCanvas("c1", "c1", 0,700,700,500);
//   TCanvas *c2 = new TCanvas("c2", "c2", 1);

  for (int i=istart;i<Nevents; i++) 
  {
    if (r->N < i) break;    
    if (i%10000 == 0) cout << i << " " << r->rel_offset << " " << passed << endl;
    if (strncmp(run_number, "084",3 )==0) {
      if (i==264000){
        i=266000;
        r->rel_offset = -300;
        cout << "HAND FIX 264000 --> 266000" << endl;
      }
      if (i==578000){
        i=582000;
        r->rel_offset = -686;
        cout << "HAND FIX 578000 --> 582000" << endl;
      }
      if (i==404000){
        i=405000;
        r->rel_offset = -477;
        cout << "HAND FIX 404000 --> 405000" << endl;
      }
      if (i==842000){
        i=845000;
        r->rel_offset = -1000;
        cout << "HAND FIX 842000 --> 845000" << endl;
      }
      if (i==918000){
        i=920000;
        r->rel_offset = -1100;
        cout << "HAND FIX 918000 --> 920000" << endl;
      }
    }
    
    if (strncmp(run_number, "056",3 )==0) {
      if (i==370000){
        i=371000;
        r->rel_offset = -452;
        cout << "HAND FIX: 370000 --> 371000" << endl;
      }
    }
    
    if (strncmp(run_number, "087",3 )==0) {
      if (i==840000){
        i=842000;
        r->rel_offset = -1004;
        cout << "HAND FIX: 840000 --> 842000" << endl;
      }
      if (i==848000){
        i=849000;
        r->rel_offset = -1010;
        cout << "HAND FIX: 848000 --> 849000" << endl;
      }
    }
    
    if (strncmp(run_number, "090",3 )==0) {
      if (i==303000){
        i=304000;
        r->rel_offset = -354;
        cout << "HAND FIX: 303000 --> 304000" << endl;
      }
    }
    
    r->ProcessTimestamp(i,verbose);    
    r->GetEvent(i);
    
//     if (r->Event.VetoPresent(2, -0.2)) continue;
//     if (r->Event.VetoPresent(2, 0.12,2)) continue;
//     if (r->Event.VetoPresent(2, 0.02,1)) continue;

    
//     cout << i << endl;
    
//     for (int it=0; it<r->Event.n; it++)
//     {
//       if( r->Event.TrackInRegion(it, whc)){
//         // Check if the track <it> goes through ROI
//         if(using_anchor)
//         {
//           r->GetEvent_s(i);  // Get anchor module data
//           r->FillResiduals(anchor_ROI);
//           int itr[2]; itr[0]=-1; itr[1]=-1;   // container for telescope / FEI4 track index through ROI
//           if (! r->CheckResiduals(itr)) continue;
//           for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfhc.wf[ipt] += r->Event.wf[1].wf[ipt];
//           nhc++;
//         }
//       }
//       else if( r->Event.TrackInRegion(it, wlc)){
//         if(using_anchor)
//         {
//           r->GetEvent_s(i);  // Get anchor module data
//           r->FillResiduals(anchor_ROI);
//           int itr[2]; itr[0]=-1; itr[1]=-1;   // container for telescope / FEI4 track index through ROI
//           if (! r->CheckResiduals(itr)) continue;
//           for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wflc.wf[ipt] += r->Event.wf[1].wf[ipt];
//           nlc++;
//         }
//       }
      /*
      if (1)
      {
        r->Event.Draw(1);
        r->Event.gr[1]->GetXaxis()->SetRangeUser(0,1111);
        r->Event.gr[1]->GetYaxis()->SetRangeUser(-0.4,0.06);
//         r->Event.Draw(1,"PL"); 
//         r->Event.Draw(2,"PL"); 
//         r->Event.Draw(3,"PL"); 
        gPad->Update();
  //       gPad->Print(Form("plots/pulse%d.png",i));
  //       cout << charge << endl;
        char a;
        cin >> a;      
      } 
      //*/
//     }
//     continue;
    
    
    /////////////////////////////////
    if(using_anchor)
    {
      r->GetEvent_s(i);  // Get anchor module data
      r->FillResiduals(anchor_ROI);
      int itr[2]; itr[0]=-1; itr[1]=-1;   // container for telescope / FEI4 track index through ROI
      if (! r->CheckResiduals(itr, maxchi2)) continue;
    }
    
    passed++;
    //////////////////////////////////
    float charge[NCH_DRS];    
    for (int ich=0; ich< r->nch; ich++)
    {
      r->histos->hN[ich]->Fill(r->GetNoise(ich, charge_mode));
      charge[ich] = r->GetEfficiency(ich, thr[ich], charge_mode, pulse_polarity);
      r->histos->hC[ich]->Fill(charge[ich]);
    }    

    r->dutcorr.ProcessEvent(r->Event, thr);   // Get correlation between DUT channels
    r->FillAnalyzedEvent(itr, charge);
    
    for (int j=0; j<r->Event.n; j++) 
    {
      r->histos->hTracks->Fill(r->Event.x[j], r->Event.y[j]);
      r->histos->hTracks_ROI->Fill(r->Event.x[j], r->Event.y[j]);	
//       r->histos->hClusterSize->Fill(r->Event.x[j], r->Event.y[j], r->Event.);	
    }
    
    /*
//     if (((charge[0] > 0.5*thr[0]) && (charge[0] < thr[0])) || ((charge[1] > 0.5*thr[1]) && (charge[1] < thr[1])))
    if (charge[0] > -1110.029)
    {
      c1->cd();
      TH2S ht("ht", "ht", nbinsROI, -1450,-1300, nbinsROI, 2450, 2600);
      ht.Fill(r->Event.x[itr[0]], r->Event.y[itr[0]]);
      ht->SetStats(0);
      ht.Draw("COLZ");
//       r->hTracks_ROI->Draw("COLZ");
//       gPad->Update();
      TBox fb[2];
      TLatex ttex; ttex.SetTextFont(43); ttex.SetTextSize(20); // tex.SetTextColor(2);
      for (int ich=0; ich<r->nch-1; ich++){
        fb[ich] = TBox(fiducial_region[ich].x1, fiducial_region[ich].y1, fiducial_region[ich].x2, fiducial_region[ich].y2);
        fb[ich].SetFillStyle(0);
        fb[ich].SetLineWidth(2);
        fb[ich].SetLineColor(1+ich);
        fb[ich].Draw();
      }
      ttex.DrawLatexNDC(0.4,0.25, Form("#splitline{Ch1 = %.1lf mV}{Ch2 = %.1lf mV}", 1000*charge[0], 1000*charge[1])); 
      gPad->Update();

      c2->cd();
      r->Event.Draw(0);
      r->Event.gr[0]->GetXaxis()->SetRangeUser(300,600);
      r->Event.gr[0]->GetYaxis()->SetRangeUser(-0.03,0.1);
//       r->Event.gr[0]->GetYaxis()->SetRangeUser(-0.55,0.55);
      r->Event.Draw(1,"PL"); 
      r->Event.Draw(2,"PL"); 
      r->Event.Draw(3,"PL"); 
      gPad->Update();
//       gPad->Print(Form("plots/pulse%d.png",i));
//       cout << charge << endl;
      char a;
      cin >> a;    
    }//*/    
  }
  
  r->IOHand->f_ana->Write();
  
  
  //
  // Generate and save plots
  //
  
  gROOT->ProcessLine(".x plotting_script_cmos.C");
   
  r->dutcorr.AddWeight(r->histos->hTracks);
  for (int i=0; i<r->nch; i++)
    for (int j=0; j<r->nch; j++){
      cout << "#hits ratio " << i << j << ": " << r->dutcorr.GetCorrelation(j,i) << endl;
    }
    
  int ctr=0;
  float wx = r->histos->heff_efficiency[0]->GetXaxis()->GetBinWidth(5);
  float wy = r->histos->heff_efficiency[0]->GetYaxis()->GetBinWidth(5);
  for (int ich=0; ich<r->nch; ich++){
    ctr=0;
    for (int ix=1; ix<=r->histos->heff_efficiency[ich]->GetNbinsX(); ix++)
      for (int iy=1; iy<=r->histos->heff_efficiency[ich]->GetNbinsY(); iy++){
        if(r->histos->heff_efficiency[ich]->GetBinContent(ix,iy) > 0.5) ctr++;
      }
    cout << ich << " ... " << ctr << " (" << ctr*wx*wy << " um2)" << endl;
  }
//       cout << "#hits ratio 21: " << r->dutcorr.GetCorrelation(0,1) << endl;
//   cout << "#hits ratio 11: " << r->dutcorr.GetCorrelation(0,0) << endl;
  
  
//   if (nlc > 0) for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wflc.wf[ipt] /= nlc;
//   if (nhc > 0) for(int ipt=0; ipt<DRS_N_POINTS; ipt++) wfhc.wf[ipt] /= nhc;
//   
//   glc = wflc.Draw();
//   ghc = wfhc.Draw();
//   
//   glc->Draw("APL");
//   ghc->Draw("PL");
  
}