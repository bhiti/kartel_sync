{
  cout << "plotting" << endl;
  TCanvas *c1 = new TCanvas("c1", "c1", 1);
  c1->SetMargin(0.12,0.05,0.12,0.02);
  TLatex tex; tex.SetTextFont(43); tex.SetTextSize(20); // tex.SetTextColor(2);
  
//   c1->SetLogy();
  TF1* fgaus[NCH_DRS];
  TF1* flangau[NCH_DRS];
  TF1* fn[NCH_DRS];
  
//   TFile f("spectra.root","RECREATE");    
//   TH1I *hs[2];
//   for (int i=0; i< r->nch; i++)
//   {
//     hs[i] = (TH1I*) r->histos->hC[i]->Clone(Form("hs[%d]",i));
//     DrawTH1 (hs[i]);
//   }
//   f.Write();
//   f.Close();

  for (int i=0; i<r->nch; i++){
//     if (fgaus[i]) delete fgaus[i];
//     if (flangau[i]) delete flangau[i];
//     if (fn[i]) delete fn[i];

//     gPad->SetLogy(1);
//     r->histos->hC[i]->GetYaxis()->SetRangeUser(0, 150);
//     r->histos->hC[i]->Draw();
//     gPad->Update();
    DrawTH1 (r->histos->hC[i]);
    fgaus[i] = new TF1(Form("gaus%d",i),"gaus", -0.05, 0.5*thr[i]);
//     fgaus[i]->SetParameters(1,1,1);
//     fgaus[i] = new TF1(Form("gaus%d",i),"gaus", -0.5*thr[i], 0.5*thr[i]);

    r->histos->hC[i]->Fit(Form("gaus%d",i),"RQ");
    flangau[i] = new TF1(Form("langau%d",i),langaufun, 1.8*thr[i], r->histos->hC[0]->GetXaxis()->GetXmax(), 4);
    flangau[i]->SetParameters(2e-3 , 0.015, 1, 2e-3);  // scale, MPV, intgral, gaus_sigma
    flangau[i]->FixParameter(3,fgaus[i]->GetParameter(2));
    r->histos->hC[i]->Fit(Form("langau%d",i),"RQ");
    r->histos->hC[i]->Draw();
    fgaus[i]->Draw("SAME");
    TLine tl(thr[i], 0, thr[i], r->histos->hC[i]->GetMaximum());
    tl.SetLineWidth(2); tl.SetLineColor(4);
    tl.Draw();
    
    tex.DrawLatexNDC(0.43,0.85, Form("#splitline{Pedestal (#sigma) = %.1lf mV}{Thr (%.0lf#sigma) = %.1lf mV}", 1000*fgaus[i]->GetParameter(2), sigma_thr, thr[i]*1000)); 
    if (flangau[i]->GetNDF() > 0){
      tex.DrawLatexNDC(0.7,0.35, Form("#splitline{#splitline{#splitline{Langau fit:}{MPV = %.0lf mV}}{#sigma_{gaus} = %.1lf mV}}{#Chi^2/NDF = %.1lf}", 1000*flangau[i]->GetParameter(1), 1000*flangau[i]->GetParameter(3), flangau[i]->GetChisquare()/flangau[i]->GetNDF())); 
    }
    c1->Print(Form("%s/spectrum_ch%d.pdf",plotsDir.Data(), i+1));
    c1->Print(Form("%s/spectrum_ch%d.png",plotsDir.Data(), i+1));
    // Draw Noise
//     gPad->SetLogy(0);
    DrawTH1 (r->histos->hN[i]);
    fn[i] = new TF1(Form("noise%d",i),"gaus", -0.1, 0.1);
    r->histos->hN[i]->Fit(Form("noise%d",i),"RQ");
    r->histos->hN[i]->Draw();
    tl.Draw();
//     tex.DrawLatexNDC(0.72,0.55, Form("#splitline{Noise = %.1lf mV}{Thr = %.1lf mV}", 1000*r->histos->hN[i]->GetRMS(), thr[i]*1000)); 
    c1->Print(Form("%s/noise_ch%d.pdf",plotsDir.Data(), i+1));
    c1->Print(Form("%s/noise_ch%d.png",plotsDir.Data(), i+1));
  }
  
  // Pedestal subtraction
  
  for (int i=0; i< r->nch; i++){
    TH1I *hSub = (TH1I *)(r->histos->hC[i]->Clone());
    int upto_bin = hSub->FindBin(3*thr[i]);
    for (int ibin=1; ibin<=upto_bin; ibin++){
      double x = hSub->GetBinCenter(ibin);
      int oldy = hSub->GetBinContent(ibin);
      hSub->SetBinContent(ibin, (int)(oldy - fgaus[i].Eval(x)));
    }
//     r->histos->hC[i]->GetYaxis()->SetRangeUser(0, 3*r->histos->hC[i]->GetBinContent( r->histos->hC[i]->FindBin(3*thr[i])));
    hSub->Draw();
    TLine tl(thr[i], hSub->GetMinimum(), thr[i], hSub->GetMaximum());
    tl.SetLineWidth(2); tl.SetLineColor(4);
    tl.Draw();
    tex.DrawLatexNDC(0.43,0.85, Form("#splitline{Pedestal (#sigma) = %.1lf mV}{Thr (%.0lf#sigma) = %.1lf mV}", 1000*fgaus[i]->GetParameter(2), sigma_thr, thr[i]*1000)); 
    if (flangau[i]->GetNDF() > 0){
      tex.DrawLatexNDC(0.7,0.35, Form("#splitline{#splitline{#splitline{Langau fit:}{MPV = %.0lf mV}}{#sigma_{gaus} = %.1lf mV}}{#Chi^2/NDF = %.1lf}", 1000*flangau[i]->GetParameter(1), 1000*flangau[i]->GetParameter(3), flangau[i]->GetChisquare()/flangau[i]->GetNDF())); 
    }

    c1->Print(Form("%s/spectrum_ch%d_0sub.pdf",plotsDir.Data(), i+1));
    c1->Print(Form("%s/spectrum_ch%d_0sub.png",plotsDir.Data(), i+1));
  }
  
  // pedestal suppresion
  
  for (int i=0; i< r->nch; i++)
  {
//     gPad->SetLogy(0);
    r->histos->hC[i]->GetYaxis()->SetRangeUser(0, 3*r->histos->hC[i]->GetBinContent( r->histos->hC[i]->FindBin(3*thr[i])));
    r->histos->hC[i]->Draw();
    TLine tl(thr[i], 0, thr[i], r->histos->hC[i]->GetMaximum());
    tl.SetLineWidth(2); tl.SetLineColor(4);
    tl.Draw();
    tex.DrawLatexNDC(0.43,0.85, Form("#splitline{Pedestal (#sigma) = %.1lf mV}{Thr (%.0lf#sigma) = %.1lf mV}", 1000*fgaus[i]->GetParameter(2), sigma_thr, thr[i]*1000)); 
    if (flangau[i]->GetNDF() > 0){
      tex.DrawLatexNDC(0.7,0.35, Form("#splitline{#splitline{#splitline{Langau fit:}{MPV = %.0lf mV}}{#sigma_{gaus} = %.1lf mV}}{#Chi^2/NDF = %.1lf}", 1000*flangau[i]->GetParameter(1), 1000*flangau[i]->GetParameter(3), flangau[i]->GetChisquare()/flangau[i]->GetNDF())); 
    }
    c1->Print(Form("%s/spectrum_ch%d_0sup.pdf",plotsDir.Data(), i+1));
    c1->Print(Form("%s/spectrum_ch%d_0sup.png",plotsDir.Data(), i+1));
  }
  

  /*
  TLine tl= TLine(thr, 0, thr, hC->GetMaximum());
  tl->SetLineWidth(2); tl->SetLineColor(2);
  tl->Draw();
  Tl.DrawLatexNDC(0.4,0.65, Form("\# events above %.1lf mV: %.0lf (%.1lf \%)", thr*1000, hC->Integral(hC->FindBin(thr),-1), 100*hC->Integral(hC->FindBin(thr),-1)/r->N)); //   Tl.DrawLatex(.5, dy,   "x^{2y}");
  cout << "Number of events above thr. " << thr*1000 << " mV: " << hC->Integral(hC->FindBin(thr),-1) << endl;
 
  */
  cout << r->rel_offset << endl;
     

//   TCanvas *c3 = new TCanvas("c3", "c3", 0,0,800,900);
  TCanvas *c3 = new TCanvas("c3", "c3", 1);
  gPad->SetMargin(0.15,0.15,0.12,0.08);
//   DrawTH2(r->heff_tracks[0], "", 1.1, 1.1, 1);
  DrawTH2(r->histos->heff_tracks[0],"", 1.05, 1.2, 1.05);
  r->histos->heff_tracks[0]->Draw("COLZ");
  c3->Print(Form("%s/tracks_ch%d.pdf",plotsDir.Data(), 1));
  c3->Print(Form("%s/tracks_ch%d.png",plotsDir.Data(), 1));
  for (int ich=0; ich<r->nch; ich++)
  {
    for (int ix=1; ix<= r->histos->heff_efficiency[ich]->GetNbinsX(); ix++)
      for (int iy=1; iy<= r->histos->heff_efficiency[ich]->GetNbinsY(); iy++)
      {
//         int ntr=TMath::Max(r->heff_tracks[ich]->GetBinContent(ix,iy),1.0);
        int ntr=r->histos->heff_tracks[ich]->GetBinContent(ix,iy);
        if (ntr==0) continue;
        float eff = 1.0*r->histos->heff_hits[ich]->GetBinContent(ix,iy) / ntr;
        if (eff==0) eff=0.001;   // distinguish non efficient bins from empty ones
        r->histos->heff_efficiency[ich]->SetBinContent(ix,iy,eff);

        double chrg = r->histos->heff_charge[ich]->GetBinContent(ix,iy) / ntr;
        r->histos->heff_charge[ich]->SetBinContent(ix,iy, chrg); 
        r->histos->hpulse_height[ich]->Fill(chrg);
      }
  
    DrawTH2(r->histos->heff_efficiency[ich],"", 1.05, 1.2, 1.05);
    r->histos->heff_efficiency[ich]->GetZaxis()->SetRangeUser(0,1);
    r->histos->heff_efficiency[ich]->Draw("COLZ");
    c3->Print(Form("%s/eff_ch%d.pdf",plotsDir.Data(), ich+1));
    c3->Print(Form("%s/eff_ch%d.png",plotsDir.Data(), ich+1));
    
    DrawTH2(r->histos->heff_charge[ich],"", 1.05, 1.2, 1.05);
    r->histos->heff_charge[ich]->GetZaxis()->SetRangeUser(TMath::Min(0., r->histos->heff_charge[ich]->GetMinimum()), r->histos->heff_charge[ich]->GetMaximum());
    r->histos->heff_charge[ich]->Draw("COLZ");
    c3->Print(Form("%s/ph_ch%d.pdf",plotsDir.Data(), ich+1));
    c3->Print(Form("%s/ph_ch%d.png",plotsDir.Data(), ich+1));
    
    DrawTH1(r->histos->hpulse_height[ich],"", 1.05, 1.2);
    r->histos->hpulse_height[ich]->GetYaxis()->SetRangeUser(0, 2*r->histos->hpulse_height[ich]->GetBinContent(r->histos->hpulse_height[ich]->FindBin(50)));
    r->histos->hpulse_height[ich]->Draw("COLZ");
    c3->Print(Form("%s/phdistr_ch%d.pdf",plotsDir.Data(), ich+1));
    c3->Print(Form("%s/phdistr_ch%d.png",plotsDir.Data(), ich+1));
  }
  
  TGraph* g=r->avwf.Draw();
  g->Draw("APL");
  g->GetXaxis()->SetRangeUser(50,950);
  gPad->Update();
  gPad->Print(Form("%s/avgwf.png",plotsDir.Data()));
  gPad->Print(Form("%s/avgwf.pdf",plotsDir.Data()));
  
  // Draw efficiency with fiducial region marked
  for (int ich=0; ich<r->nch; ich++){
    r->histos->heff_efficiency[ich]->Draw("COLZ");
    TBox fid(fiducial_region[ich].x1, fiducial_region[ich].y1, fiducial_region[ich].x2, fiducial_region[ich].y2);
    fid.SetFillStyle(0);
    fid.SetLineWidth(2);
    fid.Draw();
    
    gPad->Print(Form("%s/eff_fid_ch%d.pdf",plotsDir.Data(), ich+1));
    gPad->Print(Form("%s/eff_fid_ch%d.png",plotsDir.Data(), ich+1));
  }
  
  TCanvas *c2 = new TCanvas("c2", "c2", 1);
  c2->SetLogy();
  TH1I* heff_distribution[NCH_DRS];
  for (int i=0; i< r->nch; i++)
  {
//     if (heff_distribution[i]) delete heff_distribution[i];
    heff_distribution[i] = new TH1I(Form("Eff_distr_ch%d",i+1), Form("Efficiency distribution ch. %d ; efficiency ; N_{bins}", i+1), 101, -0.005 ,1.005);
    for (int ix=1; ix<=r->histos->heff_efficiency[i]->GetNbinsX(); ix++)
      for (int iy=1; iy<=r->histos->heff_efficiency[i]->GetNbinsY(); iy++)
      {
        float xc = r->histos->heff_efficiency[i]->GetXaxis()->GetBinCenter( ix );
        float yc = r->histos->heff_efficiency[i]->GetYaxis()->GetBinCenter( iy );
        if (! fiducial_region[i].Contains(xc, yc)) continue;
        if (r->histos->heff_tracks[i]->GetBinContent(ix,iy)) heff_distribution[i]->Fill( r->histos->heff_efficiency[i]->GetBinContent(ix,iy));  
      }
    DrawTH1(heff_distribution[i]);
    heff_distribution[i]->GetYaxis()->SetRangeUser(0.4, 8*heff_distribution[i]->GetMaximum());
    heff_distribution[i]->Draw();
    tex.DrawLatexNDC(0.2,0.65, Form("#splitline{#splitline{Fiducial region ch. %d:}{x: [%.0lf, %.0lf]}}{y: [%.0lf, %.0lf]}", i+1, fiducial_region[i].x1, fiducial_region[i].x2, fiducial_region[i].y1, fiducial_region[i].y2)); 
    c2->Print(Form("%s/eff_distr_ch%d.pdf",plotsDir.Data(), i+1));
    c2->Print(Form("%s/eff_distr_ch%d.png",plotsDir.Data(), i+1));
  }
  //*/
  
  r->IOHand->f_res->Write();
}