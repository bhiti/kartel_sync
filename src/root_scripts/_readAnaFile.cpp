{
#define NCH 2
#define NPTS_WF 150   // number of points in the waveforms = 150
  TString fname = "data/ana/ana087.root";

  // ch1 = 4x4, ch2 = 2x2
  Float_t wf[NCH][NPTS_WF];
  Float_t charge[NCH];
  
  Double_t xt;      // track position in the telescope
  Double_t yt;
  Double_t slopeX;
  Double_t slopeY;
  Double_t chi2;
  Double_t xa;      // cluster center in the anchor module
  Double_t ya;
  
  TFile *f = new TFile(fname.Data());
  TTree *tch[NCH];
  TTree *ttracks;
  
  for (int i=0; i<NCH; i++){
    tch[i] = (TTree*) f->Get(Form("ch%d", i+1));
    tch[i]->SetBranchAddress("waveform", wf[i]);
    tch[i]->SetBranchAddress("charge15ns", &charge[i]);
  }
  
  ttracks = (TTree*) f->Get("Tracks");
  ttracks->SetBranchAddress("xt", &xt);
  ttracks->SetBranchAddress("yt", &yt);
  ttracks->SetBranchAddress("chi2", &chi2);
  ttracks->SetBranchAddress("slopeX", &slopeX);
  ttracks->SetBranchAddress("slopeY", &slopeY);
  ttracks->SetBranchAddress("xa", &xa);
  ttracks->SetBranchAddress("ya", &ya);
  
  int Nevents = ttracks->GetEntries();
  
  for (int i=0; i<Nevents; i++){
   // do your magic here 
    for (int ch=0; ch<NCH; ch++) tch[ch]->GetEntry(i);
    ttracks->GetEntry(i);
    
    if (i%10000==0) cout << i << " " << xa << " " << charge[0] << endl;
  }
  
}
  
  