// #include "DUTCorrelation.hxx"

DUTCorrelation::DUTCorrelation(){};
DUTCorrelation::DUTCorrelation(window ROI, int nbinsx, int nbinsy, int nch):_ROI(ROI), _nch(nch)
{
  _ROI.x1=ROI.x1; _ROI.x2=ROI.x2; _ROI.y1=ROI.y1; _ROI.y2=ROI.y2;
//   cout << _ROI.x1 << " " << _ROI.x2 << endl;;  
//   cout << ROI.x1 << " " << ROI.x2 << endl;;

  for (int i=0; i<_nch; i++)
    for (int j=0; j<_nch; j++)
      _hHits[i][j] = new TH2S(Form("CorrHits%d_%d", i,j), Form("CorrHits%d_%d", i,j), nbinsx, _ROI.x1, _ROI.x2, nbinsy, _ROI.y1, _ROI.y2);
};
DUTCorrelation::~DUTCorrelation()
{
//   for (int i=0; i<NCH_DRS; i++)
//     for (int j=0; j<NCH_DRS; j++)
//       if (_hHits[i][j]) delete _hHits[i][j];
//   if (_hWeights) delete _hWeights;
};
  
void DUTCorrelation::ProcessEvent(event &Event, float *thr)
{
  int itr=Event.SearchTrack(_ROI); // index of the first track going through the roi
//   cout << _ROI.x1 << " " << _ROI.x2 << endl;;
  if(itr >= 0){
    for(int i=0; i<_nch; i++){
      if (Event._charge[i] > thr[i])
      {
        _hHits[i][i]->Fill(Event.x[itr], Event.y[itr]);                           // increment reference counter
        for (int j=0; j<_nch; j++){        
          // count hits in one DUT given a hit is present in the other DUT
          // loop over all possible combinations
          // counts are saved in a 2d historgram to allow for applying beam intensity weights
          if (i==j) continue;   // skip autocorrelation
          if (Event._charge[j] > thr[j]) _hHits[i][j]->Fill(Event.x[itr], Event.y[itr]);   // increment probe counter
        }
      }
    }
  }
  return;
};

TH2S *DUTCorrelation::AddWeight(TH2S* h)
{
  _hWeights = h;
  return _hWeights;
};


float DUTCorrelation::GetCorrelation(int ch_probe, int ch_ref)
{
  // Return the fraction of hits over threshold in channel [ch_probe] with respect to the total number of hits over threshold in channel [ch_ref]
  // Each bin is weighted by the inverse of the number of telescope tracks through the bin to account for beam non-uniformity
  
  float sum_probe=0, sum_ref=0;
  for (int i=1; i<= _hHits[ch_ref][ch_ref]->GetNbinsX(); i++)
    for (int j=1; j<= _hHits[ch_ref][ch_ref]->GetNbinsX(); j++){
      if (_hWeights->GetBinContent(i,j) > 0){
        sum_probe += 1.0*_hHits[ch_ref][ch_probe]->GetBinContent(i,j) / _hWeights->GetBinContent(i,j);
        sum_ref += 1.0*_hHits[ch_ref][ch_ref]->GetBinContent(i,j) / _hWeights->GetBinContent(i,j);
//         sum_probe += 1.0*_hHits[ch_ref][ch_probe]->GetBinContent(i,j);
//         sum_ref += 1.0*_hHits[ch_ref][ch_ref]->GetBinContent(i,j);
      }
    }
    
  if (sum_ref==0) return -2;  // Error: no counts!
  return sum_probe/sum_ref;
};  