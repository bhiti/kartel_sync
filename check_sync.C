#define NCH_DRS 8
#define DRS_N_POINTS 1024
#define NPTS_ANAWF 150

#include "TStyle.h"
  
#include "src/waveform.hpp"
#include "src/event.hpp"
#include "src/histograms.hpp"
#include "src/DUTCorrelation.hxx"
#include "src/IOHandler.hpp"
#include "src/run.hpp"

#include "src/waveform.cpp"
#include "src/event.cpp"
#include "src/histograms.cpp"
#include "src/DUTCorrelation.cxx"
#include "src/IOHandler.cpp"
#include "src/run.cpp"

// void check_sync(const char* fname_tel="../../2018/oct2018/data/conv/anchor1154.root", const char* fname_anchor="../../2018/oct2018/data/conv/run1154.root", const char* fname_synced="", int Nevents = 123456789)
void check_sync(const char* fname_tel="../../../2018/oct2018/data/conv/tel_final1154.root", const char* fname_anchor="../../../2018/oct2018/data/conv/ref1154.root", const char* fname_synced="", int Nevents = 123456789)
{
  gErrorIgnoreLevel = kWarning;
  gStyle->SetPalette(55);
    
//  char fname_tel[64], fname_anchor[64], fname_synced[64];
//   const char* run_number="087";
//  const char* run_number= Form("%03d", runNumber);
//  cout << "Run number: " << run_number << endl;
  
//   int Nevents = 200000;     // EDIT HERE: number of events to process
//  Nevents = 3333333;     // EDIT HERE: number of events to process
  int verbose=0;

  int sync_mode=0;    // 0... sync DRS+tel, 1...sync FEI4+tel

//   if (strlen(fname_tel) == 0 || strlen(fname_anchor) == 0 || strlen(fname_synced) == 0){
//     cout << "Syncing: paths not specified. Exiting. " << endl;
//     return;
//   }
//     sprintf(fname_tel, "../data/telescope/ref%s.root", run_number);        // Telescope file for synchronization with DUT
//     sprintf(fname_anchor, "../data/telescope/anchor%s.root", run_number);   // Anchor module data (FEI4)
//     sprintf(fname_synced, "../data/proteusFiles/synced%s.root", run_number);                      // synchronized telescope and FEI4 events, merged into one file
//   }

  run *r = new run(fname_anchor, fname_tel, fname_anchor, fname_synced);
//   r->AddTelescopePlane(0);
//   r->AddTelescopePlane(1);
//   r->AddTelescopePlane(2);
//   r->AddTelescopePlane(3);
//   r->AddTelescopePlane(4);
//   r->AddTelescopePlane(5);
  
//   if (r->CreateSyncedFile() < 0) return;
  
//   float mean_ratio = 0.5;
//   float tolerance_factor = 1.02;
  float mean_ratio = 1;      // tel+fei4=0.5, tel+drs=0.000125, drs+fei4=0.00025
//   float mean_ratio = 0.000125;      // tel+fei4=0.5, tel+drs=0.000125, drs+fei4=0.00025
  float tolerance_factor = 1.5;
  r->SetRatio(mean_ratio, tolerance_factor);          // ratio of time stamp deltas for syncing (target_ratio, tolerance_factor)
  
  int istart = 0;
//   int iend = TMath::Min(r->N, Nevents)-100;
  int iend = TMath::Min(r->N, Nevents);
  int i;
//   int istart = 14330;
//   r->rel_offset = 214;
  
  int passed=0; 
  int big=0;
  
  cout << "Run has " << r->N << " events" << endl;
  cout << "Analyzing " << iend << " events" << endl;
  cout << "Event \trel. \tpassed" << endl;
  for (i=istart; i<iend; i++) 
  {
    if (r->N < i) {
      cout << r->N << " " << i << endl;
      break;    
    }
    if (i%10000 == 0){
      cout << i << "\t" << r->rel_offset << "\t" << passed << endl;
//       cout << r->IOHand->outputEvents->GetEntries() << " " << i << endl;
    }
    
    // synchronization
    if (r->ProcessTimestamp(i,sync_mode,verbose) != 0){     // Find out if offseting is required
      // If yes, then do offseting
      int iskipped = r->ProcessOffsets(i);                   
      
      i += iskipped;
      if (iskipped > 2){
        i--;
        big++;
        continue;                           // if offseting for more than 3 events is required, then check again
      }
    }
//     r->FillSyncedFile(i);
    passed++;
  }  
  
  cout << i << "\t" << r->rel_offset << "\t" << passed << endl;
//   cout << "Telescope events: " << r->IOHand->t_tel_ts->GetEntries() << endl;
//   cout << "FEI4 events: " << r->IOHand->t_anchor_ts->GetEntries() << " (diff. " << r->IOHand->t_anchor_ts->GetEntries() - r->IOHand->t_tel_ts->GetEntries() << ")" << endl;
// //   cout << "Events in synced file: " << r->IOHand->outputEvents->GetEntries() << endl;
// //   r->WriteAll();
  
//   r->hRatio->Draw();
// //   cout << r->hRatio->Integral(6,20) << endl;
//   
// //   r->IOHand->t_anchor_ts->GetEntry(0);
// //   cout << r->Event->ts_anchor << endl;
// //   r->IOHand->t_anchor_ts->GetEntry(r->IOHand->t_anchor_ts->GetEntries()-2);
// //   cout << r->Event->ts_anchor << endl;
// //   r->IOHand->t_anchor_ts->GetEntry(r->IOHand->t_anchor_ts->GetEntries());
// //   cout << r->Event->ts_anchor << endl;
//   
// //   cout << big << " run: " << r->bigOffset << endl;
}

