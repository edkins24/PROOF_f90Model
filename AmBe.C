/*

 Usage:
 $ root -b -q Argon.C(+)

 To run in compiled mode, use the "+".

 2016-04-06:
 - modernized to SLAD v2.3
 - improved example with more flexibility
 
  
 */

#include "rootstart.h"

#include <iostream>

#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TSelector.h"
#include "TProof.h"

using namespace std;

Bool_t load_chan = false;
Bool_t load_allpulses = true;
Bool_t load_masas_xy = true;
Bool_t load_xyl_xy = true;
Bool_t load_aww_xy = true;
Bool_t load_veto = true;

//------------------------------------------------------------------------------
// Main event loop is contained here.
void event_loop(TChain* events)
{
  TProof::Open("");
  events->SetProof();
  events->Process("ProofEventDataSelector.C+");
}

//------------------------------------------------------------------------------
// open SLAD files
TChain* load_files(TString mainFileName)
{
  TFile* main = new TFile(mainFileName, "READ");
  TChain* events = new TChain("events");
  events->Add(mainFileName);
  TChain* logbook = new TChain("logbook");
  logbook->Add(mainFileName);
  events->AddFriend(logbook);

  // Construct the name of the friend tree files and open them
  // Extract the trees and friend them.

  if (load_chan) {
    TString chanFileName = mainFileName;
    chanFileName.Remove(chanFileName.Length()-5);
    chanFileName += "_chan.root";
    TFile* chanFile = TFile::Open(chanFileName);
    TChain* s2_chan = new TChain("s2_chan");
    s2_chan->Add(chanFileName);
    events->AddFriend(s2_chan);
  }
  
  if (load_allpulses) {
    TString pulseFileName = mainFileName;
    pulseFileName.Remove(pulseFileName.Length()-5);
    pulseFileName += "_allpulses.root";
    TFile* pulseFile = TFile::Open(pulseFileName);
    TChain* pulse_info = new TChain("pulse_info");
    pulse_info->Add(pulseFileName);
    events->AddFriend(pulse_info);
  }
  
  if (load_masas_xy) {
    TString masasxyFileName = mainFileName;
    masasxyFileName.Remove(masasxyFileName.Length()-5);
    masasxyFileName+="_masas_xy.root";
    TFile* masasxyFile = TFile::Open(masasxyFileName);
    TChain* masas_xy = new TChain("masas_xy");
    masas_xy->Add(masasxyFileName);
    TChain* masas_xy_allpulses = new TChain("allpulses_xy");
    masas_xy_allpulses->Add(masasxyFileName);
    events->AddFriend(masas_xy);
    events->AddFriend(masas_xy_allpulses);
  }
  
  if (load_xyl_xy) {
    TString xylxyFileName = mainFileName;
    xylxyFileName.Remove(xylxyFileName.Length()-5);
    xylxyFileName+="_xylocator_xy.root";
    TFile* xylxyFile = TFile::Open(xylxyFileName);
    TChain* xyl_xy = new TChain("xylocator_xy");
    xyl_xy->Add(xylxyFileName);
    TChain* xyl_xy_allpulses = new TChain("allpulses_xyl_xy");
    xyl_xy_allpulses->Add(xylxyFileName);
    events->AddFriend(xyl_xy);
    events->AddFriend(xyl_xy_allpulses);
  }

  if (load_aww_xy) {
    TString awwxyFileName = mainFileName;
    awwxyFileName.Remove(awwxyFileName.Length()-5);
    awwxyFileName+="_aww_xy.root";
    TFile* awwxyFile = TFile::Open(awwxyFileName);
    TChain* aww_xy = new TChain("aww_xy");
    aww_xy->Add(awwxyFileName);
    TChain* aww_xy_allpulses = new TChain("allpulses_aww_xy");
    aww_xy_allpulses->Add(awwxyFileName);
    events->AddFriend(aww_xy);
    events->AddFriend(aww_xy_allpulses);
  }
  
  if (load_veto) {
    TString vetoFileName = mainFileName;
    vetoFileName.Remove(vetoFileName.Length()-5);
    vetoFileName+="_veto_cluster.root";
    TFile* vetoFile = TFile::Open(vetoFileName);
    TChain* veto = new TChain("veto");
    veto->Add(vetoFileName);
    events->AddFriend(veto);
  }

  return events;
}//load_files()

//------------------------------------------------------------------------------
// Main method. Load DST files and invoke event_loop().
void AmBe()
{
   // Prevent canvases from being drawn.
  gROOT->SetBatch(kTRUE);
  
  SetMyStyle();

  TStopwatch* clock = new TStopwatch();
  clock->Start();

  // The main SLAD file containing the data we want.
  TString mainFileName = "/scratch/darkside/slad/AmBe_160nps_SLAD_v2.3.2.root";
  
  TChain* events = load_files(mainFileName);

  event_loop(events);
  
  cout << "Done! " << clock->RealTime() << " s." << endl;
}
