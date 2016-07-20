/*
 
 This class is derived from the ROOT class TSelector.
 $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.
 
 The following methods are defined in this file:
    SlaveBegin():   called after Begin(), when on PROOF called only on the slave
                    servers, a convenient place to create your histograms.
    Process():      called for each event, in this function you decide what to read
                    and fill your histograms.
    Terminate():    called at the end of the loop on the tree,
                    a convenient place to draw/fit your histograms.
 */

#define ProofEventDataSelector_cxx
#include "ProofEventDataSelector.h"

#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TMath.h"

using namespace std;

void ProofEventDataSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  /* defining histos */
  h_lifetime = new TH1D("h_lifetime", "Lifetime counter; ; Livetime [s]", nCuts, -0.5, nCuts-0.5);
  for (Int_t i=1; i<=nCuts; ++i) h_lifetime->GetXaxis()->SetBinLabel(i, EventCut_label[i-1]);
  
  h_s1 = new TH1D("h_s1", "; S1 [PE]; a.u.", 500, 0., 6000.);
  h_s1->Sumw2();
  h_s1->SetMarkerStyle(kFullCircle);
  
  h_s2 = new TH1D("h_s2", "; S2 [PE]; a.u.", 500, 0., 60000.);
  h_s2->Sumw2();
  h_s2->SetMarkerStyle(kFullCircle);
  
  h_f90vss1 = new TH2D("h_f90vss1", "; S1 [PE]; f90", 1000, 0., 3000., 1000, 0., 1.);
  h_f90vss1->Sumw2();
  
  h_log10s2overs1vsf90 = new TH2D("h_log10s2overs1vsf90", "; f90; Log_{10}(S2/S1)", 500, 0., 1., 500, -1., 4.);
  h_log10s2overs1vsf90->Sumw2();
  
  h_ene = new TH1D("h_ene", "; E [keV]; a.u.", 500, 0., 600.);
  h_ene->Sumw2();
  h_ene->SetMarkerStyle(kFullCircle);
  
  h_lys1vsene = new TH2D("h_lys1vsene", "; E [keV]; LY_{S1} [PE/keV]", 500, 0., 600., 500, 0., 10.);
  h_lys1vsene->Sumw2();
  
  h_f90vsene = new TH2D("h_f90vsene", "; E [keV]; f90", 500, 0., 600., 500, 0., 1.);
  h_f90vsene->Sumw2();
  
  h_rvsene = new TH2D("h_rvsene", "; E [keV]; r", 500, 0., 600., 500, 0., 1.);
  h_rvsene->Sumw2();
  
  h_veto_cluster_charge_vec0 = new TH1D("h_veto_cluster_charge_vec0", "; PE; a.u.", 100, 0., 5000.);
  h_veto_cluster_charge_vec0->Sumw2();
  h_veto_cluster_charge_vec0->SetMarkerStyle(kFullCircle);
  
  // Add to output list (needed for PROOF)
  GetOutputList()->Add(h_lifetime);
  GetOutputList()->Add(h_s1);
  GetOutputList()->Add(h_s2);
  GetOutputList()->Add(h_f90vss1);
  GetOutputList()->Add(h_log10s2overs1vsf90);
  GetOutputList()->Add(h_ene);
  GetOutputList()->Add(h_lys1vsene);
  GetOutputList()->Add(h_f90vsene);
  GetOutputList()->Add(h_rvsene);
  GetOutputList()->Add(h_veto_cluster_charge_vec0);
  
  prepareS1MF();

}

Bool_t ProofEventDataSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree to be
  // processed. The entry argument specifies which entry in the currently
  // loaded tree is to be processed.
  // It can be passed to either EventSelector::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the TTree.
  //
  // This function should contain the "body" of the analysis: select relevant
  // tree entries, run algorithms on the tree entry and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  // *** 1. *** Tell the reader to load the data for this entry:
  events_reader.SetEntry(entry);
  
  // *** 2. *** Do the actual analysis
  // variables
  Boot_t isPrescaled = *tpc_digital_sum >= 250 && *tpc_digital_sum < 1200;
  Float_t s1 = *total_s1_corr;
  Float_t xycorr = (*masas_xycorr_factor != 1. ? *masas_xycorr_factor : (*xyl_xycorr_factor != 1. ? *xyl_xycorr_factor : (*aww_xycorr_factor != 1. ? *aww_xycorr_factor : -998.)));
  Float_t s2 = *total_s2_corr * xycorr;
  Float_t log10s2overs1 = (s2>0. ? TMath::Log10(s2/s1) : -999.);
  Float_t f90 = *total_f90;  
  
  // xy reconstruction
  Float_t x_coor = (*masas_x > -99. ? *masas_x : (*xyl_x > -99. ? *xyl_x : (*aww_x > -99. ? *aww_x : -998.)));
  Float_t y_coor = (*masas_y > -99. ? *masas_y : (*xyl_y > -99. ? *xyl_y : (*aww_y > -99. ? *aww_y : -998.)));
  Float_t radius = (x_coor>0. && y_coor>0.) ? TMath::Sqrt(x_coor*x_coor + y_coor*y_coor) : -998.;  
  
  // values got from dokeplots_UAr.root
  const Float_t g2overg1_factor = 178.296;
  const Float_t alpha_factor = 0.000673624;
  // REFERENCE: T. Doke, A. Hitachi, J. Kikuchi, K. Masuda, H. Okada, and E. Shibamura, Jpn. J. Appl. Phys. 41, 1538-1545 (2002).
  const Float_t W = 19.5E-3;
  const Float_t f = 0.21;
  const Float_t g2 = W / alpha_factor;
  const Float_t g1 = g2 / g2overg1_factor;

  Float_t ene = alpha_factor * (s2 + g2overg1_factor * s1);
  
  Float_t lys1 = s1 / ene;
  Float_t ng = s1 / g1;
  Float_t ne = s2 / g2;
  Float_t y = ne / (ne + ng);
  Float_t r = 1. - (f + 1.) * y;    

  // cuts
  Bool_t CXNChannels = *nchannels == 38;
  Bool_t CXBaseline = *baseline_not_found == false;
  Bool_t CXEventDt = (*lifetime + *inhibittime) >= 1.35E-3;
  Bool_t CXLongLivetime = *lifetime < 1.;
  
  Bool_t CXBasic = CXNChannels && CXBaseline && CXEventDt && CXLongLivetime;
  
  Bool_t CXTriggerTime = ((*run_id >= -999 && *run_id < 7344) && (*s1_start_time >= -0.25 && *s1_start_time <= -0.15)) ||
                         ((*run_id >= 7344 && *run_id < 7641) && (*s1_start_time >= -4.10 && *s1_start_time <= -4.00)) ||
                         ((*run_id >= 7641 && *run_id < 999999) && (*s1_start_time >= -6.10 && *s1_start_time <= -6.00));
  Bool_t CXSingleScatter = *npulses == 2 || (*npulses == 3 && *has_s3);
  Bool_t CXS1MF = *max_s1_frac_cut_exceeds99 == 0;
  Bool_t CXS2F90 = *total_s2_f90_fixed < 0.2;
  //Bool_t CXS2Size = s2 > 100. && (radius > 0. && radius < 20.);
  Bool_t CXS2Size = s2 > 30.
  
  Bool_t CXVetoPresent = *veto_present == true;
  Bool_t CXVetoEventMatch = *run_id == *veto_run_id && *event_id == *veto_event_id;
  
  Bool_t CXAnalysis = CXSingleScatter && CXS1MF && CXS2F90 && CXVetoPresent && CXVetoEventMatch;
  
  if (CXBasic && CXAnalysis) {
    Float_t w = isPrescaled ? 33. : 1.;
    
    h_s1->Fill(s1, w);
    h_s2->Fill(s2, w);
    h_f90vss1->Fill(s1, f90, w);
    h_log10s2overs1vsf90->Fill(f90, log10s2overs1, w);
    h_ene->Fill(ene, w);
    h_lys1vsene->Fill(ene, lys1, w);
    h_f90vsene->Fill(ene, f90, w);
    h_rvsene->Fill(ene, r, w);
    
    if (veto_cluster_charge_vec->size() > 0) h_veto_cluster_charge_vec0->Fill(veto_cluster_charge_vec->at(0));
  }
  
  return kTRUE;
}
  
void ProofEventDataSelector::Terminate()
{
  // The Terminate() function is the last function to be called during the
  // analysis of a tree with a selector. It always runs on the client, it can
  // be used to present the results graphically or save the results to file.
  
  // finalize lifetime
  const Double_t total_lf = h_lifetime->GetBinContent(NOCUTS+1);
  for (Int_t i=NCHANNEL+1; i<=nCuts; ++i) h_lifetime->SetBinContent(i, total_lf - h_lifetime->GetBinContent(i));
  Double_t total_lifetime = h_lifetime->GetBinContent(LONGWAIT+1);
  
  TFile* fout = new TFile("results.root", "RECREATE");
  
  fout->cd();
  h_lifetime->Write();
  h_s1->Write();
  h_s2->Write();
  h_f90vss1->Write();
  h_log10s2overs1vsf90->Write();
  h_ene->Write();
  h_lys1vsene->Write();
  h_f90vsene->Write();
  h_rvsene->Write();
  h_veto_cluster_charge_vec0->Write();
  fout->Close();
} 
