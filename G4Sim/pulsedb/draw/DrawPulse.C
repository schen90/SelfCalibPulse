#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

TGraph *gr;

void DrawPulse(){
  //---------------------------
  Int_t simseg = 0;
  Double_t simpos[3] = {0,0,0};
  //---------------------------

  TCanvas *c = new TCanvas("c","c",1265,400);
  c->SetMargin(0.06,0.01,0.12,0.01);
  TFile *f = new TFile("pulseA.root");
  TTree *tree = (TTree *)f->Get("tree");
  Int_t seg;
  Double_t pos[3];
  Double_t core[121];
  Double_t spulse[4356];
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("pos",&pos);
  tree->SetBranchAddress("core",core);
  tree->SetBranchAddress("spulse",spulse);
  int nentries = tree->GetEntriesFast();

  tree->GetEntry(20300);
  
  Double_t x[4477];
  for(int iseg=0; iseg<37; iseg++)
    for(int i=0; i<121; i++)
      x[iseg*121+i] = iseg+i*1./121-0.5;

  double tmpspulse[4477];
  for(int i=0; i<4356; i++) tmpspulse[i] = spulse[i];
  for(int i=0; i<121; i++) tmpspulse[i+4356] = core[i];
  
  gr = new TGraph(4477,x,tmpspulse);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetTitle("Segment Number");
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitleOffset(0.5);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetTitle("Normalised Charge");
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetRangeUser(-0.5,36.5);
  gr->GetYaxis()->SetRangeUser(-0.1,1.1);
  gr->Draw("APL");
  
  f->Close();
  return;
}
