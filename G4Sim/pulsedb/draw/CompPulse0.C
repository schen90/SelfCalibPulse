#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

const int BSIZE = 56;
const int NSEGS = 36;
const int NCHAN = NSEGS+1;

TGraph *gr;
TGraph *gr2;

void CompPulse0(int ientry = 4385){
  //---------------------------
  Int_t simseg = 0;
  Double_t simpos[3] = {0,0,0};
  //---------------------------

  TCanvas *c = new TCanvas("c","c",1265,400);
  c->SetMargin(0.06,0.01,0.12,0.01);

  // ori simulation
  TFile *f = new TFile("before/LibTrap_A001.root");
  TTree *tree = (TTree *)f->Get("tree");
  Int_t seg;
  Float_t pos[3];
  //Double_t core[BSIZE];
  Float_t spulse[BSIZE*NCHAN];
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("pos",pos);
  //tree->SetBranchAddress("core",core);
  tree->SetBranchAddress("spulse",spulse);
  int nentries = tree->GetEntriesFast();

  // with XT and preAmp
  //TFile *f2 = new TFile("Xtalk/LibTrap_A001.root");
  TFile *f2 = new TFile("LibTrap_A001.root");
  TTree *tree2 = (TTree *)f2->Get("tree");
  Int_t seg2;
  Float_t pos2[3];
  //Double_t core2[BSIZE];
  Float_t spulse2[BSIZE*NCHAN];
  tree2->SetBranchAddress("seg",&seg2);
  tree2->SetBranchAddress("pos",pos2);
  //tree2->SetBranchAddress("core",core2);
  tree2->SetBranchAddress("spulse",spulse2);
  
 
  tree->GetEntry(ientry);
  tree2->GetEntry(ientry);
  
  Double_t x[BSIZE*NCHAN];
  for(int iseg=0; iseg<NCHAN; iseg++)
    for(int i=0; i<BSIZE; i++)
      x[iseg*BSIZE+i] = iseg+i*1./BSIZE-0.5;

  double tmpspulse[BSIZE*NCHAN];
  for(int i=0; i<BSIZE*NCHAN; i++) tmpspulse[i] = spulse[i];
  //for(int i=0; i<BSIZE; i++) tmpspulse[i+BSIZE*NSEGS] = core[i];
  
  double tmpspulse2[BSIZE*NCHAN];
  for(int i=0; i<BSIZE*NCHAN; i++) tmpspulse2[i] = spulse2[i];
  //for(int i=0; i<BSIZE; i++) tmpspulse2[i+BSIZE*NSEGS] = core2[i];

  gr = new TGraph(BSIZE*NCHAN,x,tmpspulse);
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
  //gr->GetXaxis()->SetRangeUser(19.5,25.5);
  gr->GetYaxis()->SetRangeUser(-0.15,1.1);

  gr2 = new TGraph(BSIZE*NCHAN,x,tmpspulse2);
  gr2->SetLineColor(2);

  cout<<Form("sim seg: %d  pos: %f  %f  %f",seg,pos[0],pos[1],pos[2])<<endl;
  cout<<Form("Lib seg: %d  pos: %f  %f  %f",seg2,pos2[0],pos2[1],pos2[2])<<endl;
  
  gr->Draw("APL");
  gr2->Draw("Lsame");
  
  f->Close();
  return;
}
