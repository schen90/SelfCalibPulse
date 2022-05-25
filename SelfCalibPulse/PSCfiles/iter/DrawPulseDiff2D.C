#include "TStyle.h"
#include "TList.h"
#include "TRandom.h"
#include "TH2.h"
#include "TProfile2D.h"

void DrawPulseDiff2D(){
  // get hist 2
  //gStyle->SetOptStat(0);
  TFile *f = new TFile("PulseDiffTree.root");
  TProfile2D *p[6];
  for(int isli=0; isli<6; isli++){
    p[isli] = (TProfile2D *)f->Get(Form("pxy%d",isli));
    p[isli]->GetZaxis()->SetRangeUser(0,20);
    p[isli]->SetTitle(Form("Y : X {slice%d}",isli+1));
  }
  TProfile2D *pxz = (TProfile2D *)f->Get("pxz");
  pxz->SetTitle("Z : X");
  pxz->GetZaxis()->SetRangeUser(0,20);


  TCanvas *c = new TCanvas("c","abc",900,600);
  c->Divide(3,2);
  c->cd(1);
  p[0]->Draw("colz");
  c->cd(2);
  p[1]->Draw("colz");
  c->cd(3);
  p[2]->Draw("colz");
  c->cd(4);
  p[3]->Draw("colz");
  c->cd(5);
  p[4]->Draw("colz");
  c->cd(6);
  p[5]->Draw("colz");


  /*
  TCanvas *c = new TCanvas("c","abc",300,300);
  pxz->Draw("colz");
  */
    
}
