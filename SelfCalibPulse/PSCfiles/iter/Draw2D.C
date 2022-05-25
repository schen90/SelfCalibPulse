#include "TStyle.h"
#include "TList.h"
#include "TRandom.h"
#include "TH2.h"
#include "TProfile2D.h"

void Draw2D(){
  // get hist 2
  TFile *f = new TFile("DiffTree.root");
  TProfile2D *p[6];
  for(int ix=0; ix<3; ix++){
    p[ix] = (TProfile2D *)f->Get(Form("pxy%d",ix));
    p[ix+3] = (TProfile2D *)f->Get(Form("pxz%d",ix));
    p[ix]->GetZaxis()->SetRangeUser(0,7);
    p[ix+3]->GetZaxis()->SetRangeUser(0,7);
  }
  
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

  
}
