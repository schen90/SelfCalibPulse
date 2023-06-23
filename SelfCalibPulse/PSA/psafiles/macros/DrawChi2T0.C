#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLegend.h"

using namespace std;

void DrawChi2T0(){
  const int NDet = 3;
  int DetId[NDet] = {0,1,2};
  int color[NDet] = {1, 2, 3};
  const int NT0 = 7;
  int T0[NT0] = { 9, 10, 11, 12, 13, 14, 15};

  float x[NT0];
  float y[NDet][NT0];
  float ymin = 30000, ymax = 0;

  for(int ip=0; ip<NT0; ip++){
    x[ip] = T0[ip];
    TFile *f = new TFile(Form("T0/Run0008_0000_T0_%02d.root",T0[ip]));
    TTree *tree = (TTree *)f->Get("tree");
    for(int id=0; id<NDet; id++){
      TH1D *htmp = new TH1D("htmp","",1000,0,30000);
      tree->Draw("chi2>>htmp",Form("CrystalId==%d",DetId[id]),"goff");
      //tree->Draw("chi2/CoreE[0]>>htmp",Form("CrystalId==%d",DetId[id]),"goff");
      y[id][ip] = htmp->GetMean();
      if(y[id][ip]<ymin) ymin=y[id][ip];
      if(y[id][ip]>ymax) ymax=y[id][ip];
    }
    f->Close();
  }

  TGraph *gr[NDet];
  for(int idet=0; idet<NDet; idet++){
    gr[idet] = new TGraph(NT0, x, y[idet]);
    gr[idet]->SetLineColor(color[idet]);
    gr[idet]->SetNameTitle(Form("gr%d",DetId[idet]),"Chi2 vs T0");
  }

  TLegend *label = new TLegend(0.15,0.6,0.35,0.89);
  gr[0]->GetYaxis()->SetRangeUser(0.9*ymin, 1.1*ymax);
  for(int idet=0; idet<NDet; idet++){
    label->AddEntry(gr[idet],Form("Det %02d",DetId[idet]),"l");
    if(idet==0) gr[idet]->Draw("APL");
    else        gr[idet]->Draw("lsame");
  }
  label->Draw("same");

}
