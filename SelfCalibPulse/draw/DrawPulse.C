#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"

TGraph *gr1;
TGraph *gr2;
TGraph *gr3;
TGraph *gr4;

void DrawPulse(){
  float xmin = -0.5, xmax = 36.5;
  
  TCanvas *c = new TCanvas("c","c",1000,600);
  //c->SetMargin(0.05,0.01,0.05,0.01);
  //TFile *f = new TFile("ComparePulse_chi2max.root");
  //TFile *f = new TFile("ComparePulse_chi2s.root");
  TFile *f = new TFile("ComparePulse_noise_chi2s.root");
  TTree *tree = (TTree *)f->Get("tree");
  int simseg;
  int nhits1, nhits2, nfired;
  vector<float> *hiteng1 = 0;
  vector<float> *hiteng2 = 0;
  float Energy1, Energy2;
  float PhiRZ1[3], PhiRZ2[3];
  Float_t pulse1[37][121], pulse2[37][121], chis[37][121];
  float chi2, dist, rdiffphi, diffphi, diffr, diffz;
  tree->SetBranchAddress("simseg",&simseg);
  tree->SetBranchAddress("nhits1",&nhits1);
  tree->SetBranchAddress("nhits2",&nhits2);
  tree->SetBranchAddress("Energy1",&Energy1);
  tree->SetBranchAddress("Energy2",&Energy2);
  tree->SetBranchAddress("chi2",&chi2);
  tree->SetBranchAddress("nfired",&nfired);
  tree->SetBranchAddress("dist",&dist);
  tree->SetBranchAddress("rdiffphi",&rdiffphi);
  tree->SetBranchAddress("diffphi",&diffphi);
  tree->SetBranchAddress("diffr",&diffr);
  tree->SetBranchAddress("diffz",&diffz);

  tree->SetBranchAddress("PhiRZ1",&PhiRZ1);
  tree->SetBranchAddress("PhiRZ2",&PhiRZ2);
  tree->SetBranchAddress("hiteng1",&hiteng1);
  tree->SetBranchAddress("hiteng2",&hiteng2);
  tree->SetBranchAddress("pulse1",pulse1);
  tree->SetBranchAddress("pulse2",pulse2);
  tree->SetBranchAddress("chis",chis);

  Float_t x[4477];
  Float_t y1[4477], y2[4477], y3[4477], y4[4477];
  for(int iseg=0; iseg<37; iseg++)
    for(int i=0; i<121; i++)
      x[iseg*121+i] = iseg+i*1./121-0.5;

  int nentries = tree->GetEntriesFast();
  int ientry;
  for(ientry = 0; ientry<nentries; ientry++){
    tree->GetEntry(ientry);
    if(simseg==10 && chi2>0 && Energy1<700 && Energy2<700 && Energy1>600 && Energy2>600)
      //if(chi2<2 && Energy1>1000 && Energy2>1000 && diffr<1 && diffz<1 && rdiffphi>5)
      //if(chi2>1 && Energy1>1000 && Energy2>1000 && dist<5)
    {
      cout<<"entry "<<ientry<<endl
	  <<" seg    = "<<simseg<<endl
	  <<" nhits1 = "<<nhits1<<" : ";
      for(int ii=0; ii<hiteng1->size(); ii++) cout<<Form("%.2f  ",hiteng1->at(ii));
      cout<<endl
	  <<" nhits2 = "<<nhits2<<" : ";
      for(int ii=0; ii<hiteng2->size(); ii++) cout<<Form("%.2f  ",hiteng2->at(ii));
      cout<<endl
	  <<" nfired = "<<nfired<<endl
	  <<" Eng1   = "<<Energy1<<endl
	  <<" Eng2   = "<<Energy2<<endl
	  <<" chi2   = "<<chi2<<endl
	  <<" dist   = "<<dist<<endl
	  <<" rdifphi= "<<rdiffphi<<endl
	  <<" diffr  = "<<diffr<<endl
	  <<" diffz  = "<<diffz<<endl
	  <<endl
	  <<" PhiRZ1 = "<<Form("%.2f  %.2f  %.2f",PhiRZ1[0],PhiRZ1[1],PhiRZ1[2])<<endl
	  <<" PhiRZ2 = "<<Form("%.2f  %.2f  %.2f",PhiRZ2[0],PhiRZ2[1],PhiRZ2[2])<<endl
	  <<" difphi = "<<diffphi<<" degree"<<endl
	  <<endl;
      break;
    }
  }

  float ymin3 = 0, ymax3 = 0, ymin4 = -0.01, ymax4 = 0;
  float chi2s[37];
  for(int iseg=0; iseg<37; iseg++){
    chi2s[iseg] = 0;
    for(int isig=0; isig<121; isig++){
      y1[iseg*121+isig] = pulse1[iseg][isig];
      y2[iseg*121+isig] = pulse2[iseg][isig];
      y3[iseg*121+isig] = chis[iseg][isig];
      if(y3[iseg*121+isig]>ymax3) ymax3 = y3[iseg*121+isig];
      chi2s[iseg] += chis[iseg][isig];
    }
    if(chi2s[iseg]==0) chi2s[iseg] = -1;
  }

  float chi2ave = 0;
  float chi2max = 0;
  int nfired2 = 0;
  for(int iseg=0; iseg<37; iseg++){
    for(int isig=0; isig<121; isig++)
      y4[iseg*121+isig] = chi2s[iseg];

    if(chi2s[iseg]>ymax4) ymax4 = chi2s[iseg];
    if(chi2s[iseg]>0){
      chi2ave += chi2s[iseg];
      nfired2++;
    }
    if(chi2s[iseg]>chi2max) chi2max=chi2s[iseg];
  }
  if(nfired2>0) chi2ave = chi2ave/nfired2;
  cout<<"calc:"<<endl
      <<" nfired = "<<nfired2<<endl
      <<" chi2   = "<<chi2max<<endl
      <<endl;

  ymin3 = ymin3 - 0.001;
  ymax3 = 1.1*ymax3 + 0.001;
  ymax4 = 1.1*ymax4 + 0.001;

  TPad *p1 = new TPad("p1","pulse",0,0,1.,0.7);
  p1->SetMargin(0.05,0.05,0.05,0.01);
  p1->Draw();
  p1->cd();
  gr1 = new TGraph(4477,x,y1);
  gr1->GetXaxis()->SetRangeUser(xmin,xmax);
  gr1->GetYaxis()->SetRangeUser(-0.3,1.1);
  gr1->SetLineWidth(2);
  gr1->SetTitle("");
  

  gr2 = new TGraph(4477,x,y2);
  gr2->SetLineColor(2);
  gr2->SetLineWidth(2);

  gr1->Draw("APL");
  gr2->Draw("Lsame");

  c->cd();
  TPad *p2 = new TPad("p2","chis",0,0.7,1.,1.);
  p2->SetMargin(0.05,0.05,0.15,0.05);
  p2->Draw();
  p2->cd();
  gr3 = new TGraph(4477,x,y3);
  gr3->GetXaxis()->SetRangeUser(xmin,xmax);
  gr3->GetYaxis()->SetRangeUser(ymin3,ymax3);
  gr3->GetXaxis()->SetLabelSize(0.08);
  gr3->GetYaxis()->SetLabelSize(0.08);
  gr3->GetYaxis()->SetNdivisions(505);
  gr3->SetLineWidth(2);
  gr3->SetTitle("");
  gr3->Draw("APL");

  float scale = (ymax3 - ymin3)/(ymax4-ymin4);
  for(int ii=0; ii<4477; ii++) y4[ii] = scale*y4[ii];
  gr4 = new TGraph(4477,x,y4);
  gr4->SetLineColor(4);
  gr4->SetMarkerColor(4);
  gr4->Draw("PLsame");
  
  TGaxis *axis = new TGaxis(xmax,ymin3,xmax,ymax3,ymin4,ymax4,502,"-L");
  axis->SetLineColor(4);
  axis->SetLabelColor(4);
  axis->SetLabelSize(0.08);
  axis->SetLabelOffset(-0.025);
  axis->Draw();
  
}
