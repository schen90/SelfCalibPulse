#include "TCanvas.h"
#include "TRandom.h"
#include "iostream"
#include "TVirtualFFT.h"
#include "TH1.h"
#include "TGraph.h"

using namespace std;

const int NSegCore = 37;
const int NSig = 121;
const int npoint = 4477;

void pulse(){
  const int nevts = 1;

  TH1D *h = new TH1D("h","h",npoint,0,npoint);
  TFile *f = new TFile("pulse.root");
  TGraph *gr = (TGraph *)f->Get("Graph");
  
  int ievt;
    
  // get random noise
  float noise[NSegCore][NSig],noise2[NSegCore][NSig];
  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++){
      int ipoint = iseg*NSig+isig;
      noise2[iseg][isig] = noise[iseg][isig] = gr->GetPointY(ipoint);
    }
  }


  // add noise
  for(int iseg=0; iseg<NSegCore; iseg++)
    for(int isig=0; isig<NSig; isig++)
      h->SetBinContent(iseg*NSig+isig,noise2[iseg][isig]);



  TCanvas *c = new TCanvas("c","c",500,800);
  c->Divide(1,2);
  c->cd(1);
  h->Draw("h");

  c->cd(2);
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = h->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");
  hm->Draw();

  
  return;
}
