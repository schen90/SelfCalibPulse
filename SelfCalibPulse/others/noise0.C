#include "TCanvas.h"
#include "TRandom.h"
#include "iostream"
#include "TVirtualFFT.h"

using namespace std;

const int NSig = 121;

void noise0(){
  const int nevts = 1;

  TH1D *h = new TH1D("h","h",NSig-1,0,600);

  int ievt;
  for(ievt=0; ievt<nevts; ievt++){
    if(ievt%100==0) cout<<"\r finish "<<ievt<<" / "<<nevts<<" evts..."<<flush;
    
    // get random noise
    float noise[NSig],noise2[NSig];
    float tmpnoise;
    for(int isig=0; isig<NSig; isig++){
      tmpnoise = gRandom->Gaus(0,3);
      noise2[isig] = noise[isig] = tmpnoise;
    }

    // add noise
    for(int isig=0; isig<NSig; isig++)
      h->SetBinContent(isig,noise2[isig]);

  }
  cout<<"\r finish "<<ievt<<" / "<<nevts<<" evts..."<<endl;

  TCanvas *c = new TCanvas("c","c",600,900);
  c->Divide(1,3);
  c->cd(1);
  h->Draw("h");

  c->cd(2);
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = h->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");
  hm->Draw();

  
  c->cd(3);
  TH1 *hp =0;
  hp = h->FFT(hp, "PH");
  hp->SetTitle("Phase of the 1st transform");
  hp->Draw();

  //return;

  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  Double_t *re_full = new Double_t[NSig];
  Double_t *im_full = new Double_t[NSig];
  fft->GetPointsComplex(re_full,im_full);


  // filter ---------------------------
  /*
  for(int isig=0; isig<NSig; isig++){
    if(isig<20) continue;
    re_full[isig] = 0;
  }
  for(int isig=0; isig<NSig; isig++){
    if(isig<20) continue;
    im_full[isig] = 0;
  }
  */
  //-----------------------------------

  
  c->cd(2);
  TH1D *h2 = new TH1D("h2","re",NSig,0,NSig);
  for(int isig=0; isig<NSig; isig++) h2->SetBinContent(isig,re_full[isig]);
  h2->Draw();

  c->cd(3);
  TH1D *h3 = new TH1D("h3","im",NSig,0,NSig);
  for(int isig=0; isig<NSig; isig++) h3->SetBinContent(isig,im_full[isig]);
  h3->Draw();

  return;
  
  // backward transform
  c->cd(3);
  int nsig = NSig;
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &nsig, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  hb = TH1::TransformHisto(fft_back,hb,"Re");
  hb->SetTitle("The backward transform");
  c->cd(3);
  hb->Draw("h");
  

  
  return;
}
