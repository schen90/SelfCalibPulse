#include "TCanvas.h"
#include "TRandom.h"
#include "iostream"
#include "TVirtualFFT.h"
#include "time.h"

using namespace std;

const int NSig = 121;

TCanvas *c = new TCanvas("c","c",600,400);

void Dev(){
  const int nevts = 1;

  TH1D *h = new TH1D("h","h",NSig-1,0,600);
  TH1 *hb = 0;
  
  time_t t;
  gRandom->SetSeed(time(&t));
  
  int ievt;
  for(ievt=0; ievt<nevts; ievt++){
    if(ievt%100==0) cout<<"\r finish "<<ievt<<" / "<<nevts<<" evts..."<<flush;
    
    // get random noise

    // FFT
    Double_t re_full[NSig];
    Double_t im_full[NSig];
    double k = -0.5;
    for(int isig=0; isig<NSig; isig++){
      if(isig==0 || isig>6){
	re_full[isig] = im_full[isig] = 0;
	continue;
      }
      re_full[isig] = gRandom->Uniform(-1,1)*exp(k*isig)/40;
      im_full[isig] = gRandom->Uniform(-1,1)*exp(k*isig)/40;

      if(isig==1){
	re_full[isig] = re_full[isig]/20;
	im_full[isig] = im_full[isig]/20;
      }else if(isig==2 || isig==3){
	re_full[isig] = re_full[isig]*2;
	im_full[isig] = im_full[isig]*2;
      }

      re_full[0] += -2*re_full[isig];
      im_full[0] += -2*im_full[isig];
    }

    // FFT backward transform
    int nsig = NSig;
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &nsig, "C2R M K");
    fft_back->SetPointsComplex(re_full,im_full);
    fft_back->Transform();
    hb = TH1::TransformHisto(fft_back,hb,"Re");    

    for(int isig=0; isig<NSig; isig++){
      double bincontent = hb->GetBinContent(isig);
      h->SetBinContent(isig,bincontent);
    }
  }
  cout<<"\r finish "<<ievt<<" / "<<nevts<<" evts..."<<endl;

  h->GetYaxis()->SetRangeUser(-0.1,0.1);
  h->Draw("h");

  return;
}
