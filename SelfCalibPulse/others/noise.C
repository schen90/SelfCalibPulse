#include "TCanvas.h"
#include "TRandom.h"
#include "iostream"

using namespace std;

const int NSegCore = 37;
const int NSig = 121;

void noise(){
  int nevts = 10000;

  TH1D *h = new TH1D("h","h",200,-20,20);

  int ievt;
  for(ievt=0; ievt<nevts; ievt++){
    if(ievt%100==0) cout<<"\r finish "<<ievt<<" / "<<nevts<<" evts..."<<flush;
    
    // get random noise
    float noise[NSegCore][NSig],noise2[NSegCore][NSig];
    float tmpnoise;
    for(int iseg=0; iseg<NSegCore; iseg++){
      tmpnoise = gRandom->Uniform(-1,1);
      for(int isig=0; isig<NSig; isig++){
	tmpnoise += -0.2*(tmpnoise+gRandom->Uniform(-10,10));
	noise2[iseg][isig] = noise[iseg][isig] = tmpnoise;
      }
    }

    // smooth noise
    int nsmooth = 3;
    for(int iseg=0; iseg<NSegCore; iseg++)
      for(int isig=nsmooth; isig<NSig-nsmooth; isig++){
	for(int i=1; i<=nsmooth; i++)
	  noise2[iseg][isig] += noise[iseg][isig-i] + noise[iseg][isig+i];

	noise2[iseg][isig] = noise2[iseg][isig] / (2*nsmooth+1) * 2;
      }

    // add noise
    for(int iseg=0; iseg<NSegCore; iseg++)
      for(int isig=0; isig<NSig; isig++)
	h->Fill(noise2[iseg][isig]);

  }
  cout<<"\r finish "<<ievt<<" / "<<nevts<<" evts..."<<endl;

  TCanvas *c = new TCanvas("c","c");
  h->Draw();
  
  return;
}
