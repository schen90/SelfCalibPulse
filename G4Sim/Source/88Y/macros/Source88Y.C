//./Source88Y irun statistics

#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

TString Path="Source/88Y";

using namespace std;

const double Pi = TMath::Pi();

double W3D(double *x, double *par){
  double theta = x[0];
  double alpha2 = par[0];
  double alpha4 = par[1];

  double costheta  = cos(theta);
  double costheta2 = costheta * costheta;
  double costheta4 = costheta2 * costheta2;
  double sintheta  = sin(theta);

  double W = ( 1 + alpha2*costheta2 + alpha4*costheta4 ) * fabs(sintheta);
  return W;
}


void Source88Y(int irun, int stat){
  float Egamma1 = 1836; // keV, 5% direct population
  float Egamma2 = 898;  // keV, 95% direct population
  TVector3 vsource(0,0,0);
  
  TRandom3 r(0);
  TF1 *fThetaCor = new TF1("fThetaCor",W3D,0,Pi,2);
  fThetaCor->SetParameter(0, -3./29);
  fThetaCor->SetParameter(1, 0);
  
  cout<<"run "<<irun<<" simulate "<<stat<<" events..."<<endl;

  ofstream fout(Form("%s/88YEvents%04d",Path.Data(),irun));
  fout<<"FORMAT 4 0"<<endl
      <<"REACTION 1 1 6 12 0.0"<<endl
      <<"EMITTED 3 1 1 1"<<endl;

  for(int i=0; i<stat; i++){

    if(i%10000==0) cout<<"\r finish "<<i<<" / "<<stat<<" ..."<<flush;

    fout<<"$"<<endl
	<<"-101"<<endl;

    // population
    double population = r.Uniform(0,100);

    if(population<5){ // 2+ -> 0+

      // 1836 keV gamma
      double ThetaGamma1 = acos(r.Uniform(-1,1));
      double PhiGamma1 = r.Uniform(-Pi,Pi);
      TVector3 vec1;
      vec1.SetMagThetaPhi(1., ThetaGamma1, PhiGamma1);
      fout<< " 1 "
	  << setw(7) << Form("%.1f ",Egamma1)
	  << setw(10) << Form("%.6f ",vec1.X())
	  << setw(10) << Form("%.6f ",vec1.Y())
	  << setw(10) << Form("%.6f ",vec1.Z())
	  << Form("%.2f ",vsource.X())
	  << Form("%.2f ",vsource.Y())
	  << Form("%.2f ",vsource.Z())
	  << endl;

    }else{ // 3- -> 2+ -> 0+

      // 898 keV gamma
      double ThetaGamma2 = acos(r.Uniform(-1,1));
      double PhiGamma2 = r.Uniform(-Pi,Pi);
      TVector3 vec2;
      vec2.SetMagThetaPhi(1., ThetaGamma2, PhiGamma2);
      fout<< " 1 "
	  << setw(7) << Form("%.1f ",Egamma2)
	  << setw(10) << Form("%.6f ",vec2.X())
	  << setw(10) << Form("%.6f ",vec2.Y())
	  << setw(10) << Form("%.6f ",vec2.Z())
	  << Form("%.2f ",vsource.X())
	  << Form("%.2f ",vsource.Y())
	  << Form("%.2f ",vsource.Z())
	  << endl;

      // angular correlation
      double ThetaCor = fThetaCor->GetRandom();
      double PhiCor = r.Uniform(-Pi,Pi);
      TVector3 vec1;
      vec1.SetMagThetaPhi(1., ThetaCor, PhiCor);
      vec1.RotateY(vec2.Theta());
      vec1.RotateZ(vec2.Phi());

      // 1836 keV gamma
      fout<< " 1 "
	  << setw(7) << Form("%.1f ",Egamma1)
	  << setw(10) << Form("%.6f ",vec1.X())
	  << setw(10) << Form("%.6f ",vec1.Y())
	  << setw(10) << Form("%.6f ",vec1.Z())
	  << Form("%.2f ",vsource.X())
	  << Form("%.2f ",vsource.Y())
	  << Form("%.2f ",vsource.Z())
	  << endl;

    }

  }
  cout<<endl;
  
  fout.close();

  return;
}

#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc<3){ cout<<"missing parameters!!!"<<endl; return 1;}

  Source88Y(atoi(argv[1]),atoi(argv[2]));
  
  return 0;
}
#endif
