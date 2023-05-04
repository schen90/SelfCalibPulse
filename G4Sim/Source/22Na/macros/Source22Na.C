//./Source22Na irun statistics

#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

TString Path="Source/22Na";

using namespace std;

void Source22Na(int irun, int stat){
  float Egamma1 = 1274.6; // keV
  float Egamma2 = 511;    // keV
  TVector3 vsource(0,0,0);
  
  TRandom3 r(0);
  const Float_t Pi = TMath::Pi();
  
  cout<<"run "<<irun<<" simulate "<<stat<<" events..."<<endl;

  ofstream fout(Form("%s/22NaEvents%04d",Path.Data(),irun));
  fout<<"FORMAT 4 0"<<endl
      <<"REACTION 1 1 6 12 0.0"<<endl
      <<"EMITTED 3 1 1 1"<<endl;

  for(int i=0; i<stat; i++){

    if(i%10000==0) cout<<"\r finish "<<i<<" / "<<stat<<" ..."<<flush;

    fout<<"$"<<endl
	<<"-101"<<endl;

    // 1274.6 keV gamma
    double ThetaGamma1 = acos(r.Uniform(-1,1));
    double PhiGamma1 = r.Uniform(-Pi,Pi);
    TVector3 vec1;
    vec1.SetMagThetaPhi(1, ThetaGamma1, PhiGamma1);
    fout<< " 1 "
	<< setw(7) << Form("%.1f ",Egamma1)
	<< setw(10) << Form("%.6f ",vec1.X())
	<< setw(10) << Form("%.6f ",vec1.Y())
	<< setw(10) << Form("%.6f ",vec1.Z())
	<< Form("%.2f ",vsource.X())
	<< Form("%.2f ",vsource.Y())
	<< Form("%.2f ",vsource.Z())
	<< endl;

    // 511 keV gammas
    double ThetaGamma2 = acos(r.Uniform(-1,1));
    double PhiGamma2 = r.Uniform(-Pi,Pi);
    TVector3 vec2;
    vec2.SetMagThetaPhi(1, ThetaGamma2, PhiGamma2);
    fout<< " 1 "
	<< setw(7) << Form("%.1f ",Egamma2)
	<< setw(10) << Form("%.6f ",vec2.X())
	<< setw(10) << Form("%.6f ",vec2.Y())
	<< setw(10) << Form("%.6f ",vec2.Z())
	<< Form("%.2f ",vsource.X())
	<< Form("%.2f ",vsource.Y())
	<< Form("%.2f ",vsource.Z())
	<< endl;

    fout<< " 1 "
	<< setw(7) << Form("%.1f ",Egamma2)
	<< setw(10) << Form("%.6f ",-vec2.X())
	<< setw(10) << Form("%.6f ",-vec2.Y())
	<< setw(10) << Form("%.6f ",-vec2.Z())
	<< Form("%.2f ",vsource.X())
	<< Form("%.2f ",vsource.Y())
	<< Form("%.2f ",vsource.Z())
	<< endl;
  }
  cout<<endl;
  
  fout.close();

  return;
}

#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc<3){ cout<<"missing parameters!!!"<<endl; return 1;}

  Source22Na(atoi(argv[1]),atoi(argv[2]));
  
  return 0;
}
#endif
