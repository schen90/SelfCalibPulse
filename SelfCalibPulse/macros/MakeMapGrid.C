// ./macros/MakeMapGrid ichi2 PSfile mapname

#include <TRint.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TProfile.h>
#include "TSystem.h"

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

double maxdist = 1; // 1mm
double mindist = 0.1;

// segmentation
const Int_t NCryType = 3;
double griddist = 2; // 2mm grid
double range[NCryType][3][2]; // range of XYZ
const Int_t MaxSteps = 50;

double Map[NCryType][MaxSteps][MaxSteps][MaxSteps][3]; // chi2s limit map (chi2s_r,chi2s_phi,chi2s_z)
double MapPos[NCryType][MaxSteps][MaxSteps][MaxSteps][3]; // position of map (x,y,z)
int    Seg[NCryType][MaxSteps][MaxSteps][MaxSteps]; // segment at grid point

// read map
void LoadMapGrid(string mapfilename){
  ifstream fin(mapfilename.c_str());
  if(!fin){
    gROOT->ProcessLine(Form(".!cp Map/MapPointsGrid.dat %s",mapfilename.c_str()));
    fin.open(mapfilename.c_str());
  }

  cout<<"\e[1;32m read position resolution map \e[0m"<<endl;

  // init
  for(int itype=0; itype<NCryType; itype++){
    for(int iaxis=0; iaxis<3; iaxis++){
      range[itype][iaxis][0] = 1;
      range[itype][iaxis][1] = 0;
    }
    for(int ix=0; ix<MaxSteps; ix++)
      for(int iy=0; iy<MaxSteps; iy++)
	for(int iz=0; iz<MaxSteps; iz++){
	  Map[itype][ix][iy][iz][0] = -100;
	  Map[itype][ix][iy][iz][1] = -100;
	  Map[itype][ix][iy][iz][2] = -100;
	  Seg[itype][ix][iy][iz] = -1;
	}
  }
  
  // read Map
  const int kMaxBufLen = 500;
  char *buffer = new char[kMaxBufLen];
  while(!fin.eof()){
    fin.getline(buffer,kMaxBufLen);

    if(strncmp(buffer,"#range",6)==0){
      cout<<"reading range for each crystal type"<<endl;
      fin.getline(buffer,kMaxBufLen);
      int itype;
      while(1){
        fin >> itype;
        if(itype==-1) break;
	int tmp;
        for(int iaxis=0; iaxis<3; iaxis++){
	  fin >> range[itype][iaxis][0] >> range[itype][iaxis][1] >> tmp;
	  if(tmp>MaxSteps) {cerr<<"change MaxSteps to >= "<<tmp<<endl; return;}
	}
      }
    }

    if(strncmp(buffer,"#Map",4)==0){
      cout<<"reading Map"<<endl;
      fin.getline(buffer,kMaxBufLen);
      int itype, iseg, ipos[3];
      double pos[3],res[3];
      while(1){
        fin >> itype >> iseg >> pos[0] >> pos[1] >> pos[2];
        if(itype==-1) break;

	for(int iaxis=0; iaxis<3; iaxis++){
	  if(pos[iaxis]-range[itype][iaxis][0]<-mindist ||
	     pos[iaxis]-range[itype][iaxis][1]>mindist){
	    cout<<Form("axis%d: %.3f outside range %.3f ~ %.3f",iaxis,pos[iaxis],range[itype][iaxis][0],range[itype][iaxis][1])<<endl;
	    ipos[iaxis] = -1;
	  }else{
	    ipos[iaxis] = (int)((pos[iaxis] - range[itype][iaxis][0]) / griddist + 0.5);
	    if(fabs(pos[iaxis]-(range[itype][iaxis][0]+griddist*ipos[iaxis]))>mindist){
	      cout<<Form("axis%d: %.3f not a grid point",iaxis,pos[iaxis])<<endl;
	      ipos[iaxis] = -1;
	    }
	  }
	}
	
        for(int i=0; i<3; i++) fin>>res[i];
	if(ipos[0]<0 || ipos[1]<0 || ipos[2]<0) continue;
        for(int i=0; i<3; i++){
	  Map[itype][ipos[0]][ipos[1]][ipos[2]][i] = res[i];
	  MapPos[itype][ipos[0]][ipos[1]][ipos[2]][i] = pos[i];
	}
	Seg[itype][ipos[0]][ipos[1]][ipos[2]] = iseg;
      }
    }

  }
  cout<<"finish load Map"<<endl;

  fin.close();
  return;
}

// write map
void WriteMapGrid(string mapfilename){

  cout<<"Write map..."<<endl;
  ofstream fout("Map/tmp.dat",ios::out);
  
  //output range
  fout<<"#range #####################"<<endl;
  fout<<"# itype  xmin  xmax  xsteps  ymin  ymax  ysteps  zmin  zmax  zsteps #####################"<<endl;
  for(int itype=0; itype<NCryType; itype++){
    fout<<Form("  %d",itype);
    for(int iaxis=0; iaxis<3; iaxis++){
      int steps = (int)((range[itype][iaxis][1]-range[itype][iaxis][0])/griddist+0.5)+1;
      fout<<Form("  %.3f  %.3f  %d",range[itype][iaxis][0],range[itype][iaxis][1],steps);
    }
    fout<<endl;
  }
  fout<<" -1  -1  -1  -1  -1  -1  -1  -1  -1  -1"<<endl;

  // output map
  fout<<"#Map #####################"<<endl;
  fout<<"# itype seg  x  y  z  chi2s[0]  chi2s[1]  chi2s[2]#####################"<<endl;
  for(int itype=0; itype<NCryType; itype++)
    for(int iz=0; iz<MaxSteps; iz++)
      for(int iy=0; iy<MaxSteps; iy++)
	for(int ix=0; ix<MaxSteps; ix++){

            if(Seg[itype][ix][iy][iz]<0) continue;

	    fout<<Form("  %d  %d   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",
		       itype,
		       Seg[itype][ix][iy][iz],
		       MapPos[itype][ix][iy][iz][0],
		       MapPos[itype][ix][iy][iz][1],
		       MapPos[itype][ix][iy][iz][2],
		       Map[itype][ix][iy][iz][0],
		       Map[itype][ix][iy][iz][1],
		       Map[itype][ix][iy][iz][2])<<endl;
          }
  fout<<" -1  -1  -1  -1  -1  -1  -1  -1"<<endl;

  fout.close();  
  gROOT->ProcessLine(Form(".!mv -f Map/tmp.dat %s",mapfilename.c_str()));
  return;
}


// make Map
void MakeMapGrid(int itype, int ichi2, string PSfile){
  if(itype<0 || itype>2) return;
  if(ichi2<0 || ichi2>2) return;
  
  // get chi2s limit
  TFile *f = new TFile(PSfile.c_str());
  TTree *tree = (TTree *)f->Get(Form("tree%d",itype));
  
  Int_t simseg;
  Float_t SimPos[3];
  Float_t chi2s[3];
  Float_t diffr;
  Float_t diffz;
  Float_t rdiffphi;
  tree->SetBranchAddress("simseg",&simseg);
  tree->SetBranchAddress("SimPos",SimPos);
  tree->SetBranchAddress("chi2s",chi2s);
  tree->SetBranchAddress("diffr",&diffr);
  tree->SetBranchAddress("diffz",&diffz);
  tree->SetBranchAddress("rdiffphi",&rdiffphi);

  int nentries = tree->GetEntriesFast();
  cout<<"type "<<itype<<": read "<<nentries<<" events from tree"<<endl;

  int nbins = MaxSteps*MaxSteps*MaxSteps;
  TProfile *pchi2;
  if(ichi2==0) pchi2 = new TProfile("pchi2","chi2s limit [r]"  ,nbins,0.5,nbins+0.5,0,2);
  if(ichi2==1) pchi2 = new TProfile("pchi2","chi2s limit [phi]",nbins,0.5,nbins+0.5,0,2);
  if(ichi2==2) pchi2 = new TProfile("pchi2","chi2s limit [z]"  ,nbins,0.5,nbins+0.5,0,2);
  TH1D *htmp = new TH1D("htmp","ibin",nbins,0.5,nbins+0.5);
  
  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)
      cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;
    tree->GetEntry(ievt);

    if(ichi2==0) if(diffr<1)    continue;
    if(ichi2==1) if(rdiffphi<1) continue;
    if(ichi2==2) if(diffz<1)    continue;
    
    for(int iz=0; iz<MaxSteps; iz++){
      int nextz = 0;
      for(int iy=0; iy<MaxSteps; iy++){
	if(nextz){ nextz=0; break;}
	for(int ix=0; ix<MaxSteps; ix++){

	  if(Seg[itype][ix][iy][iz]<0) continue;

	  if(fabs(SimPos[2]-MapPos[itype][ix][iy][iz][2])>maxdist){ nextz=1; break;} //z
	  if(fabs(SimPos[1]-MapPos[itype][ix][iy][iz][1])>maxdist) break; //y
	  if(fabs(SimPos[0]-MapPos[itype][ix][iy][iz][0])>maxdist) continue; //x

	  int ibin =
	    iz*MaxSteps*MaxSteps +
	    iy*MaxSteps +
	    ix;

	  htmp->Fill(ibin);
	  pchi2->Fill(ibin,chi2s[ichi2]);
	} // end of loop x
      } // end of loop y
    } // end of loop z
    
  }// end of loop evts


  // update map for itype, ichi2
  cout<<"update map for type "<<itype<<" chi2s["<<ichi2<<"] ..."<<endl;
  for(int iz=0; iz<MaxSteps; iz++)
    for(int iy=0; iy<MaxSteps; iy++)
      for(int ix=0; ix<MaxSteps; ix++){

	if(Seg[itype][ix][iy][iz]<0) continue;

	int ibin =
	  iz*MaxSteps*MaxSteps +
	  iy*MaxSteps +
	  ix;

	if(pchi2->GetBinContent(ibin)>0)
	  Map[itype][ix][iy][iz][ichi2] = pchi2->GetBinContent(ibin);
	    
      }

  TFile *ftmp = new TFile("htmp.root","RECREATE");
  htmp->Write();
  pchi2->Write();
  ftmp->Close();

  return;
}


#ifndef __CINT__
int main(int argc, char *argv[]){
  // clock
  time_t start, stop;
  time(&start);

  int ichi2 = 1;
  string PSfile = "ComparePS3_Noise_D2.root";
  string Mapfile = "Map/MapGrid.dat";
  
  if(argc>3) Mapfile = string(argv[3]);
  if(argc>2) PSfile = string(argv[2]);
  if(argc>1) ichi2 = atoi(argv[1]);


  LoadMapGrid(Mapfile);

  for(int itype=0; itype<3; itype++){
    MakeMapGrid( itype, ichi2, PSfile);
  }

  WriteMapGrid(Mapfile);

  
  time(&stop);
  cout<<Form("============= Elapsed time: %.1f seconds =============",difftime(stop,start))<<endl;
  
  return 0;
}
#endif
