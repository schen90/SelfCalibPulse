// ./macros/CalcMapPointsGrid

#include <TRint.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

const Int_t NCryType = 3;
double mindist = 0.1;
double griddist = 2; //2mm grid
double range[NCryType][3][2]; // range of XYZ

void swap(vector<TVector3> &v, int m, int l){
  TVector3 temp;

  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
}

void swap(vector<int> &v, int m, int l){
  int temp;

  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
}

void CalcMapPointsGrid(){

  // clock
  time_t start, stop;
  time(&start);

  string filename[NCryType] = {"LibTrap_A001.root","LibTrap_B001.root","LibTrap_C001.root"};
    
  string inputfiledir = "G4Sim/pulsedb/";
  string outputfile = "Map/MapPointsGrid.dat";

  ofstream fout(outputfile.c_str(),ios::out);
  fout<<"#range #####################"<<endl;
  fout<<"# itype  xmin  xmax  xsteps  ymin  ymax  ysteps  zmin  zmax  zsteps #####################"<<endl;
  for(int itype=0; itype<NCryType; itype++){
    fout<<Form("  %d",itype);
    for(int iaxis=0; iaxis<3; iaxis++){
      fout<<Form("  min%d_%d  max%d_%d  steps%d_%d",iaxis,itype,iaxis,itype,iaxis,itype);
      range[itype][iaxis][0] = 1000;
      range[itype][iaxis][1] = -1000;
    }
    fout<<endl;
  }
  fout<<" -1  -1  -1  -1  -1  -1  -1  -1  -1  -1"<<endl;
  
  fout<<"#Map #####################"<<endl;
  fout<<"# itype seg  x  y  z  chi2s[0]  chi2s[1]  chi2s[2]#####################"<<endl;

  for(int itype=0; itype<NCryType; itype++){
    string inputfile = inputfiledir + filename[itype];
    TFile *fin = new TFile(inputfile.c_str());
    TTree *tree = (TTree *)fin->Get("tree");
    int seg;
    float pos[3];
    tree->SetBranchAddress("seg",&seg);
    tree->SetBranchAddress("pos",pos);
    int nentries = tree->GetEntriesFast();

    vector<int> gseg;
    vector<TVector3> gpos;
    int ievt;
    for(ievt=0; ievt<nentries; ievt++){
      if(ievt%1000==0)  cout<<"\r type "<<itype<<": finish "<<ievt<<"/"<<nentries<<" entries... "<<flush;
      tree->GetEntry(ievt);

      gseg.push_back(seg);
      TVector3 tmppos(pos[0],pos[1],pos[2]);
      gpos.push_back(tmppos);
      for(int iaxis=0; iaxis<3; iaxis++){
	if(pos[iaxis]<range[itype][iaxis][0]) range[itype][iaxis][0] = pos[iaxis];
	if(pos[iaxis]>range[itype][iaxis][1]) range[itype][iaxis][1] = pos[iaxis];
      }
    }
    cout<<"\r type "<<itype<<": finish "<<ievt<<"/"<<nentries<<" entries... "<<endl;

    // sort
    cout<<"sort type "<<itype<<" points..."<<endl;
    for(int ipoint=0; ipoint<gpos.size(); ipoint++){
      for(int jpoint=ipoint+1; jpoint<gpos.size(); jpoint++){
	if(gpos[ipoint].Z()-gpos[jpoint].Z()<-mindist)
	  continue;
	else if(gpos[ipoint].Z()-gpos[jpoint].Z()<mindist)
	  if(gpos[ipoint].Y()-gpos[jpoint].Y()<-mindist)
	    continue;
	  else if(gpos[ipoint].Y()-gpos[jpoint].Y()<mindist)
	    if(gpos[ipoint].X()-gpos[jpoint].X()<-mindist)
	      continue;
	  
	swap(gseg,ipoint,jpoint);	
	swap(gpos,ipoint,jpoint);	
      }
    }

    //output
    for(int ipoint=0; ipoint<gpos.size(); ipoint++){
      double sigma_phi = 6;
      double sigma_r = 1;
      double sigma_z = 1.5;
      fout<<Form("  %d  %d   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f",itype,gseg[ipoint],gpos[ipoint].X(),gpos[ipoint].Y(),gpos[ipoint].Z(),-1.,-1.,-1.)<<endl;
    }
    
    fin->Close();
  }
  fout<<" -1  -1  -1  -1  -1  -1  -1  -1"<<endl;
  fout.close();

  // output range
  for(int itype=0; itype<NCryType; itype++){
    for(int iaxis=0; iaxis<3; iaxis++){
      int steps = (int)((range[itype][iaxis][1]-range[itype][iaxis][0])/griddist+0.5)+1;
      gROOT->ProcessLine(Form(".!sed -i 's/min%d_%d/%.3f/g' %s",iaxis,itype,range[itype][iaxis][0],outputfile.c_str()));
      gROOT->ProcessLine(Form(".!sed -i 's/max%d_%d/%.3f/g' %s",iaxis,itype,range[itype][iaxis][1],outputfile.c_str()));
      gROOT->ProcessLine(Form(".!sed -i 's/steps%d_%d/%d/g' %s",iaxis,itype,steps,outputfile.c_str()));
    }
  }

  time(&stop);
  cout<<Form("============= Elapsed time: %.1f seconds =============",difftime(stop,start))<<endl;
  return;
}

#ifndef __CINT__
int main(int argc, char *argv[]){
    
  CalcMapPointsGrid();
  
  return 0;
}
#endif
