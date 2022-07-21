#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TH1.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

bool kextrapol = true;

const int nsig = 56;
const int nseg = 36;
const int nsegcore = 37;

// selfcalib
vector<int> calibseg;
vector<TMatrixD> calibpos;
vector<TMatrixD> calibpulse;


// main
void MakePSBasePSC(string PSCfile, string Basefile){

  // load calib pscfile
  TChain *tree = new TChain();
  for(int iseg=0; iseg<nseg; iseg++)
    tree->AddFile(PSCfile.c_str(),0,Form("tree%d",iseg));

  int det, seg;
  float cadpos2[3];
  float spulse[nsegcore][nsig];
  int npaths;
  tree->SetBranchAddress("det",&det);
  tree->SetBranchAddress("seg",&seg);
  //tree->SetBranchAddress("detpos",cadpos2);
  tree->SetBranchAddress("cadpos2",cadpos2);
  tree->SetBranchAddress("spulse",spulse);
  tree->SetBranchAddress("npaths",&npaths);

  int nentries = tree->GetEntriesFast();

  // read PSCfile
  int ientry;
  for(ientry=0; ientry<nentries; ientry++){
    if(ientry%1000==0) cout<<"\r finish "<<ientry<<" / "<<nentries<<" points from PSC"<<flush;
    tree->GetEntry(ientry);
    //seg = seg+1; // seg in db start from 1...
    
    if(npaths<=100) continue;

    TMatrixD tmppos(3,1);
    for(int ix=0; ix<3; ix++){
      tmppos(ix,0)=cadpos2[ix];
    }
    TMatrixD tmppulse(nsegcore*nsig,1);
    for(int iseg=0; iseg<nsegcore; iseg++)
      for(int isig=0; isig<nsig; isig++)
	tmppulse(iseg*nsig+isig,0) = spulse[iseg][isig];

    calibseg.push_back(seg);
    calibpos.push_back(tmppos);
    calibpulse.push_back(tmppulse);

  }
  cout<<"\r finish "<<ientry<<" / "<<nentries<<" points from PSC"<<endl;

  // output
  float pos[3];
  float ospulse[nsig*nsegcore];
  
  TFile *fout = new TFile(Basefile.c_str(),"RECREATE");
  TTree *outtree = new TTree("tree","SelfCalib Base");

  outtree->Branch("det",&det);
  outtree->Branch("seg",&seg);
  outtree->Branch("pos",pos,"pos[3]/F");
  outtree->Branch("spulse",ospulse,Form("spulse[%d]/F",nsig*nsegcore));

  // output at PSC position
  int npos = calibpos.size();
  int ipos;
  for(ipos=0; ipos<npos; ipos++){
    if(ipos%1000==0) cout<<"\r finish "<<ipos<<" / "<<npos<<" PSC points"<<flush;
    seg = calibseg[ipos];

    for(int ix=0; ix<3; ix++){
      pos[ix] = calibpos[ipos](ix,0);
    }

    for(int isig=0; isig<nsig*nsegcore; isig++) ospulse[isig] = calibpulse[ipos](isig,0);
    
    outtree->Fill();
  }
  cout<<"\r finish "<<ipos<<" / "<<npos<<" PSC points"<<endl;
    
  fout->cd();
  outtree->Write();
  fout->Close();

  return;
}


#ifndef __CINT__
int main(int argc, char *argv[]){
  string PSCfilename = "PSCfiles/Det0000.root";
  string PSBasename = "PSBase/Det0000.root";

  if(argc>2) PSBasename = string(argv[2]);
  if(argc>1) PSCfilename = string(argv[1]);

  MakePSBasePSC(PSCfilename, PSBasename);
  return 0;
}
#endif
