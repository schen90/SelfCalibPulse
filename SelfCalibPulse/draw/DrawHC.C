#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

// draw selection
int detid = 0;
int segid = 2;
int nhitslimit = 10;
int Number = 1;


struct HitCollection{
  int det;
  int seg;
  int nhits;

  float labpos[3];
  float detpos[3];

  vector<int> hid;
};

//typedef struct HitCollection HitCollection;

struct Hit{
  int det;
  int seg;
  float depE;

  float labpos[3];
  float detpos[3];
};

using namespace std;

// transform matrix
Int_t NDets = 0;
TMatrixD Rt[180];
TMatrixD Tr[180];

void LoadMatrix(string LookUpTable){
  // find input LookUpTable
  ifstream fin;
  int dummy_i;
  fin.open(LookUpTable.c_str());
  if(!fin){ cerr<<"Cannot find "<<LookUpTable<<endl; return;}
  cout<<"\e[1;32m find CrystalPosition from "<<LookUpTable<<"... \e[0m"<<endl;

  int ir;
  for(int i=0; i<180; i++){
    ir = -1;
    fin >> ir >> dummy_i;
    if(ir<0 || ir>=180){ ir=i-1; break;}

    Tr[ir].ResizeTo(3,1); Tr[ir].Zero();
    Rt[ir].ResizeTo(3,3); Rt[ir].Zero();
    for(int it=0; it<3; it++) fin >> Tr[ir](it,0);
    for(int it=0; it<3; it++){
      fin >> dummy_i;
      for(int it2=0; it2<3; it2++) fin >> Rt[ir](it,it2);
    }
    Rt[ir].Invert();  // change to rot from world frame -> detector frame
  }
  fin.close();
  NDets = ir+1;
  cout<<"read position for "<<NDets<<" detectors"<<endl;

  return;
}


TMatrixD Lab2DetPos(Int_t idet, TMatrixD LabPos){
  if(idet<0 || idet>=NDets){
    cerr<<"cannot find matrix for idet = "<<idet<<endl;
    return LabPos;
  }
  TMatrixD DetPos(3,1);
  DetPos = Rt[idet]*(LabPos-Tr[idet]);
  return DetPos;
}


void DrawHC(){
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");

  LoadMatrix("LookUp/CrystalPositionLookUpTable");
  
  vector<HitCollection> fAllHCs; // Hit Collection storage
  vector<Hit> fHits; //fHits in Events for tracking

  int det;
  int seg;

  float labpos[3];

  // HitCollection
  TFile *hcfile = new TFile("HCfiles/AllHCs.root");
  TTree *hctree = (TTree *)hcfile->Get("hctree");

  hctree->SetBranchAddress("det",&det);
  hctree->SetBranchAddress("seg",&seg);

  hctree->SetBranchAddress("labpos",labpos); // real position in lab frame

  int Nhcs = hctree->GetEntriesFast();
  int ihc;
  for(ihc=0; ihc<Nhcs; ihc++){
    if(ihc%1000==0) cout<<"\r load "<<ihc<<" / "<<Nhcs<<" HitCollections..."<<flush;
    hctree->GetEntry(ihc);

    HitCollection ahc;
    ahc.det = det;
    ahc.seg = seg;
    ahc.nhits = 0;

    TMatrixD LabPos(3,1);
    for(int it=0; it<3; it++) LabPos(it,0) = labpos[it];
    TMatrixD DetPos= Lab2DetPos(det,LabPos);

    for(int it=0; it<3; it++){
      ahc.labpos[it] = LabPos(it,0);
      ahc.detpos[it] = DetPos(it,0);
    }

    fAllHCs.push_back(ahc);
  }
  cout<<"\r load "<<ihc<<" / "<<Nhcs<<" HitCollections..."<<endl;

  // Hit
  TFile *hfile = new TFile("HCfiles/EventHits.root");
  TTree *htree = (TTree *)hfile->Get("htree");
  vector<int>            *vdet = 0;
  vector<int>            *vseg = 0;
  vector<vector<int>>    *vhcid = 0;
  vector<float>          *vdepE = 0;
  vector<vector<double>> *vlabpos = 0;

  htree->SetBranchAddress("det",&vdet);
  htree->SetBranchAddress("seg",&vseg);
  htree->SetBranchAddress("hcid",&vhcid);
  htree->SetBranchAddress("depE",&vdepE);
  htree->SetBranchAddress("labpos",&vlabpos);

  int Nevts = htree->GetEntriesFast();
  int ievt, hid=0;
  for(ievt=0; ievt<Nevts; ievt++){
    if(ievt%1000==0) cout<<"\r load "<<ievt<<" / "<<Nevts<<" EventHits..."<<flush;
    htree->GetEntry(ievt);

    for(int i=0; i<vdet->size(); i++){
      det = vdet->at(i);
      seg = vseg->at(i);

      float depE = vdepE->at(i); // keV
      TMatrixD LabPos(3,1);
      for(int it=0; it<3; it++) LabPos(it,0) = vlabpos->at(i)[it];
      TMatrixD DetPos= Lab2DetPos(det,LabPos);
      
      Hit ahit;
      ahit.det = det;
      ahit.seg = seg;
      ahit.depE = depE;

      for(int it=0; it<3; it++){
	ahit.labpos[it] = LabPos(it,0);
	ahit.detpos[it] = DetPos(it,0);
      }
      fHits.push_back(ahit);

      // connect with HCs
      for(int ii=0; ii<vhcid->at(i).size(); ii++){
	ihc = vhcid->at(i)[ii];
	fAllHCs[ihc].hid.push_back(hid);
	fAllHCs[ihc].nhits++;
      }
      hid++;
    }
  }
  cout<<"\r load "<<ievt<<" / "<<Nevts<<" EventHits..."<<endl;


  TH1D *hnhits = new TH1D("h","nhits",50,-0.5,49.5);
  TH2D *hxy0 = new TH2D("hxy0","HCs Y:X {nhits>3}",200,-50,50,200,-50,50);
  TH2D *hxz0 = new TH2D("hxz0","HCs Z:X {nhits>3}",200,-50,50,200,-5,95);
  TH2D *hxy = new TH2D("hxy","hxy",100,-50,50,100,-50,50);
  TH2D *hxz = new TH2D("hxz","hxz",100,-50,50,100,-5,95);

  int dihc = -Number;
  for(ihc=0; ihc<fAllHCs.size(); ihc++){
    if(fAllHCs[ihc].det!=detid) continue;
    hnhits->Fill(fAllHCs[ihc].nhits);
    if(fAllHCs[ihc].nhits<4) continue;

    double x = fAllHCs[ihc].detpos[0];
    double y = fAllHCs[ihc].detpos[1];
    double z = fAllHCs[ihc].detpos[2];

    hxy0->Fill(x,y);
    hxz0->Fill(x,z);

    if(dihc<0){
      if(fAllHCs[ihc].seg!=segid) continue;
      if(fAllHCs[ihc].nhits<nhitslimit) continue;
      if(sqrt(x*x+y*y)<25) continue;
      if(sqrt(x*x+y*y)>35) continue;
      dihc++;
      if(dihc==0) dihc = ihc;
    }
  }
    

  for(int i=0; i<fAllHCs[dihc].hid.size();i++){
    hid = fAllHCs[dihc].hid[i];
    double x = fHits[hid].detpos[0];
    double y = fHits[hid].detpos[1];
    double z = fHits[hid].detpos[2];

    hxy->Fill(x,y);
    hxz->Fill(x,z);
  }

  hxy->SetTitle(Form("Hits Y:X {Det%d, Seg%d, nhits-%d}",fAllHCs[dihc].det,fAllHCs[dihc].seg,fAllHCs[dihc].nhits));
  hxz->SetTitle(Form("Hits Z:X {Det%d, Seg%d, nhits-%d}",fAllHCs[dihc].det,fAllHCs[dihc].seg,fAllHCs[dihc].nhits));
  
  cout<<"Draw det"<<fAllHCs[dihc].det<<" seg"<<fAllHCs[dihc].seg<<" nhits-"<<fAllHCs[dihc].nhits<<endl;

  //hnhits->Draw();
  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(2,2);
  c->cd(1);
  hxy0->Draw("colz");
  c->cd(2);
  hxy->Draw("colz");
  c->cd(3);
  hxz0->Draw("colz");
  c->cd(4);
  hxz->Draw("colz");

}
