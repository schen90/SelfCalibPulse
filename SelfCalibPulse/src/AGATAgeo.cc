#ifndef AGATAGEO_CC
#define AGATAGEO_CC

#include <fstream>
#include <iostream>
#include <algorithm>
#include <x86intrin.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>

#include "AGATAgeo.hh"

using namespace std;

AGATAgeo::AGATAgeo(){
  string gridfile[3] = {"G4Sim/pulsedb/pulseA.root",
			"G4Sim/pulsedb/pulseB.root",
			"G4Sim/pulsedb/pulseC.root"};
  for(int itype=0; itype<NType; itype++){
    GridRange[itype][0][0] = -40.250;    GridRange[itype][0][1] = 39.750;
    GridRange[itype][1][0] = -40.250;    GridRange[itype][1][1] = 39.750;
    GridRange[itype][2][0] = 2.250;      GridRange[itype][2][1] = 90.250;

    LoadGrid(itype, gridfile[itype]);
  }

  NDets = 0;
  LoadMatrix("LookUp/CrystalPositionLookUpTable");
  LoadNextSegTable("LookUp/NextSegTable");
  LoadSegPos("LookUp/SegPosTable");

}

AGATAgeo::~AGATAgeo(){
}

void AGATAgeo::LoadGrid(Int_t itype, string gridfile){
  if(itype<0 || itype>=NType) return;

  TFile *fgrid = new TFile(gridfile.c_str());
  if(!fgrid->IsOpen()){
    cerr<<"cannot find gridfile "<<gridfile<<endl;
    return;
  }

  for(int ix=0; ix<GridMaxSteps; ix++)
    for(int iy=0; iy<GridMaxSteps; iy++)
      for(int iz=0; iz<GridMaxSteps; iz++)
	gridimap[itype][ix][iy][iz] = -1;
  
  TTree *gridtree = (TTree *)fgrid->Get("tree");

  Int_t gridsegi;
  Double_t gridposi[3];

  gridtree->SetBranchAddress("seg",&gridsegi);
  gridtree->SetBranchAddress("pos",gridposi);
  int npoint = gridtree->GetEntriesFast();

  TMatrixD tmppos(3,1);
  int idx[3];
  for(int ipoint=0; ipoint<npoint; ipoint++){
    gridtree->GetEntry(ipoint);

    for(int i=0; i<3; i++){
      tmppos(i,0)=gridposi[i];
      idx[i] = (int)((gridposi[i]-GridRange[itype][i][0]) / GridDist + 0.5);
      if(idx[i]<0 || idx[i]>=GridMaxSteps){
	cerr<<"grid point outside Map range!!!"<<endl;
	return;
      }
    }

    GridSeg[itype].push_back(gridsegi-1); //seg in db start from 1...
    gridimap[itype][idx[0]][idx[1]][idx[2]] = ipoint;
  }

  cout<<"load "<<npoint<<" points from "<<gridfile<<" for type "<<itype<<endl;
  fgrid->Close();

  return;
}


void AGATAgeo::LoadMatrix(string LookUpTable){
  // find input LookUpTable
  ifstream fin;
  int dummy_i;
  fin.open(LookUpTable.c_str());
  if(!fin){ cerr<<"Cannot find "<<LookUpTable<<endl; return;}
  cout<<"\e[1;32m find CrystalPosition from "<<LookUpTable<<"... \e[0m"<<endl;

  int ir;
  for(int i=0; i<MaxNDets; i++){
    ir = -1;
    fin >> ir >> dummy_i;
    if(ir<0 || ir>=MaxNDets){ ir=i-1; break;}

    Tr[ir].ResizeTo(3,1); Tr[ir].Zero();
    Rt[ir].ResizeTo(3,3); Rt[ir].Zero();
    Rt2[ir].ResizeTo(3,3); Rt2[ir].Zero();
    for(int it=0; it<3; it++) fin >> Tr[ir](it,0);
    for(int it=0; it<3; it++){
      fin >> dummy_i;
      for(int it2=0; it2<3; it2++) fin >> Rt[ir](it,it2);
      for(int it2=0; it2<3; it2++) Rt2[ir](it,it2) =  Rt[ir](it,it2);
    }
    Rt[ir].Invert();  // change to rot from world frame -> detector frame
  }
  fin.close();
  NDets = ir+1;
  cout<<"read position for "<<NDets<<" detectors"<<endl;
  
  return;
}

TMatrixD AGATAgeo::Lab2DetPos(Int_t idet, TMatrixD LabPos){
  if(idet<0 || idet>=NDets){
    cerr<<"cannot find matrix for idet = "<<idet<<endl;
    return LabPos;
  }
  TMatrixD DetPos(3,1);
  DetPos = Rt[idet]*(LabPos-Tr[idet]);
  return DetPos;
}


TMatrixD AGATAgeo::Det2LabPos(Int_t idet, TMatrixD DetPos){
  if(idet<0 || idet>=NDets){
    cerr<<"cannot find matrix for idet = "<<idet<<endl;
    return DetPos;
  }
  TMatrixD LabPos(3,1);
  LabPos = Rt2[idet]*DetPos+Tr[idet];
  return LabPos;
}


void AGATAgeo::LoadSegPos(string SegPosTable){
  // find input SegPosTable
  ifstream fin;
  int dummy_i;
  fin.open(SegPosTable.c_str());
  if(!fin){ cerr<<"Cannot find "<<SegPosTable<<endl; return;}
  cout<<"\e[1;32m find Segment Position from "<<SegPosTable<<"... \e[0m"<<endl;

  int detid, segid, nhits;
  for(int i=0; i<NDets; i++){
    for(int iseg=0; iseg<NSeg; iseg++){
      fin >> detid >> segid >> nhits;
      if(nhits<1){ cout<<"cannot find position for Det "<<detid<<", Seg "<<segid<<endl;}
      SegPos[detid][segid].ResizeTo(3,1);
      for(int it=0; it<3; it++) fin >> SegPos[detid][segid](it,0);
    }
  }

  fin.close();
  cout<<"read position for "<<NDets<<" detectors"<<endl;
  
  return;
}


void AGATAgeo::LoadNextSegTable(string NextSegTable){
  // find input NextSegTable
  ifstream fin;
  int dummy_i;
  fin.open(NextSegTable.c_str());
  if(!fin){ cerr<<"Cannot find "<<NextSegTable<<endl; return;}
  cout<<"\e[1;32m find Segment neighbor from "<<NextSegTable<<"... \e[0m"<<endl;

  int segid, nxsec1, nxsec2, nxsli1, nxsli2;
  int iseg;
  for(iseg=0; iseg<NSeg; iseg++){
    segid = -1;
    fin >> segid >> nxsec1 >> nxsec2 >> nxsli1 >> nxsli2;
    if(segid<0) break;
    NextSec[segid][0] = nxsec1;   NextSec[segid][1] = nxsec2;
    NextSli[segid][0] = nxsli1;   NextSli[segid][1] = nxsli2;
  }

  fin.close();
  cout<<"read segment neighbor for "<<iseg<<" segments"<<endl;
  
  // assign weight for comparison
  //cout<<"\e[1;32m assign weight for comparison hitseg-1, core-1, nxsec-1, nxsli-1 \e[0m"<<endl;
  for(iseg=0; iseg<NSeg; iseg++){
    /*
    for(int iiseg=0; iiseg<NSeg; iiseg++)
      SegWeight[iseg][iiseg] = 0;

    SegWeight[iseg][iseg] = 1;
    SegWeight[iseg][NSegCore-1] = 1;
    SegWeight[iseg][NextSec[iseg][0]] = 1;    SegWeight[iseg][NextSec[iseg][1]] = 1;
    SegWeight[iseg][NextSli[iseg][0]] = 1;    SegWeight[iseg][NextSli[iseg][1]] = 1;
    */

    for(int iiseg=0; iiseg<NSegCore; iiseg++)
      SegWeight[iseg][iiseg] = 1;

    //SegWeight[iseg][iseg] = 0; // remove direct hit seg
    //SegWeight[iseg][NSegCore-1] = 0; // remove core
  }

  return;
}


void AGATAgeo::GetNextSegs(Int_t iseg, Int_t *fseg){
  fseg[0] = iseg;
  fseg[1] = NSegCore-1;
  fseg[2] = NextSec[iseg][0];
  fseg[3] = NextSec[iseg][1];
  fseg[4] = NextSli[iseg][0];
  fseg[5] = NextSli[iseg][1];

#if NSeg_comp == 37
  int idx = 6;
  int uflg[NSeg_comp]; for(int i=0; i<NSeg_comp; i++) uflg[i]=0;
  for(int i=0; i<idx; i++) uflg[fseg[i]]=1;
  for(int i=0; i<NSegCore && idx<NSeg_comp; i++){
    if(uflg[i]==1) continue;
    fseg[idx] = i;
    idx++;
    uflg[i]=1;
  }
#endif
  
  return;
}


bool AGATAgeo::CheckBounds(int detid, int segid, const double *lpos){
  int itype = detid%3;
  TMatrixD LPos(3,1);
  for(int ix=0; ix<3; ix++) LPos(ix,0) = lpos[ix];
  TMatrixD DPos = Lab2DetPos(detid, LPos);

  int idx[3];
  for(int ix=0; ix<3; ix++){
    idx[ix] = (int)((DPos(ix,0)-GridRange[itype][ix][0]) / GridDist + 0.5);
    if(idx[ix]<0 || idx[ix]>GridMaxSteps-1) return false;
  }

  int ipoint = gridimap[itype][idx[0]][idx[1]][idx[2]];
  if(ipoint<0) return false;

  if(segid!=GridSeg[itype][ipoint]) return false;

  return true;
}



#endif