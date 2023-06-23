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
  string gridfile[3] = {"pulsedb/LibTrap_A001.root",
			"pulsedb/LibTrap_B001.root",
			"pulsedb/LibTrap_C001.root"};
  for(int itype=0; itype<NType; itype++){

    LoadGrid(itype, gridfile[itype]);
  }

  NDets = 0;
  LoadMatrix("LookUp/CrystalPositionLookUpTable");
  MakeSegmentMap();
  MakeSegPos();

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

  for(int iseg=0; iseg<NSEGS; iseg++){
    NSegGrid[itype][iseg] = 0;
    LocalSegPos[itype][iseg].ResizeTo(3,1); LocalSegPos[itype][iseg].Zero();
  }
  
  TTree *gridtree = (TTree *)fgrid->Get("tree");

  Int_t gridsegi;
  Float_t gridposi[3];

  gridtree->SetBranchAddress("seg",&gridsegi);
  gridtree->SetBranchAddress("pos",gridposi);
  int npoint = gridtree->GetEntriesFast();

  // find grid range
  GridRange[itype][0][0] = GridRange[itype][1][0] = GridRange[itype][2][0] = 999;
  GridRange[itype][0][1] = GridRange[itype][1][1] = GridRange[itype][2][1] = -999;
  for(int ipoint=0; ipoint<npoint; ipoint++){
    gridtree->GetEntry(ipoint);
    if(gridposi[0]<GridRange[itype][0][0]) GridRange[itype][0][0] = gridposi[0];
    if(gridposi[0]>GridRange[itype][0][1]) GridRange[itype][0][1] = gridposi[0];

    if(gridposi[1]<GridRange[itype][1][0]) GridRange[itype][1][0] = gridposi[1];
    if(gridposi[1]>GridRange[itype][1][1]) GridRange[itype][1][1] = gridposi[1];

    if(gridposi[2]<GridRange[itype][2][0]) GridRange[itype][2][0] = gridposi[2];
    if(gridposi[2]>GridRange[itype][2][1]) GridRange[itype][2][1] = gridposi[2];
  }
  cout<<"type "<<itype<<" : "
      <<Form("x %.3f ~ %.3f ; ",GridRange[itype][0][0],GridRange[itype][0][1])
      <<Form("y %.3f ~ %.3f ; ",GridRange[itype][1][0],GridRange[itype][1][1])
      <<Form("z %.3f ~ %.3f",GridRange[itype][2][0],GridRange[itype][2][1])<<endl;

  
  TMatrixD tmppos(3,1);
  int idx[3];
  for(int ipoint=0; ipoint<npoint; ipoint++){
    gridtree->GetEntry(ipoint);

    for(int i=0; i<3; i++){
      tmppos(i,0)=gridposi[i];
      idx[i] = (int)((gridposi[i]-GridRange[itype][i][0]) / fstep + 0.5);
      if(idx[i]<0 || idx[i]>=GridMaxSteps){
	cerr<<"grid point outside Map range!!!"<<endl;
	return;
      }
    }

    GridSeg[itype].push_back(gridsegi); //seg in db start from 0...
    gridimap[itype][idx[0]][idx[1]][idx[2]] = ipoint;

    LocalSegPos[itype][gridsegi] = LocalSegPos[itype][gridsegi] + tmppos;
    NSegGrid[itype][gridsegi]++;
  }

  cout<<"load "<<npoint<<" points from "<<gridfile<<" for type "<<itype<<endl;

  for(int iseg=0; iseg<NSEGS; iseg++){
    if(NSegGrid[itype][iseg]>0)
      LocalSegPos[itype][iseg] = 1./NSegGrid[itype][iseg] * LocalSegPos[itype][iseg];
  }

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


void AGATAgeo::MakeSegPos(){
  cout<<"Make Segment Center Position..";

  int detid, segid;
  for(int detid=0; detid<NDets; detid++){
    int itype = detid%3;
    for(int iseg=0; iseg<NSEGS; iseg++){
      SegPos[detid][iseg].ResizeTo(3,1);
      SegPos[detid][iseg] = Det2LabPos(detid,LocalSegPos[itype][iseg]);
    }
  }
  
  cout<<endl;
  return;
}


void AGATAgeo::MakeSegmentMap(){
  cout<<"Make Segment Map..";
  //cout<<"\e[1;32m assign weight for comparison hitseg-1, core-1, nxsec-1, nxsli-1 \e[0m"<<endl;
  for(int iseg=0; iseg<NSEGS; iseg++){
    NextSec[iseg][0] = NextSec[iseg][1] = -1;
    NextSli[iseg][0] = NextSli[iseg][1] = -1;
    int NextSli2 = -1;
    
    int isec = iseg/NSLIC;
    int isli = iseg%NSLIC;

    for(int jseg=0; jseg<NSEGS; jseg++){

      if(jseg==iseg) continue;

      int jsec = jseg/NSLIC;
      int jsli = jseg%NSLIC;

      int distV = abs(jsli - isli);
      int distH = abs(jsec - isec);
      distH = min(distH, abs(jsec - isec + NSECT));
      distH = min(distH, abs(jsec - isec - NSECT));

      if(distV==0 && distH==1){
	if(NextSec[iseg][0]<0) NextSec[iseg][0] = jseg;
	else                   NextSec[iseg][1] = jseg;
      }

      if(distH==0 && distV==1){
	if(NextSli[iseg][0]<0) NextSli[iseg][0] = jseg;
	else                   NextSli[iseg][1] = jseg;
      }
      if(distH==0 && distV==2) NextSli2 = jseg;
    }

    if(NextSli[iseg][1]<0) NextSli[iseg][1] = NextSli2;
  }

  cout<<endl;

  cout<<"assign segment weight..";
  // assign weight for comparison
  //cout<<"\e[1;32m assign weight for comparison hitseg-1, core-1, nxsec-1, nxsli-1 \e[0m"<<endl;
  for(int iseg=0; iseg<NSEGS; iseg++){
    /*
    for(int iiseg=0; iiseg<NSEGS; iiseg++)
      SegWeight[iseg][iiseg] = 0;

    SegWeight[iseg][iseg] = 1;
    SegWeight[iseg][INDCC] = 1;
    SegWeight[iseg][NextSec[iseg][0]] = 1;    SegWeight[iseg][NextSec[iseg][1]] = 1;
    SegWeight[iseg][NextSli[iseg][0]] = 1;    SegWeight[iseg][NextSli[iseg][1]] = 1;
    */

    for(int iiseg=0; iiseg<NCHAN; iiseg++)
      SegWeight[iseg][iiseg] = 1;

    //SegWeight[iseg][iseg] = 0; // remove direct hit seg
    //SegWeight[iseg][INDCC] = 0; // remove core
  }
  cout<<endl;

  return;
}


void AGATAgeo::GetNextSegs(Int_t iseg, Int_t *fseg){
  fseg[0] = iseg;
  fseg[1] = INDCC;
  fseg[2] = NextSec[iseg][0];
  fseg[3] = NextSec[iseg][1];
  fseg[4] = NextSli[iseg][0];
  fseg[5] = NextSli[iseg][1];

  if( NCOMP == 37){
    int idx = 6;
    int uflg[NCOMP]; for(int i=0; i<NCOMP; i++) uflg[i]=0;
    for(int i=0; i<idx; i++) uflg[fseg[i]]=1;
    for(int i=0; i<NCHAN && idx<NCOMP; i++){
      if(uflg[i]==1) continue;
      fseg[idx] = i;
      idx++;
      uflg[i]=1;
    }
  }
  
  return;
}


bool AGATAgeo::CheckBounds(int detid, int segid, const double *lpos){
  int itype = detid%3;
  TMatrixD LPos(3,1);
  for(int ix=0; ix<3; ix++) LPos(ix,0) = lpos[ix];
  TMatrixD DPos = Lab2DetPos(detid, LPos);

  if(DPos(2,0)<0 || DPos(2,0)>90) return false;

  int idx[3];
  for(int ix=0; ix<3; ix++){
    idx[ix] = (int)((DPos(ix,0)-GridRange[itype][ix][0]) / fstep + 0.5);
    if(idx[ix]<0 || idx[ix]>GridMaxSteps-1) return false;
  }

  int ipoint = gridimap[itype][idx[0]][idx[1]][idx[2]];
  if(ipoint<0) return false;

  if(segid!=GridSeg[itype][ipoint]) return false;

  return true;
}


void AGATAgeo::GetChi2sLimit(int detid, const double *dpos, float chi2slimit[]){
  int itype = detid%3;

  int idx[3];
  for(int ix=0; ix<3; ix++){
    idx[ix] = (int)((dpos[ix]-GridRange[itype][ix][0]) / fstep + 0.5);
    if(idx[ix]<0 || idx[ix]>GridMaxSteps-1) return;
  }

  for(int ix=0; ix<3; ix++) chi2slimit[ix] = gridchi2smap[itype][idx[0]][idx[1]][idx[2]][ix];

  return;
}



#endif
