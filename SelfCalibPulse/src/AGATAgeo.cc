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
  string gridfile[3] = {"G4Sim/pulsedb/LibTrap_A001.root",
			"G4Sim/pulsedb/LibTrap_B001.root",
			"G4Sim/pulsedb/LibTrap_C001.root"};
  for(int itype=0; itype<NType; itype++){
    LoadGrid(itype, gridfile[itype]);
  }

  NDets = 0;
  LoadMatrix("LookUp/CrystalPositionLookUpTable");
  MakeSegPos();
  MakeSegMap(2);

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

  for(int iseg=0; iseg<NSeg; iseg++){
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
      idx[i] = (int)((gridposi[i]-GridRange[itype][i][0]) / GridDist + 0.5);
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

  for(int iseg=0; iseg<NSeg; iseg++){
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
    for(int iseg=0; iseg<NSeg; iseg++){
      SegPos[detid][iseg].ResizeTo(3,1);
      SegPos[detid][iseg] = Det2LabPos(detid,LocalSegPos[itype][iseg]);
    }
  }
  
  cout<<endl;
  return;
}


void AGATAgeo::MakeSegMap(int neighbours){
  cout<<"Make Segment Map...";
  memset(hmask, '0', sizeof(hmask));
  // '1' self, '2' neighbour, '9' CC
  for(int iseg=0; iseg<NSeg; iseg++){
    int isec = iseg/NSli;
    int isli = iseg%NSli;

    for(int jseg=0; jseg<NSeg; jseg++){
      int jsec = jseg/NSli;
      int jsli = jseg%NSli;

      int distV = abs(jsli - isli);
      int distH = abs(jsec - isec);
      distH = min(distH, abs(jsec - isec + NSec));
      distH = min(distH, abs(jsec - isec - NSec));
      int mdist = distV + distH;
      if(mdist<=neighbours && distH<neighbours && distV<neighbours){
        hmask[iseg][jseg] = (iseg==jseg) ? '1' : '2';
      }
    }

    hmask[iseg][NSegCore-1] = '9';
    hmask[iseg][NSegCore] = 0; // to close each line as a string    
  }
  cout<<endl;
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
