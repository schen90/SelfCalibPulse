#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "stdlib.h"
#include <vector>

using namespace std;

#define GridDist 2. // 2mm grid
#define GridMaxSteps 50

void DrawGridPos(){
  string gridfile = "pulseA.root";
  double GridRange[3][2];
  GridRange[0][0] = -40.250;    GridRange[0][1] = 39.750;
  GridRange[1][0] = -40.250;    GridRange[1][1] = 39.750;
  GridRange[2][0] = 2.250;      GridRange[2][1] = 90.250;

  TFile *fgrid = new TFile(gridfile.c_str());
  if(!fgrid->IsOpen()){
    cerr<<"cannot find gridfile "<<gridfile<<endl;
    return;
  }

  int gridimap[GridMaxSteps][GridMaxSteps][GridMaxSteps];
  vector<Int_t> GridSeg;
  for(int ix=0; ix<GridMaxSteps; ix++)
    for(int iy=0; iy<GridMaxSteps; iy++)
      for(int iz=0; iz<GridMaxSteps; iz++)
        gridimap[ix][iy][iz] = -1;

  TTree *gridtree = (TTree *)fgrid->Get("tree");

  Int_t gridsegi;
  Double_t gridposi[3];

  gridtree->SetBranchAddress("seg",&gridsegi);
  gridtree->SetBranchAddress("pos",gridposi);
  int npoint = gridtree->GetEntriesFast();

  int idx[3];
  for(int ipoint=0; ipoint<npoint; ipoint++){
    gridtree->GetEntry(ipoint);

    for(int i=0; i<3; i++){
      idx[i] = (int)((gridposi[i]-GridRange[i][0]) / GridDist + 0.5);
      if(idx[i]<0 || idx[i]>=GridMaxSteps){
        cerr<<"grid point outside Map range!!!"<<endl;
        return;
      }
    }

    GridSeg.push_back(gridsegi-1); //seg in db start from 1...
    gridimap[idx[0]][idx[1]][idx[2]] = ipoint;
  }
  cout<<"load "<<npoint<<" points from "<<gridfile<<endl;

  TFile *fout = new TFile("draw/CheckGrid.root","RECREATE");
  TTree *tree = new TTree("tree","check grid");
  Float_t pos[3];
  int seg;
  tree->Branch("pos",pos,"pos[3]/F");
  tree->Branch("seg",&seg,"seg/I");

  int N = 10000000;
  for(int i=0; i<N; i++){
    if(i%1000==0) cout<<"\r finish "<<i<<" / "<<N<<flush;
    for(int ix=0; ix<3; ix++) pos[ix] = gRandom->Uniform(GridRange[ix][0]-0.9,GridRange[ix][1]+0.9);

    int idx[3];
    for(int ix=0; ix<3; ix++){
      idx[ix] = (int)((pos[ix]-GridRange[ix][0]) / GridDist + 0.5);
      if(idx[ix]<0 || idx[ix]>GridMaxSteps-1) return;
    }

    int ipoint = gridimap[idx[0]][idx[1]][idx[2]];
    if(ipoint<0) seg = -1;
    else         seg = GridSeg[ipoint];
    
    tree->Fill();
  }

  fout->cd();
  tree->Write();
  fout->Close();
  return;
}
