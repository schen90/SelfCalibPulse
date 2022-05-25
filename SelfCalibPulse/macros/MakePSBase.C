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

const int nsig = 121;
const int nseg = 36;
const int nsegcore = 37;

double griddist  = 2; // 2mm grid
double maxdist = 2.0; // max dist in one axis to assign a selfcalib point to grid
double range[3][2];
const int MaxSteps = 50;
int imap[MaxSteps][MaxSteps][MaxSteps];

// pulsedb
vector<int>      dbseg;
vector<TMatrixD> dbpos;

// selfcalib
vector<vector<int>> calibpoint;
vector<TMatrixD> calibpos;
vector<TMatrixD> calibpulse;

void swap(vector<TMatrixD> &v, int m, int l){
  TMatrixD temp;
  temp.ResizeTo( v[m] );
  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
}

void swap(vector<double> &v, int m, int l){
  double temp;
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

void LoadDBPos(string DBPosfile);
TMatrixD GetPulse(int igrid, double *cpos);

// main
void MakePSBase(string DBPosfile, string PSCfile, string Basefile){

  // load dbpos
  LoadDBPos(DBPosfile);

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
    seg = seg+1; // seg in db start from 1...
    
    if(npaths<=30) continue;

    TMatrixD tmppos(3,1);
    int idx[3];
    for(int ix=0; ix<3; ix++){
      tmppos(ix,0)=cadpos2[ix];
      idx[ix] = (int)((cadpos2[ix]-range[ix][0]) / griddist + 0.5);
    }
    TMatrixD tmppulse(nsegcore*nsig,1);
    for(int iseg=0; iseg<nsegcore; iseg++)
      for(int isig=0; isig<nsig; isig++)
	tmppulse(iseg*nsig+isig,0) = spulse[iseg][isig];

    calibpos.push_back(tmppos);
    calibpulse.push_back(tmppulse);

    // find grid pos in range
    int idxrange = 1;
    for(int ix=idx[0]-idxrange; ix<=idx[0]+idxrange; ix++)
      for(int iy=idx[1]-idxrange; iy<=idx[1]+idxrange; iy++)
	for(int iz=idx[2]-idxrange; iz<=idx[2]+idxrange; iz++){
	  if(ix<0 || ix>=MaxSteps) continue;
	  if(iy<0 || iy>=MaxSteps) continue;
	  if(iz<0 || iz>=MaxSteps) continue;

	  int ipoint = imap[ix][iy][iz];
	  if(ipoint<0) continue;

	  if(dbseg[ipoint]!=seg) continue;
	  if(fabs(tmppos(0,0)-dbpos[ipoint](0,0))>maxdist) continue;
	  if(fabs(tmppos(1,0)-dbpos[ipoint](1,0))>maxdist) continue;
	  if(fabs(tmppos(2,0)-dbpos[ipoint](2,0))>maxdist) continue;

	  calibpoint[ipoint].push_back(calibpos.size()-1);
	}
  }
  cout<<"\r finish "<<ientry<<" / "<<nentries<<" points from PSC"<<endl;

  // output
  double pos[3];
  double core[121];
  double segpulse[4356];
  int npoints;
  double cpos[3];
  double dist;
  
  TFile *fout = new TFile(Basefile.c_str(),"RECREATE");
  TTree *outtree = new TTree("tree","SelfCalib Grid Base");

  outtree->Branch("det",&det);
  outtree->Branch("seg",&seg);
  outtree->Branch("pos",pos,"pos[3]/D");
  outtree->Branch("core",core,"core[121]/D");
  outtree->Branch("spulse",segpulse,"spulse[4356]/D");
  outtree->Branch("npoints",&npoints);
  outtree->Branch("cpos",cpos,"cpos[3]/D");
  outtree->Branch("dist",&dist);

  // combine pulse shape at grid position
  int ngrid = dbpos.size();
  int igrid;
  for(igrid=0; igrid<ngrid; igrid++){
    if(igrid%1000==0) cout<<"\r finish "<<igrid<<" / "<<ngrid<<" grid points"<<flush;
    seg = dbseg[igrid];

    for(int ix=0; ix<3; ix++){
      pos[ix] = dbpos[igrid](ix,0);
      cpos[ix] = 0;
    }
    for(int isig=0; isig<121; isig++) core[isig] = 0;
    for(int isig=0; isig<4356; isig++) segpulse[isig] = 0;
    dist = -1;
    npoints = calibpoint[igrid].size();

    if(npoints>0){
      TMatrixD tmppulse = GetPulse(igrid,cpos);
      dist = 0;
      for(int ix=0; ix<3; ix++) dist += pow(cpos[ix]-pos[ix],2);
      dist = sqrt(dist);
      
      for(int isig=0; isig<121; isig++) core[isig] = tmppulse(4356+isig,0);
      for(int isig=0; isig<4356; isig++) segpulse[isig] = tmppulse(isig,0);
    }
    
    outtree->Fill();
  }
  cout<<"\r finish "<<igrid<<" / "<<ngrid<<" grid points"<<endl;
    
  fout->cd();
  outtree->Write();
  fout->Close();

  return;
}


void LoadDBPos(string DBPosfile){
  range[0][0] = -40.250;    range[0][1] = 39.750;
  range[1][0] = -40.250;    range[1][1] = 39.750;
  range[2][0] = 2.250;      range[2][1] = 90.250;
  for(int ix=0; ix<MaxSteps; ix++)
    for(int iy=0; iy<MaxSteps; iy++)
      for(int iz=0; iz<MaxSteps; iz++)
	imap[ix][iy][iz] = -1;

  TFile *fdb = new TFile(DBPosfile.c_str());
  if(!fdb->IsOpen()){
    cerr<<"cannot find dbpos "<<DBPosfile<<endl;
    return;
  }

  TTree *dbtree = (TTree *)fdb->Get("tree");

  Int_t dbsegi;
  Double_t dbposi[3];

  dbtree->SetBranchAddress("seg",&dbsegi);
  dbtree->SetBranchAddress("pos",dbposi);

  int npoint = dbtree->GetEntriesFast();

  TMatrixD tmppos(3,1);
  int idx[3];
  for(int ipoint=0; ipoint<npoint; ipoint++){
    dbtree->GetEntry(ipoint);

    for(int i=0; i<3; i++){
      tmppos(i,0)=dbposi[i];
      idx[i] = (int)((dbposi[i]-range[i][0]) / griddist + 0.5);
      if(idx[i]<0 || idx[i]>=MaxSteps){
	cerr<<"grid point outside Map range!!!"<<endl;
	return;
      }
    }

    dbseg.push_back(dbsegi);
    dbpos.push_back(tmppos);
    imap[idx[0]][idx[1]][idx[2]] = ipoint;

    calibpoint.push_back(vector<int>());
  }

  cout<<"load "<<npoint<<" points from "<<DBPosfile<<endl;
  fdb->Close();
  
  return;
}


TMatrixD GetPulse(int igrid, double *cpos){
  TMatrixD tmppulse(nsegcore*nsig,1);
  for(int iseg=0; iseg<nsegcore; iseg++)
    for(int isig=0; isig<nsig; isig++)
      tmppulse(iseg*nsig+isig,0)=0;

  if(calibpoint[igrid].size()==0) return tmppulse;

  double pos[3];
  for(int ix=0; ix<3; ix++) pos[ix]=dbpos[igrid](ix,0);

  vector<TMatrixD> poslist;
  vector<TMatrixD> pulselist;
  vector<double> dist;
  for(int ip=0; ip<calibpoint[igrid].size(); ip++){
    int ipoint = calibpoint[igrid][ip];
    poslist.push_back(calibpos[ipoint]);
    pulselist.push_back(calibpulse[ipoint]);

    double tmpdist = sqrt((calibpos[ipoint]-dbpos[igrid]).Sqr().Sum());
    dist.push_back(tmpdist);
  }

  // sort according to dist, small to large
  for(int i=0; i<poslist.size(); i++)
    for(int j=i+1; j<poslist.size(); j++)
      if(dist[i]>dist[j]){
	swap(poslist,i,j);
	swap(pulselist,i,j);
	swap(dist,i,j);
      }

  //***************************************************//
  // use the closest point
  //for(int ix=0; ix<3; ix++) cpos[ix] = poslist[0](ix,0);
  //tmppulse = pulselist[0];
  //return tmppulse;
  //***************************************************//

  // interpolation
  for(int repeat=0; repeat<3; repeat++){

    int nsize = poslist.size();
    for(int i=0; i<nsize; i++){
      for(int j=i+1; j<nsize; j++){
	for(int ix=0; ix<3; ix++){
	  if(fabs(poslist[i](ix,0)-pos[ix])<0.1 || fabs(poslist[j](ix,0)-pos[ix])<0.1) continue;
	  double factori = (pos[ix]-poslist[j](ix,0)) / (poslist[i](ix,0)-poslist[j](ix,0));
	  double factorj = (pos[ix]-poslist[i](ix,0)) / (poslist[j](ix,0)-poslist[i](ix,0));

	  TMatrixD tmppos(3,1);
	  tmppos = factori*poslist[i] + factorj*poslist[j];

	  // too far away
	  double tmpdist = sqrt((tmppos-dbpos[igrid]).Sqr().Sum());
	  if(fabs(tmppos(0,0)-dbpos[igrid](0,0))>maxdist) continue;
	  if(fabs(tmppos(1,0)-dbpos[igrid](1,0))>maxdist) continue;
	  if(fabs(tmppos(2,0)-dbpos[igrid](2,0))>maxdist) continue;
	  if(tmpdist>dist[i] || tmpdist>dist[j]) continue;

	  // already exist
	  int skip = 0;
	  for(int k=0; k<poslist.size(); k++){
	    if(sqrt((tmppos-poslist[k]).Sqr().Sum())<0.1){
	      skip = 1;
	      break;
	    }
	  }
	  if(skip==1) continue;

	  // add to list
	  tmppulse = factori*pulselist[i] + factorj*pulselist[j];
	  poslist.push_back(tmppos);
	  pulselist.push_back(tmppulse);
	  dist.push_back(tmpdist);
	}
      }
    }
    //cout<<"repeat "<<repeat<<" poslist.size() = "<<poslist.size()<<endl;

    // sort according to dist, small to large
    for(int i=0; i<poslist.size(); i++)
      for(int j=i+1; j<poslist.size(); j++)
	if(dist[i]>dist[j]){
	  swap(poslist,i,j);
	  swap(pulselist,i,j);
	  swap(dist,i,j);
	}

    if(dist[0]<0.1) break;
  }

  for(int ix=0; ix<3; ix++) cpos[ix] = poslist[0](ix,0);
  tmppulse = pulselist[0];
  return tmppulse;
}


#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc>3){
    MakePSBase(string(argv[1]), string(argv[2]), string(argv[3]));
  }else{
    MakePSBase("G4Sim/pulsedb/pulseA.root", "PSCfiles/runMix/run31/it1/Det0000_fit4.root", "PSBase/Det0000.root");
  }
  return 0;
}
#endif
