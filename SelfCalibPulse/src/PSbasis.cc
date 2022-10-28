#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TRandom.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

#include "PSbasis.hh"

using namespace std;

PSbasis::PSbasis(){
  ReadPSbasis();
}

PSbasis::~PSbasis() {
}


void PSbasis::ReadPSbasis(){

  cout<<"\e[1;31m Read basis for linear interpolation PS: \e[0m"<<endl;
  string dbfile[3] = {"G4Sim/pulsedb/LibTrap_A001.root",
                      "G4Sim/pulsedb/LibTrap_B001.root",
                      "G4Sim/pulsedb/LibTrap_C001.root"};

  for(int itype=0; itype<NType; itype++){

    for(int ix=0; ix<3; ix++){
      range[itype][ix][0] = 1000;
      range[itype][ix][1] = -1000;
    }
    for(int ix=0; ix<GridMaxSteps; ix++)
      for(int iy=0; iy<GridMaxSteps; iy++)
        for(int iz=0; iz<GridMaxSteps; iz++)
          imap[itype][ix][iy][iz] = -1;

    TFile *fdb = new TFile(dbfile[itype].c_str());
    if(!fdb->IsOpen()){
      cerr<<"cannot find dbfile "<<dbfile[itype]<<endl;
      return;
    }

    TTree *dbtree = (TTree *)fdb->Get("tree");

    Int_t dbsegi;
    Float_t dbposi[3];
    Float_t dbspulsei[NSig*NSegCore];

    dbtree->SetBranchAddress("seg",&dbsegi);
    dbtree->SetBranchAddress("pos",dbposi);
    dbtree->SetBranchAddress("spulse",dbspulsei);
    int npoint = dbtree->GetEntriesFast();

    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);
      for(int ix=0; ix<3; ix++){
        if(dbposi[ix]<range[itype][ix][0]) range[itype][ix][0] = dbposi[ix];
        if(dbposi[ix]>range[itype][ix][1]) range[itype][ix][1] = dbposi[ix];
      }
    }

    TMatrixD tmppos(3,1);
    int idx[3];
    TMatrixD tmpspulse(NSig*NSegCore,1);
    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);

      for(int i=0; i<3; i++){
        tmppos(i,0)=dbposi[i];
        idx[i] = (int)((dbposi[i]-range[itype][i][0]) / GridDist + 0.5);
        if(idx[i]<0 || idx[i]>=GridMaxSteps){
          cerr<<"grid point outside Map range!!!"<<endl;
          return;
        }
      }

      for(int iseg=0; iseg<NSeg; iseg++){
        for(int i=0; i<NSig; i++)
          tmpspulse(iseg*NSig+i,0)=dbspulsei[iseg*NSig+i];
      }
      for(int i=0; i<NSig; i++) tmpspulse(NSeg*NSig+i,0)=dbspulsei[NSeg*NSig+i];

      dbseg[itype].push_back(dbsegi); // start from 0
      dbpos[itype].push_back(tmppos);
      imap[itype][idx[0]][idx[1]][idx[2]] = ipoint;

      dbspulse[itype].push_back(tmpspulse);
    }

    cout<<"load "<<npoint<<" points from "<<dbfile[itype]<<endl;
    fdb->Close();
  }

  return;
}


Int_t PSbasis::GetPS(int itype, TMatrixD pos, double energy, int &seg, TMatrixD &spulse) {

  int ngrid = 0;

  int idx[3];
  for(int ix=0; ix<3; ix++){
    idx[ix] = (int)((pos(ix,0)-range[itype][ix][0]) / GridDist + 0.5);
  }

  // find grid within 2mm cube
  vector <int> ip; // id list of close grid point for interpolation
  vector <int> ip2; // id list of close grid point 4mm for extrapolation
  vector <double> ipdist;
  vector <double> ip2dist;
  double mindist = 5;
  int minip = -1;
  int tmpseg = -1;
  int npoint = dbpos[itype].size();

  int idxrange = 2;
  for(int ix=idx[0]-idxrange; ix<=idx[0]+idxrange; ix++)
    for(int iy=idx[1]-idxrange; iy<=idx[1]+idxrange; iy++)
      for(int iz=idx[2]-idxrange; iz<=idx[2]+idxrange; iz++){
        if(ix<0 || ix>=GridMaxSteps) continue;
        if(iy<0 || iy>=GridMaxSteps) continue;
        if(iz<0 || iz>=GridMaxSteps) continue;

        int ipoint = imap[itype][ix][iy][iz];
        if(ipoint<0) continue;

	if(kextrapol){
          if(fabs(pos(0,0)-dbpos[itype][ipoint](0,0))>4) continue;
          if(fabs(pos(1,0)-dbpos[itype][ipoint](1,0))>4) continue;
          if(fabs(pos(2,0)-dbpos[itype][ipoint](2,0))>4) continue;
          double disttmp = sqrt((pos-dbpos[itype][ipoint]).Sqr().Sum());
          ip2.push_back(ipoint);
          ip2dist.push_back(disttmp);
        }

	if(fabs(pos(0,0)-dbpos[itype][ipoint](0,0))>2) continue;
        if(fabs(pos(1,0)-dbpos[itype][ipoint](1,0))>2) continue;
        if(fabs(pos(2,0)-dbpos[itype][ipoint](2,0))>2) continue;
        double disttmp = sqrt((pos-dbpos[itype][ipoint]).Sqr().Sum());
        ip.push_back(ipoint);
        ipdist.push_back(disttmp);
        if(disttmp<mindist){
          mindist = disttmp;
          minip = ipoint;
          tmpseg = dbseg[itype][ipoint]; // seg of closest grid point
        }

        if(kextrapol){ ip2.pop_back(); ip2dist.pop_back();}
      }

  if(tmpseg<0){ // remove evt if interaction point outside grid map
    ngrid = 0;
    return ngrid;
  }
  seg = tmpseg;

  // remove grid from different segment
  for(int i=0; i<ip.size();){
    if(dbseg[itype][ip[i]]!=tmpseg){
      ip.erase(ip.begin()+i);
      ipdist.erase(ipdist.begin()+i);
    }else{
      i++;
    }
  }
  // sort according to dist
  for(int i=0; i<ip.size(); i++){
    for(int j=i+1; j<ip.size(); j++){
      if(ipdist[i] > ipdist[j]){
        swap(ipdist, i, j);
        swap(ip, i, j);
      }
    }
  }

  // extrapolation ---------------------------------------------------
  if(kextrapol){
    if(ip.size()<8){
      // remove grid from different segment
      for(int i=0; i<ip2.size();){
        if(dbseg[itype][ip2[i]]!=tmpseg){
          ip2.erase(ip2.begin()+i);
          ip2dist.erase(ip2dist.begin()+i);
        }else{
          i++;
        }
      }
      // sort according to dist
      for(int i=0; i<ip2.size(); i++){
        for(int j=i+1; j<ip2.size(); j++){
          if(ip2dist[i] > ip2dist[j]){
            swap(ip2dist, i, j);
            swap(ip2, i, j);
          }
        }
      }
      // add extrapolation points
      for(int i=0; i<ip2.size();i++){
        ip.push_back(ip2[i]);
      }
    }
  }
  // -----------------------------------------------------------------

  if(mindist==0){

    ngrid = 1;
    spulse = energy*dbspulse[itype][minip];

  }else{// more than 1 grid

    // prepare spulse list
    vector <TMatrixD> spulselist;
    vector <TMatrixD> poslist;
    vector <vector<int>> igridlist;
    for(int i=0; i<ip.size(); i++){
      spulselist.push_back(dbspulse[itype][ip[i]]);
      poslist.push_back(dbpos[itype][ip[i]]);
      igridlist.push_back(vector<int>(1));
      igridlist[i][0] = i;
    }

    // linear combine pulse
    for(int iaxis=0; iaxis<3; iaxis++){ //combine iaxis direction
      int iaxis1 = (iaxis+1)%3;
      int iaxis2 = (iaxis+2)%3;

      for(int i=0; i<poslist.size(); ){
        for(int j=i+1; j<poslist.size(); ){
          bool kcombine = true;
          if(fabs(poslist[i](iaxis1,0)-poslist[j](iaxis1,0))>0.1) kcombine = false;
          if(fabs(poslist[i](iaxis2,0)-poslist[j](iaxis2,0))>0.1) kcombine = false;

          if(kcombine){
            double factori = (pos(iaxis,0)-poslist[j](iaxis,0)) / (poslist[i](iaxis,0)-poslist[j](iaxis,0));
            double factorj = (pos(iaxis,0)-poslist[i](iaxis,0)) / (poslist[j](iaxis,0)-poslist[i](iaxis,0));
            spulselist[i] = factori*spulselist[i] + factorj*spulselist[j];
            poslist[i](iaxis,0) = pos(iaxis,0);

            if(factori==0){
              igridlist[i].clear();
            }

            if(factorj!=0)
              for(int ig=0; ig<igridlist[j].size(); ig++){
                igridlist[i].push_back(igridlist[j][ig]);
              }

            spulselist.erase(spulselist.begin()+j);
            poslist.erase(poslist.begin()+j);
            igridlist.erase(igridlist.begin()+j);

          }else{
            j++;
          }
        }// end of loop j
	
        if(!kextrapol) poslist[i](iaxis,0) = pos(iaxis,0); // reduced interpolation
        i++;
      }// end of loop i
    }// end of loop iaxis

    if(poslist.size()>1){

      if(!kextrapol){
        cerr<<"something wrong..."<<endl;

      }else{ // if use extrapolation
	
        while(poslist.size()>0){
          if(fabs(pos(0,0)-poslist[0](0,0))<0.01 &&
             fabs(pos(1,0)-poslist[0](1,0))<0.01 &&
             fabs(pos(2,0)-poslist[0](2,0))<0.01)
            break;

          spulselist.erase(spulselist.begin());
          poslist.erase(poslist.begin());
          igridlist.erase(igridlist.begin());
        }

        if(poslist.size()==0){ // cannot get simulated position in XYZ combine
          ngrid = -1;
          return ngrid;
        }// end of cannot get simulated position in XYZ combine
	
      }// end of kextrapol
      
    }// end of poslist>1

    ngrid = igridlist[0].size();
    spulse = energy*spulselist[0];
    
  } // end of more than 1 grid

  return ngrid;
}
