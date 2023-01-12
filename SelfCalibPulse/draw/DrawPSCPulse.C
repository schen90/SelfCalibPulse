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

const int nsig = 60;
const int nsigdb = 56;
const int nseg = 36;
const int nsegcore = 37;

int nskip = 2;

double griddist  = 2; // 2mm grid
double maxdist = 2; // max dist in one axis to assign a selfcalib point to grid
double range[3][2];
const int MaxSteps = 50;
int imap[MaxSteps][MaxSteps][MaxSteps];

// pulsedb
vector<int>      dbseg;
vector<TMatrixD> dbpos;
vector<TMatrixD> dbspulse;


//
TMatrixD realspulse(nsig*nsegcore,1);
TMatrixD calibspulse(nsig*nsegcore,1);


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
void GetDBPulse(float* pos, int &seg, int &ngrid, int &extrpl);


// main
void DrawPSCPulse(string DBPosfile = "pulsedb/LibTrap_A001.root", 
		  string PSCfile   = "PSCfiles/Det0000.root"){

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
  tree->SetBranchAddress("cadpos2",cadpos2);
  tree->SetBranchAddress("spulse",spulse);
  tree->SetBranchAddress("npaths",&npaths);

  int nentries = tree->GetEntriesFast();
  string sout;

  // read PSCfile
  int ientry;
  for(ientry=0; ientry<nentries; ientry++){
    tree->GetEntry(ientry);
    
    if(npaths<=50) continue;

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

    // calc Chi2------------------------------------
    int dbsegi, ngrid, extrpl;
    GetDBPulse(cadpos2, dbsegi, ngrid, extrpl);

    for(int iseg=0; iseg<nsegcore; iseg++){
      for(int isig=0; isig<nsig; isig++){
	calibspulse(iseg*nsig+isig,0) = spulse[iseg][isig];
      }
    }

    double chi2 = 0;
    for(int iseg=0; iseg<nsegcore; iseg++){
      for(int isig=0; isig<nsig; isig++){
	int i = iseg*nsig+isig;
	float chis = calibspulse(i,0) - realspulse(i,0);
	float sigma = realspulse(i,0)>0.01? realspulse(i,0) : 0.01;
	chi2 += chis*chis/sigma;
      }
    }

    if(dbsegi!=seg || ngrid!=8) continue;
    //if(chi2>10) continue;
    // ---------------------------------------------
    if(nskip>0){ nskip--; continue;}

    sout = Form("ientry = %d seg = %d npaths = %d cadpos = %.2f %.2f %.2f Chi2 = %.3f", ientry, seg, npaths, cadpos2[0], cadpos2[1], cadpos2[2], chi2);
    cout<<sout<<endl;

    break;
  }

  Float_t x[2220],y1[2220], y2[2220];  
  Float_t ymin = 0, ymax = 0;
  for(int iseg=0; iseg<nsegcore; iseg++)
    for(int i=0; i<nsig; i++)
      x[iseg*nsig+i] = iseg+i*1./nsig-0.5;

  for(int iseg=0; iseg<nsegcore; iseg++){
    for(int i=0; i<nsig; i++){
      y1[iseg*nsig+i] = calibspulse(iseg*nsig+i,0);
      y2[iseg*nsig+i] = realspulse(iseg*nsig+i,0);

      if(y1[iseg*nsig+i]>ymax) ymax=y1[iseg*nsig+i];
      if(y1[iseg*nsig+i]<ymin) ymin=y1[iseg*nsig+i];
      if(y2[iseg*nsig+i]>ymax) ymax=y2[iseg*nsig+i];
      if(y2[iseg*nsig+i]<ymin) ymin=y2[iseg*nsig+i];
    }
  }

  ymin = 1.1*ymin - 0.001;
  ymax = 1.1*ymax + 0.001;

  TCanvas *c = new TCanvas("c","c",1265,400);
  c->SetMargin(0.06,0.01,0.12,0.01);

  TGraph *gr1 = new TGraph(2220,x,y1);
  gr1->SetTitle("");
  gr1->GetXaxis()->SetTitleOffset(0.9);
  gr1->GetXaxis()->SetTitleSize(0.06);
  gr1->GetXaxis()->SetLabelSize(0.05);
  gr1->GetXaxis()->SetTitle("Segment Number");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->SetTitleOffset(0.5);
  gr1->GetYaxis()->SetTitleSize(0.06);
  gr1->GetYaxis()->SetLabelSize(0.05);
  gr1->GetYaxis()->SetTitle("Normalised Charge");
  gr1->GetYaxis()->CenterTitle();
  gr1->GetXaxis()->SetRangeUser(-0.5,36.5);
  gr1->GetYaxis()->SetRangeUser(ymin,ymax);
  gr1->SetLineColor(2);
  gr1->SetLineWidth(2);

  TGraph *gr2 = new TGraph(2220,x,y2);
  gr2->SetLineColor(1);
  gr2->SetLineWidth(2);

  TLegend *legend = new TLegend(0.42,0.8,0.58,0.95);
  legend->AddEntry(gr1,"SelfCalib","l");
  legend->AddEntry(gr2,"ADL_Basis","l");

  gr1->Draw("APL");
  gr2->Draw("Lsame");
  gr1->Draw("Lsame");
  legend->Draw("same");

  TLatex *tex = new TLatex();
  tex->SetTextAlign(22);
  tex->DrawLatex(18,0.5,sout.c_str());

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
  Float_t dbposi[3];
  Float_t dbspulsei[2072];

  dbtree->SetBranchAddress("seg",&dbsegi);
  dbtree->SetBranchAddress("pos",dbposi);
  dbtree->SetBranchAddress("spulse",dbspulsei);

  int npoint = dbtree->GetEntriesFast();

  TMatrixD tmppos(3,1);
  int idx[3];
  TMatrixD tmpspulse(nsig*nsegcore,1);
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

    for(int iseg=0; iseg<nsegcore; iseg++){
      int nsigoff = nsig - nsigdb;
      for(int i=0; i<nsigoff; i++)
	tmpspulse(iseg*nsig+i,0) = 0;
      for(int i=0; i<nsigdb; i++)
        tmpspulse(iseg*nsig+nsigoff+i,0)=dbspulsei[iseg*nsigdb+i];
    }

    dbseg.push_back(dbsegi);
    dbpos.push_back(tmppos);
    imap[idx[0]][idx[1]][idx[2]] = ipoint;

    dbspulse.push_back(tmpspulse);
    
    calibpoint.push_back(vector<int>());
  }

  cout<<"load "<<npoint<<" points from "<<DBPosfile<<endl;
  fdb->Close();
  
  return;
}



void GetDBPulse(float* pos, int &seg, int &ngrid, int &extrpl){
  TMatrixD tmppos(3,1);
  int idx[3];
  for(int iaxis=0; iaxis<3; iaxis++){
    tmppos(iaxis,0)=pos[iaxis];
    idx[iaxis] = (int)((pos[iaxis]-range[iaxis][0]) / griddist + 0.5);
  }

  // find grid within 2mm cube
  vector <int> ip; // id list of close grid point for interpolation
  vector <int> ip2; // id list of close grid point 4mm for extrapolation
  vector <double> ipdist;
  vector <double> ip2dist;
  vector <int> exflag;
  double mindist = 5;
  int minip = -1;
  int tmpseg = -1;
  int npoint = dbpos.size();

  int idxrange = 2;
  for(int ix=idx[0]-idxrange; ix<=idx[0]+idxrange; ix++)
    for(int iy=idx[1]-idxrange; iy<=idx[1]+idxrange; iy++)
      for(int iz=idx[2]-idxrange; iz<=idx[2]+idxrange; iz++){
	if(ix<0 || ix>=MaxSteps) continue;
	if(iy<0 || iy>=MaxSteps) continue;
	if(iz<0 || iz>=MaxSteps) continue;

	int ipoint = imap[ix][iy][iz];
	if(ipoint<0) continue;

	if(kextrapol){
	  if(fabs(tmppos(0,0)-dbpos[ipoint](0,0))>4) continue;
	  if(fabs(tmppos(1,0)-dbpos[ipoint](1,0))>4) continue;
	  if(fabs(tmppos(2,0)-dbpos[ipoint](2,0))>4) continue;
	  double disttmp = sqrt((tmppos-dbpos[ipoint]).Sqr().Sum());
	  ip2.push_back(ipoint);
	  ip2dist.push_back(disttmp);
	}

	if(fabs(tmppos(0,0)-dbpos[ipoint](0,0))>2) continue;
	if(fabs(tmppos(1,0)-dbpos[ipoint](1,0))>2) continue;
	if(fabs(tmppos(2,0)-dbpos[ipoint](2,0))>2) continue;
	double disttmp = sqrt((tmppos-dbpos[ipoint]).Sqr().Sum());
	ip.push_back(ipoint);
	ipdist.push_back(disttmp);
	if(disttmp<mindist){
	  mindist = disttmp;
	  minip = ipoint;
	  tmpseg = dbseg[ipoint]; // seg of closest grid point
	}

	if(kextrapol){ ip2.pop_back(); ip2dist.pop_back();}
      }

  if(tmpseg<0){ // remove evt if interaction point outside grid map
    seg = -1;
    return;
  }
  seg = tmpseg; // segment

  // remove grid from different segment
  for(int i=0; i<ip.size();){
    if(dbseg[ip[i]]!=tmpseg){
      ip.erase(ip.begin()+i);
    }else{
      exflag.push_back(0);
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
    extrpl = 0;
    if(ip.size()<8){
      // remove grid from different segment
      for(int i=0; i<ip2.size();){
	if(dbseg[ip2[i]]!=tmpseg){
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
	exflag.push_back(1);
      }
    }
  }
  // -----------------------------------------------------------------

  if(mindist==0){
    ngrid = 1;
    realspulse = dbspulse[minip];

  }else{// more than 1 grid

    // prepare spulse list
    vector <TMatrixD> spulselist;
    vector <TMatrixD> poslist;
    vector <vector<int>> igridlist;
    vector <vector<float>> gridwgtlist;
    for(int i=0; i<ip.size(); i++){
      spulselist.push_back(dbspulse[ip[i]]);
      poslist.push_back(dbpos[ip[i]]);
      igridlist.push_back(vector<int>(1));
      igridlist[i][0] = i;
      gridwgtlist.push_back(vector<float>(1));
      gridwgtlist[i][0] = 1;
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
	    double factori = (tmppos(iaxis,0)-poslist[j](iaxis,0)) / (poslist[i](iaxis,0)-poslist[j](iaxis,0));
	    double factorj = (tmppos(iaxis,0)-poslist[i](iaxis,0)) / (poslist[j](iaxis,0)-poslist[i](iaxis,0));
	    spulselist[i] = factori*spulselist[i] + factorj*spulselist[j];
	    poslist[i](iaxis,0) = tmppos(iaxis,0);

	    if(factori==0){
	      igridlist[i].clear();
	      gridwgtlist[i].clear();
	    }else{
	      for(int ig=0; ig<gridwgtlist[i].size(); ig++)
		gridwgtlist[i][ig] = gridwgtlist[i][ig]*factori;
	    }

	    if(factorj!=0)
	      for(int ig=0; ig<igridlist[j].size(); ig++){
		igridlist[i].push_back(igridlist[j][ig]);
		gridwgtlist[i].push_back(gridwgtlist[j][ig]*factorj);
	      }

	    spulselist.erase(spulselist.begin()+j);
	    poslist.erase(poslist.begin()+j);
	    igridlist.erase(igridlist.begin()+j);
	    gridwgtlist.erase(gridwgtlist.begin()+j);
	    
	  }else{
	    j++;
	  }
	}// end of loop j

	if(!kextrapol) poslist[i](iaxis,0) = tmppos(iaxis,0); // reduced interpolation
	i++;
      }// end of loop i
    }// end of loop iaxis

    if(poslist.size()>1){

      if(!kextrapol){
	cerr<<"something wrong..."<<endl;

      }else{ // if use extrapolation

	while(poslist.size()>0){
	  if(fabs(tmppos(0,0)-poslist[0](0,0))<0.01 &&
	     fabs(tmppos(1,0)-poslist[0](1,0))<0.01 &&
	     fabs(tmppos(2,0)-poslist[0](2,0))<0.01)
	    break;

	  spulselist.erase(spulselist.begin());
	  poslist.erase(poslist.begin());
	  igridlist.erase(igridlist.begin());
	  gridwgtlist.erase(gridwgtlist.begin());
	}

	if(poslist.size()==0){ // cannot get simulated position in XYZ combine
	  ngrid = 0;
	  return;
	}// end of cannot get simulated position in XYZ combine

      }// end of kextrapol

    }// end of poslist>1

    ngrid = igridlist[0].size();
    for(int ig=0; ig<igridlist[0].size(); ig++){
      if(kextrapol){
	if(exflag[igridlist[0][ig]]==1) extrpl++;
      }
    }
    realspulse = spulselist[0];
  }// end of more than 1 grid

  return;
}

