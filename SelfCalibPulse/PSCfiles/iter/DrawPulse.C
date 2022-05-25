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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

bool kextrapol = true;
double chi2l = 9, chi2h = 11;
int skip = 0;
int knext = 0;

const int nsig = 121;
const int nseg = 36;
const int nsegcore = 37;

double griddist  = 2; // 2mm grid
double range[3][2];
const int MaxSteps = 50;
int imap[MaxSteps][MaxSteps][MaxSteps];

// pulsedb
string dbfile = "../../G4Sim/pulsedb/pulseA.root";
vector<int>      dbseg;
vector<TMatrixD> dbpos;
vector<TMatrixD> dbspulse;

// calib pscfile
string pscfile = "it1/Det0000_fit4.root";

// 
TMatrixD realspulse(nsig*nsegcore,1);
TMatrixD calibspulse(nsig*nsegcore,1);


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

void LoadDB();
void GetPulse(float* pos, int &seg, int &ngrid, int &extrpl);

// main
void DrawPulse(){
  float xmin = -0.5, xmax = 36.5;

  // load db
  LoadDB();

  // load calib pscfile
  TChain *tree = new TChain();
  for(int iseg=0; iseg<nseg; iseg++)
    tree->AddFile(pscfile.c_str(),0,Form("tree%d",iseg));

  int det, seg;
  float detpos[3];
  float cadpos2[3];
  float spulse[nsegcore][nsig];
  int npaths;
  tree->SetBranchAddress("det",&det);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("detpos",detpos);
  tree->SetBranchAddress("cadpos2",cadpos2);
  tree->SetBranchAddress("spulse",spulse);
  tree->SetBranchAddress("npaths",&npaths);

  int nentries = tree->GetEntriesFast();

  float phi0, r0, z0; // real
  float phi1, r1, z1; // calib
  float dist;

  int dbsegi;
  int ngrid, extrpl;
  
  Float_t x[4477];
  Float_t y1[4477], y2[4477], y3[4477], y4[4477];

  float ymin3 = 0, ymax3 = 0, ymin4 = -0.01, ymax4 = 0;
  float chi2s[nsegcore];
  float chi2;

  TCanvas *c = new TCanvas("c","c",1000,600);
  TPad *p1 = new TPad("p1","pulse",0,0,1.,0.7);
  p1->SetMargin(0.05,0.05,0.05,0.01);
  p1->Draw();
  c->cd();
  TPad *p2 = new TPad("p2","diff",0,0.7,1.,1.);
  p2->SetMargin(0.05,0.05,0.15,0.05);
  p2->Draw();

  int ientry = 0;
 nextpulse:
  for(; ientry<nentries; ientry++){
    tree->GetEntry(ientry);
    //if(npaths<50) continue;

    TVector3 rpos(detpos[0],detpos[1],0);
    phi0 = rpos.Phi()/TMath::Pi()*180;
    r0 = rpos.Mag();
    z0 = detpos[2];

    TVector3 fpos(cadpos2[0],cadpos2[1],0);
    phi1 = fpos.Phi()/TMath::Pi()*180;
    r1 = fpos.Mag();
    z1 = cadpos2[2];

    dist = (fpos-rpos).Mag();
    
    for(int iseg=0; iseg<nsegcore; iseg++){
      for(int isig=0; isig<nsig; isig++){
	calibspulse(iseg*nsig+isig,0) = spulse[iseg][isig];
      }
    }

  
    //GetPulse(detpos, dbsegi, ngrid, extrpl);
    GetPulse(cadpos2, dbsegi, ngrid, extrpl);

    chi2 = 0;
    ymin3 = 0; ymax3 = 0; ymin4 = -0.01; ymax4 = 0;
    for(int iseg=0; iseg<nsegcore; iseg++){
      chi2s[iseg] = 0;
      for(int isig=0; isig<nsig; isig++){
	int i = iseg*nsig+isig;

	x[i] = i*1./nsig-0.5;
	y1[i] = realspulse(i,0);
	y2[i] = calibspulse(i,0);

	y3[i] = y2[i] - y1[i];
	float sigma = y1[i]>0.01? y1[i] : 0.01;
	//y3[i] = y3[i]*y3[i]/sigma;
	chi2 += y3[i]*y3[i]/sigma;

	if(fabs(y3[i])>ymax3) ymax3 = fabs(y3[i]);
	chi2s[iseg] += y3[i];
      }
      if(chi2s[iseg]==0) chi2s[iseg] = -1;
    }

    for(int iseg=0; iseg<nsegcore; iseg++){
      for(int isig=0; isig<nsig; isig++){
	y4[iseg*nsig+isig] = chi2s[iseg];

	if(chi2s[iseg]>ymax4) ymax4 = chi2s[iseg];
      }
    }

    if(chi2>chi2l&&chi2<chi2h) skip--;
    if(skip>-1) continue;
    else break;

  }

  //ymin3 = ymin3 - 0.001;
  //ymax3 = 1.1*ymax3 + 0.001;
  //ymax4 = 1.1*ymax4 + 0.001;
  ymax3 = 1.1*ymax3;
  ymin3 = -ymax3;

  p1->cd();
  TGraph *gr1 = new TGraph(4477,x,y1);
  gr1->GetXaxis()->SetRangeUser(xmin,xmax);
  gr1->GetYaxis()->SetRangeUser(-0.5,1.1);
  gr1->SetLineWidth(2);
  gr1->SetTitle("");

  TGraph *gr2 = new TGraph(4477,x,y2);
  gr2->SetLineColor(2);
  gr2->SetLineWidth(2);
  
  gr1->Draw("APL");
  gr2->Draw("Lsame");

  TLegend *label = new TLegend(0.7,0.85,0.9,0.99);
  label->AddEntry(gr1,"real pulse shape","l");
  label->AddEntry(gr2,"calib pulse shape","l");
  label->Draw("same");

  p1->Modified();
  p1->Update();

  p2->cd();
  TGraph *gr3 = new TGraph(4477,x,y3);
  gr3->GetYaxis()->SetTitle("diff");
  gr3->GetYaxis()->SetTitleSize(0.1);
  gr3->GetYaxis()->SetTitleOffset(0.2);
  gr3->GetYaxis()->CenterTitle();
  gr3->GetXaxis()->SetRangeUser(xmin,xmax);
  gr3->GetYaxis()->SetRangeUser(ymin3,ymax3);
  gr3->GetXaxis()->SetLabelSize(0.08);
  gr3->GetYaxis()->SetLabelSize(0.08);
  gr3->GetYaxis()->SetNdivisions(505);
  gr3->SetLineWidth(2);
  gr3->SetTitle("");
  gr3->Draw("APL");

  p2->Modified();
  p2->Update();

  c->Modified();
  c->Update();

  cout<<endl<<"ientry "<<ientry<<": det "<<det<<" seg "<<seg<<" npaths "<<npaths<<endl
      <<Form(" realpos: %.1f %.1f %.1f ( PhiRZ: %.1f %.1f %.1f )",detpos[0],detpos[1],detpos[2],phi0,r0,z0)<<endl
      <<Form("calibpos: %.1f %.1f %.1f ( PhiRZ: %.1f %.1f %.1f )",cadpos2[0],cadpos2[1],cadpos2[2],phi1,r1,z1)<<endl
      <<Form("    dist: %.1f",dist)<<endl
      <<"    chi2: "<<chi2<<endl;

  
  cout<<endl<<"dbseg "<<dbsegi<<" ngrid "<<ngrid<<" extrpl "<<extrpl<<endl<<endl;

  cout<<"next pulse?"<<endl;

  cin>>knext;
  if(knext>0){
    ientry++;
    skip = knext-1;
    goto nextpulse;
  }

  return;
}


void LoadDB(){
  range[0][0] = -40.250;    range[0][1] = 39.750;
  range[1][0] = -40.250;    range[1][1] = 39.750;
  range[2][0] = 2.250;      range[2][1] = 90.250;
  for(int ix=0; ix<MaxSteps; ix++)
    for(int iy=0; iy<MaxSteps; iy++)
      for(int iz=0; iz<MaxSteps; iz++)
	imap[ix][iy][iz] = -1;

  TFile *fdb = new TFile(dbfile.c_str());
  if(!fdb->IsOpen()){
    cerr<<"cannot find dbfile "<<dbfile<<endl;
    return;
  }

  TTree *dbtree = (TTree *)fdb->Get("tree");

  Int_t dbsegi;
  Double_t dbposi[3];
  Double_t dbcorei[121];
  Double_t dbspulsei[4356];

  dbtree->SetBranchAddress("seg",&dbsegi);
  dbtree->SetBranchAddress("pos",dbposi);
  dbtree->SetBranchAddress("core",dbcorei);
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

    for(int i=0; i<nsig; i++) tmpspulse(nseg*nsig+i,0)=dbcorei[i];
    for(int iseg=0; iseg<nseg; iseg++){
      for(int i=0; i<nsig; i++)
	tmpspulse(iseg*nsig+i,0)=dbspulsei[iseg*nsig+i];
    }

    dbseg.push_back(dbsegi-1);
    dbpos.push_back(tmppos);
    imap[idx[0]][idx[1]][idx[2]] = ipoint;

    dbspulse.push_back(tmpspulse);
  }

  cout<<"load "<<npoint<<" points from "<<dbfile<<endl;
  fdb->Close();
  
  return;
}


void GetPulse(float* pos, int &seg, int &ngrid, int &extrpl){
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
