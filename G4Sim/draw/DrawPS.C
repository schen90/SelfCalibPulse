#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include <TROOT.h>
#include "TInterpreter.h"

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#define NOISE 1000000
#define NSig 56
#define NSeg 36
#define NSegCore 37

#define NType 3
#define GridDist 2. // 2mm grid
#define GridMaxSteps 50

using namespace std;

#ifdef NOISE
float noise[NOISE]; // base for noise

void MakeNoise(){
  cout<<"AddNoise: create noise base size "<<NOISE<<endl;
  // get random noise
  float noise0[NOISE];
  float tmpnoise;
  for(int i=0; i<NOISE; ){
    tmpnoise = gRandom->Uniform(-1,1);
    for(int isig=0; isig<NSig; isig++){
      if(i>=NOISE) break;
      tmpnoise += -0.2*(tmpnoise+gRandom->Uniform(-10,10));
      noise[i] = noise0[i] = tmpnoise;
      i++;
    }
  }

  // smooth noise
  int nsmooth = 3;
  for(int i=nsmooth; i<NOISE-nsmooth; i++){
    for(int j=1; j<=nsmooth; j++)
      noise[i] += noise0[i-j] + noise0[i+j];

    noise[i] = noise[i] / (2*nsmooth+1) * 2; // scale noise
  }

  return;
}
#endif


// structure for Pules Shape
struct PS{
  int   det;
  int   seg;
  int   nhits;
  int   interid;  // interaction id in a event
  vector<float> hiteng; // hit energy
  float energy; // core energy

  float labpos[3]; // lab position
  float detpos[3]; // det position

  float opulse[NSegCore][NSig]; // original pulse shape
};


struct OBJ{
  // Declaration of leaf types
  int                      ievent;     // event id
  vector<int>             *ndet = 0;   // detector id
  vector<int>             *g4seg = 0;  // segment id
  vector<float>           *energy = 0; // deposit energy

  vector<vector<float>>   *posa = 0;   // absolute/global interaction position vector<float(3)>
  vector<vector<float>>   *posr = 0;   // relative/local interaction position vector<float(3)>

  // pulse shape vector<>
  vector<int>             *pdet = 0;   // detector id for pulse shape
  vector<float>           *ecore = 0;  // core energy
  vector<vector<int>>     *inter = 0;  // interaction id in G4 vector<>

  vector<vector<int>>     *pseg = 0;   // segment id in pulsedb
  vector<vector<int>>     *ngrid = 0;  // number of grid found around ppos
  vector<vector<int>>     *extrpl = 0;  // number of grid for extrapolation
  vector<vector<float>>   *core = 0;   // core pulse shape vector<float(56)>
  vector<vector<float>>   *spulse = 0; // segment pulse shape vector<float(2016)>

  int                      category; // 1: max 1 seg fired in a det, >1 det fired; 2: max >1 seg fired in a det
};

OBJ obj;


// PSbasis
bool kextrapol = true;

double range[NType][3][2];
int imap[NType][GridMaxSteps][GridMaxSteps][GridMaxSteps];

vector<Int_t>    dbseg[NType];
vector<TMatrixD> dbpos[NType];
vector<TMatrixD> dbspulse[NType]; // core at the end

void swap(vector<double> &v,int m, int l){
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

void ReadPSbasis(){
  cout<<"\e[1;31m Read basis for linear interpolation PS: \e[0m"<<endl;
  string dbfile[3] = {"../pulsedb/LibTrap_A001.root",
                      "../pulsedb/LibTrap_B001.root",
                      "../pulsedb/LibTrap_C001.root"};

  for(int itype=0; itype<NType; itype++){

    if(itype!=0) continue;

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


Int_t GetPS(int itype, TMatrixD pos, double energy, int &seg, TMatrixD &spulse){
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


PS GetAPS(int idet, int nidx, int nidxshift){
  PS aps;
  aps.det = -1;
  // calc average position
  vector<int>           simseg;
  vector<int>           siminterid;
  vector<int>           simnhits;
  vector<vector<float>> simhiteng;
  vector<float>         simeng;

  vector<vector<float>> simlabpos;
  vector<vector<float>> simdetpos;

  TMatrixD              simspulse(NSig*NSegCore,1);

  for(int i=0; i<obj.inter->at(idet).size(); i++){
    int interid = obj.inter->at(idet)[i];

    simseg.push_back(obj.pseg->at(idet)[i]); //pseg start from 0...
    siminterid.push_back(interid);
    simnhits.push_back(1);
    vector<float> tmphiteng;
    tmphiteng.push_back(obj.energy->at(interid));
    simhiteng.push_back(tmphiteng);
    simeng.push_back(obj.energy->at(interid));

    vector<float> tmplabpos;
    vector<float> tmpdetpos;
    for(int ii=0; ii<3; ii++){
      tmplabpos.push_back(obj.posa->at(interid)[ii]); // lab position
      tmpdetpos.push_back(obj.posr->at(interid)[ii]); // det position
    }
    simlabpos.push_back(tmplabpos);
    simdetpos.push_back(tmpdetpos);

    int itype = obj.pdet->at(idet)%3;
    TMatrixD tmppos(3,1);
    for(int ix=0; ix<3; ix++) tmppos(ix,0) = tmpdetpos[ix];
    double tmpenergy = tmphiteng[0];
    int tmpseg;
    TMatrixD tmpspulse(NSig*NSegCore,1);
    int tmpngrid = GetPS(itype, tmppos, tmpenergy, tmpseg, tmpspulse);
    if(tmpseg!=obj.pseg->at(idet)[i]){
      cerr<<"segment not match!!! tmpseg = "<<tmpseg
	  <<", simseg = "<<obj.pseg->at(idet)[i]<<endl;
      return aps;
    }
    if(tmpngrid<1) return aps; // cannot get PS
    simspulse = simspulse + tmpspulse;
  }

  // sum up in a segment

  for(int i=0; i<simeng.size(); i++){
    for(int j=i+1; j<simeng.size(); j++){
      if(simseg[i]==simseg[j]){
        siminterid[i] = simeng[i] > simeng[j] ? siminterid[i] : siminterid[j];
        simnhits[i] = simnhits[i] + simnhits[j];
        for(int jj=0; jj<simhiteng[j].size(); jj++) simhiteng[i].push_back(simhiteng[j][jj]);
        double tmpe = simeng[i]+simeng[j];

        for(int ix=0; ix<3; ix++){
          simlabpos[i][ix] = simeng[i]/tmpe*simlabpos[i][ix] + simeng[j]/tmpe*simlabpos[j][ix];
          simdetpos[i][ix] = simeng[i]/tmpe*simdetpos[i][ix] + simeng[j]/tmpe*simdetpos[j][ix];
        }

        simeng[i] = tmpe;

        simseg.erase(simseg.begin()+j);
        siminterid.erase(siminterid.begin()+j);
        simnhits.erase(simnhits.begin()+j);
        simhiteng[j].clear();
        simhiteng.erase(simhiteng.begin()+j);
        simeng.erase(simeng.begin()+j);

        simlabpos.erase(simlabpos.begin()+j);
        simdetpos.erase(simdetpos.begin()+j);

        j--;
      }
    }
  }

  if(simeng.size()!=1){
    return aps;// require only 1 seg fired in a det
  }

  int idx = 0;

  aps.det = obj.pdet->at(idet);
  aps.seg = simseg[idx];
  aps.interid = siminterid[idx];
  aps.nhits = simnhits[idx];
  for(int jj=0; jj<simhiteng[idx].size(); jj++) aps.hiteng.push_back(simhiteng[idx][jj]);
  aps.energy = simeng[idx];

  for(int ix=0; ix<3; ix++){
    aps.labpos[ix] = simlabpos[idx][ix];
    aps.detpos[ix] = simdetpos[idx][ix];
  }

  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++){
      double tmpamp = simspulse(iseg*NSig+isig,0);
      aps.opulse[iseg][isig] = tmpamp;
    }
  }


#ifdef NOISE
  if(nidx>=0){
    // add noise
    for(int iseg=0; iseg<NSegCore; iseg++){
      for(int isig=0; isig<NSig; isig++){
        nidx = nidx%NOISE;
        aps.opulse[iseg][isig] += noise[nidx];
        nidx++;
      }
      nidx+=nidxshift;
    }
  }
#endif

  return aps;
}


void DrawPS(){
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");

  ReadPSbasis();

#ifdef NOISE
  MakeNoise();
#endif

  TChain *tree = new TChain();
  for(int i=0; i<60; i++)
    tree->AddFile(Form("noPS_22Na/G4SimData%04d.root",i),0,"tree");

  // set branch addresses and branch pointers
  tree->SetBranchAddress("ievent",&obj.ievent);
  tree->SetBranchAddress("ndet",&obj.ndet);
  tree->SetBranchAddress("g4seg",&obj.g4seg);
  tree->SetBranchAddress("energy",&obj.energy);

  tree->SetBranchAddress("posa",&obj.posa);
  tree->SetBranchAddress("posr",&obj.posr);

  tree->SetBranchAddress("pdet",&obj.pdet);
  tree->SetBranchAddress("ecore",&obj.ecore);
  tree->SetBranchAddress("inter",&obj.inter);

  tree->SetBranchAddress("pseg",&obj.pseg);
  tree->SetBranchAddress("ngrid",&obj.ngrid);

  tree->SetBranchAddress("category",&obj.category);

  TH2D *hPS0 = new TH2D("hPS0","Core PS",60,-0.5,59.5,1500,-100,1400);
  TH2D *hPS0c = new TH2D("hPS0c","Core PS normal",60,-0.5,59.5,500,-0.25,1.25);

  TH2D *hPS1 = new TH2D("hPS1","Seg0 PS",60,-0.5,59.5,1500,-100,1400);
  TH2D *hPS1c = new TH2D("hPS1c","Seg0 PS normal",60,-0.5,59.5,500,-0.25,1.25);

  Long64_t nentries = tree->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;

  for(int i=0; i<nentries; i++){
    if(i%1000==0) cout<<"\r finish "<<i<<" / "<<nentries<<flush;

    tree->GetEntry(i);

    for(int idet=0; idet<obj.pdet->size(); idet++){ // loop dets
      int detid = obj.pdet->at(idet);

      if(detid!=0) continue;

      int tmpnidx = -1, tmpnidxshift = 0;
#ifdef NOISE
      tmpnidx = (int)gRandom->Uniform(0,NOISE);
      tmpnidxshift = (int)gRandom->Uniform(0,NOISE/(NSig*NSegCore));
#endif
      PS aps = GetAPS(idet,tmpnidx,tmpnidxshift); // aps w/ PS

      if(aps.det<0) continue;
      if(aps.seg!=0) continue;

      for(int isig=0; isig<NSig; isig++){
	hPS0->Fill(isig+60-NSig, aps.opulse[NSegCore-1][isig]);
        hPS0c->Fill(isig+60-NSig, aps.opulse[NSegCore-1][isig]/aps.energy);

        hPS1->Fill(isig+60-NSig, aps.opulse[0][isig]);
        hPS1c->Fill(isig+60-NSig, aps.opulse[0][isig]/aps.energy);
      }
    }
    
  }

  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c","c",600,600);
  c->SetMargin(0.,0.,0.,0.);
  c->Divide(2,2);

  c->cd(1)->SetMargin(0.18,0.01,0.15,0.01);
  c->cd(1)->SetLogz(1);

  hPS0->SetTitle("");
  hPS0->GetXaxis()->CenterTitle();
  hPS0->GetXaxis()->SetTitleSize(0.07);
  hPS0->GetXaxis()->SetTitleOffset(0);
  hPS0->GetXaxis()->SetLabelSize(0.06);
  hPS0->GetXaxis()->SetTitle("Core0");
  hPS0->GetYaxis()->CenterTitle();
  hPS0->GetYaxis()->SetTitleSize(0.07);
  hPS0->GetYaxis()->SetTitleOffset(1.2);
  hPS0->GetYaxis()->SetLabelSize(0.06);
  hPS0->GetYaxis()->SetTitle("Amp");
  hPS0->Draw("col");


  c->cd(2)->SetMargin(0.18,0.01,0.15,0.01);
  c->cd(2)->SetLogz(1);

  hPS0c->SetTitle("");
  hPS0c->GetXaxis()->CenterTitle();
  hPS0c->GetXaxis()->SetTitleSize(0.07);
  hPS0c->GetXaxis()->SetTitleOffset(0);
  hPS0c->GetXaxis()->SetLabelSize(0.06);
  hPS0c->GetXaxis()->SetTitle("Core0");
  hPS0c->GetYaxis()->CenterTitle();
  hPS0c->GetYaxis()->SetTitleSize(0.07);
  hPS0c->GetYaxis()->SetTitleOffset(1.2);
  hPS0c->GetYaxis()->SetLabelSize(0.06);
  hPS0c->GetYaxis()->SetTitle("Amp");
  hPS0c->Draw("col");


  c->cd(3)->SetMargin(0.18,0.01,0.15,0.01);
  c->cd(3)->SetLogz(1);

  hPS1->SetTitle("");
  hPS1->GetXaxis()->CenterTitle();
  hPS1->GetXaxis()->SetTitleSize(0.07);
  hPS1->GetXaxis()->SetTitleOffset(0);
  hPS1->GetXaxis()->SetLabelSize(0.06);
  hPS1->GetXaxis()->SetTitle("Seg0");
  hPS1->GetYaxis()->CenterTitle();
  hPS1->GetYaxis()->SetTitleSize(0.07);
  hPS1->GetYaxis()->SetTitleOffset(1.2);
  hPS1->GetYaxis()->SetLabelSize(0.06);
  hPS1->GetYaxis()->SetTitle("Amp");
  hPS1->Draw("col");


  c->cd(4)->SetMargin(0.18,0.01,0.15,0.01);
  c->cd(4)->SetLogz(1);

  hPS1c->SetTitle("");
  hPS1c->GetXaxis()->CenterTitle();
  hPS1c->GetXaxis()->SetTitleSize(0.07);
  hPS1c->GetXaxis()->SetTitleOffset(0);
  hPS1c->GetXaxis()->SetLabelSize(0.06);
  hPS1c->GetXaxis()->SetTitle("Seg0");
  hPS1c->GetYaxis()->CenterTitle();
  hPS1c->GetYaxis()->SetTitleSize(0.07);
  hPS1c->GetYaxis()->SetTitleOffset(1.2);
  hPS1c->GetYaxis()->SetLabelSize(0.06);
  hPS1c->GetYaxis()->SetTitle("Amp");
  hPS1c->Draw("col");


  return;
}
