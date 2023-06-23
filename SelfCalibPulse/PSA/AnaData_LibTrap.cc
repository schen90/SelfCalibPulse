#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <TMath.h>
#include <TVector3.h>
#include "TRandom.h"
#include "TSystem.h"
#include <chrono>
#include <time.h>

#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TMatrixD.h"
#include "TInterpreter.h"

#include "Global.hh"
#include "PSAFilter.cc"

using namespace std;

// read tracefile
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

// structure for Pules Shape
struct PS{
  Int_t CrystalId;
  float SegE[NSEGS]; // segment energy
  float CoreE; // core energy

  float opulse[NCHAN][BSIZE];   // original pulse shape
};

const int fstep = 2; // mm fine grid
static const int GridMaxSteps = 50;
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

void ReadPSbasis();
Int_t GetPSbasis(int itype, TMatrixD pos, double energy, int &seg, TMatrixD &spulse);
PS GetAPS(int idet);

void AnaData(string tracefile, string psafile){

  ReadPSbasis();

  TFile *fin = new TFile(tracefile.c_str());
  if(!fin->IsOpen()){
    cerr<<"cannot find tracefile "<<tracefile<<endl;
    return;
  }
  TTree *intree = (TTree *)fin->Get("tree");
  // set branch addresses and branch pointers
  intree->SetBranchAddress("ievent",&obj.ievent);
  intree->SetBranchAddress("ndet",&obj.ndet);
  intree->SetBranchAddress("g4seg",&obj.g4seg);
  intree->SetBranchAddress("energy",&obj.energy);

  intree->SetBranchAddress("posa",&obj.posa);
  intree->SetBranchAddress("posr",&obj.posr);

  intree->SetBranchAddress("pdet",&obj.pdet);
  intree->SetBranchAddress("ecore",&obj.ecore);
  intree->SetBranchAddress("inter",&obj.inter);

  intree->SetBranchAddress("pseg",&obj.pseg);
  intree->SetBranchAddress("ngrid",&obj.ngrid);

  intree->SetBranchAddress("category",&obj.category);
  int nentries = intree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from "<<tracefile<<endl;


  // output psafile
  TFile *fout = new TFile(psafile.c_str(),"RECREATE");
  TTree *anatree = new TTree("tree","PSA tree");
  Float_t SegE[36], CoreE;
  Int_t   CrystalId;
  Int_t   numNetCharges;
  Int_t   seg;
  Float_t NetCh;
  Float_t pos[3];
  Float_t segpos[3];
  Float_t chi2;
  anatree->Branch("EntryID",   &obj.ievent,   "EntryID/I");
  anatree->Branch("SegE",       SegE,         "SegE[36]/F");
  anatree->Branch("CoreE",     &CoreE,        "CoreE/F");
  anatree->Branch("CrystalId", &CrystalId,    "CrystalId/I");

  anatree->Branch("numNetCharges", &numNetCharges, "numNetCharges/I");
  anatree->Branch("seg",   &seg,   "seg/I");
  anatree->Branch("NetCh", &NetCh, "NetCh/F");
  anatree->Branch("pos",    pos,   "pos[3]/F");
  anatree->Branch("segpos", segpos,"segpos[3]/F");

  anatree->Branch("chi2",  &chi2,  "chi2/F");


  // PSA
  PSAFilter *fpsa[NDets];
  for(int idet=0; idet<NDets; idet++){
    cout<<"PSA for CrystalId "<<DetId[idet]<<" : ";
    fpsa[idet] = new PSAFilter(ADLpath+ADLfile[idet]);
  }

  // clock
  time_t start, stop;
  time(&start);

  // loop all entries
  int ientry;
  for(ientry=0; ientry<nentries; ientry++){
    if(ientry%1000==0)
      cout<<"\r finish "<<ientry<<" / "<<nentries<<" entries..."<<flush;

    intree->GetEntry(ientry);

    for(int iidet=0; iidet<obj.pdet->size(); iidet++){ // loop dets

      CrystalId = obj.pdet->at(iidet);

      // check the basis loaded
      int idet = -1;
      for(int id=0; id<NDets; id++){
	if(DetId[id]==CrystalId){
	  idet = id;
	  break;
	}
      }
      if(idet<0) continue;

      /////////////////////////////////////////
      // read data
      /////////////////////////////////////////

      PS aps = GetAPS(iidet);

      if(aps.CrystalId<0) continue;
      CrystalId = aps.CrystalId;
      for(int iseg=0; iseg<NSEGS; iseg++){
	SegE[iseg] = aps.SegE[iseg];
      }
      CoreE = aps.CoreE;

      pointExp pE;
      // read SegTrace
      for(int iseg=0; iseg<NCHAN; iseg++){
	float last = 0;
	int ann = 0;
	int bnn = BZERO;
	for(; ann<BSIZE; ann++, bnn++){
	  if(bnn < BSIZE){
	    last = aps.opulse[iseg][bnn];
	  }
	  pE.tAmp[iseg][ann] = last;
	}
      }

      // find NetCharge segments
      int numsegs = 0;
      for(int iseg=0; iseg<NSEGS; iseg++){
	if(aps.SegE[iseg]>15){
	  pE.netChargeSegnum[numsegs] = iseg;
	  pE.netChargeEnergy[numsegs] = aps.SegE[iseg];
	  numsegs++;
	}
      }
      pE.numNetCharges = numsegs;
      pE.netChSeg = -1;
      pE.bestPt = -1;
      pE.chi2min = float(1.e30);
      pE.isValid = false;
      pE.isInitialized = false;

      ////////////////////////////////////
      // grid search
      ////////////////////////////////////
      fpsa[idet]->ProcessOneEvent(pE);

      ////////////////////////////////////
      // output results
      ////////////////////////////////////
      numNetCharges = pE.numNetCharges;
      chi2 = fpsa[idet]->GetChi2(pE);
      for(int snum=0; snum<numNetCharges; snum++){
	seg   = pE.netChargeSegnum[snum];
	NetCh = pE.netChargeEnergy[snum];

	int bestPt = pE.resPt[snum];
	fpsa[idet]->GetPtPos(seg, bestPt, pos);
	fpsa[idet]->GetPtPos(seg, -1, segpos); //segment center

	anatree->Fill();
      }
    }
  }// end of loop entries
  cout<<"\r finish "<<ientry<<" / "<<nentries<<" entries.  "<<endl;

  time(&stop);
  printf("============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  cout<<"write to "<<psafile<<" ..."<<endl;
  fout->cd();
  anatree->Write();
  fout->Close();
  fin->Close();

  return;
}


void ReadPSbasis(){

  cout<<"\e[1;31m Read basis for linear interpolation PS: \e[0m"<<endl;
  string dbfile[3] = {"../G4Sim/pulsedb/LibTrap_A001.root",
		      "../G4Sim/pulsedb/LibTrap_B001.root",
		      "../G4Sim/pulsedb/LibTrap_C001.root"};

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
    Float_t dbspulsei[BSIZE*NCHAN];

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
    TMatrixD tmpspulse(BSIZE*NCHAN,1);
    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);

      for(int i=0; i<3; i++){
        tmppos(i,0)=dbposi[i];
        idx[i] = (int)((dbposi[i]-range[itype][i][0]) / fstep + 0.5);
        if(idx[i]<0 || idx[i]>=GridMaxSteps){
          cerr<<"grid point outside Map range!!!"<<endl;
          return;
        }
      }

      for(int iseg=0; iseg<NSEGS; iseg++){
        for(int i=0; i<BSIZE; i++)
          tmpspulse(iseg*BSIZE+i,0)=dbspulsei[iseg*BSIZE+i];
      }
      for(int i=0; i<BSIZE; i++) tmpspulse(NSEGS*BSIZE+i,0)=dbspulsei[NSEGS*BSIZE+i];

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


Int_t GetPSbasis(int itype, TMatrixD pos, double energy, int &seg, TMatrixD &spulse) {

  int ngrid = 0;

  int idx[3];
  for(int ix=0; ix<3; ix++){
    idx[ix] = (int)((pos(ix,0)-range[itype][ix][0]) / fstep + 0.5);
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


PS GetAPS(int idet){
  PS aps;
  for(int iseg=0; iseg<NSEGS; iseg++) aps.SegE[iseg]=0;
  aps.CoreE = 0;
  aps.CrystalId = -1;
  // calc average position
  vector<int>           simseg;
  vector<int>           siminterid;
  vector<int>           simnhits;
  vector<vector<float>> simhiteng;
  vector<float>         simeng;

  vector<vector<float>> simlabpos;
  vector<vector<float>> simdetpos;

  TMatrixD              simspulse(BSIZE*NCHAN,1);

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
    TMatrixD tmpspulse(BSIZE*NCHAN,1);

    int tmpngrid = GetPSbasis(itype, tmppos, tmpenergy, tmpseg, tmpspulse);

    if(tmpseg!=obj.pseg->at(idet)[i]){
      cerr<<"segment not match!!! tmpseg = "<<tmpseg
	  <<", simseg = "<<obj.pseg->at(idet)[i]<<endl;
      return aps;
    }
    if(tmpngrid<1) return aps; // cannot get PS
    simspulse = simspulse + tmpspulse;
  }

  aps.CrystalId = obj.pdet->at(idet);

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


  for(int i=0; i<simeng.size(); i++){
    int seg = simseg[i];
    aps.SegE[seg] = simeng[i];
    aps.CoreE += simeng[i];
  }

  
  for(int iseg=0; iseg<NCHAN; iseg++){
    for(int isig=0; isig<BSIZE; isig++){
      double tmpamp = simspulse(iseg*BSIZE+isig,0);
      aps.opulse[iseg][isig] = tmpamp;
    }
  }
  

  // add noise
  float noise[NCHAN*BSIZE],noise0[NCHAN*BSIZE];
  float tmpnoise;
  for(int i=0; i<NCHAN*BSIZE; ){
    tmpnoise = gRandom->Uniform(-1,1);
    for(int isig=0; isig<BSIZE; isig++){
      if(i>=NCHAN*BSIZE) break;
      tmpnoise += -0.2*(tmpnoise+gRandom->Uniform(-10,10));
      noise[i] = noise0[i] = tmpnoise;
      i++;
    }
  }

  // smooth noise
  int nsmooth = 3;
  for(int i=nsmooth; i<NCHAN*BSIZE-nsmooth; i++){
    for(int j=1; j<=nsmooth; j++)
      noise[i] += noise0[i-j] + noise0[i+j];

    noise[i] = noise[i] / (2*nsmooth+1) * 2; // scale noise
  }

  for(int iseg=0; iseg<NCHAN; iseg++){
    for(int isig=0; isig<BSIZE; isig++){
      aps.opulse[iseg][isig] += noise[iseg*BSIZE+isig];
    }
  }


  return aps;
}



#ifndef __CINT__
int main(int argc, char *argv[]){
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");


  if(argc!=3){
    cout << "Usage: \n"
	 << argv[0]
	 << " trace-file psa-file\n";
    return 0;
  }

  AnaData(string(argv[1]), string(argv[2]));

  return 0;
}
#endif
