#ifndef TREEREADERPULSE_CC
#define TREEREADERPULSE_CC

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <time.h>
#include <thread>

#include <TMath.h>
#include <TVector3.h>

#include "Global.hh"
#include "AGATA.hh"
#include "AGATAgeo.hh"
#include "TreeReaderPulse.hh"
#include "Tracker.hh"

using namespace std;

TreeReaderPulse::TreeReaderPulse(int detid){

  Detid = detid;
  
#ifdef ADDPS
  apsb = new PSbasis(detid);
#endif

#ifdef NOISE
  LoadNoise();
#endif
  
  for(int i=0; i<NChain; i++){
    fChain[i] = new TChain();
  }

  ievt = 0;
  kcout = true;
  minchi2 = 1e9;
  kInterrupt = 0;
}


TreeReaderPulse::~TreeReaderPulse(){
}


void TreeReaderPulse::Load(string configfile){
  nConfig = 0;
  cout<<"\e[1;31m Configure file: "<<configfile<<"\e[0m"<<endl;

  // read configure file
  char *buffer = new char[500];
  ifstream fin(configfile.c_str());

  float sE; // source energy keV
  float spos[3]; // source position mm
  string pathtmp;
  int run[2];
  int nevt;
  while(!fin.eof()){
    fin.getline(buffer,500);
    if(strncmp(buffer,"#input",6)==0){
      fin >> buffer >> sE >> spos[0] >> spos[1] >> spos[2];
      fin >> buffer >> pathtmp;
      fin >> buffer >> run[0] >> run[1];
      fin >> buffer >> nevt;

      fSourceE.push_back(sE);
      fSourcePos.push_back(TVector3(spos[0],spos[1],spos[2]));
      path.push_back(pathtmp);
      MinRun.push_back(run[0]);
      MaxRun.push_back(run[1]);
      Nevts.push_back(nevt);
      nConfig++;
    }
  }
  fin.close();

  // read input configures
  if(nConfig<1){
    cerr<<"cannot find input from "<<configfile<<endl;
    return;
  }
  
  cout<<"\e[1;33m find "<<nConfig<<" inputs:\e[0m"<<endl;
  for(int i=0; i<nConfig; i++){
    cout<<"#input "<<i<<": "
	<<Form("source %.1fkeV at %.2f %.2f %.2f",
	       fSourceE[i],
	       fSourcePos[i].X(),fSourcePos[i].Y(),fSourcePos[i].Z())
	<<endl;

    for(int run=MinRun[i]; run<=MaxRun[i]; run++){
      fChain[0]->AddFile(Form("%s/G4SimData%04d.root",path[i].c_str(),run),0,"tree");
    }

    int nentries = fChain[0]->GetEntriesFast();
    cout<<" find \e[1m"<<nentries<<"\e[0m events from rootfiles "
	<<Form("%s/G4SimData%04d ~ %04d",path[i].c_str(),MinRun[i],MaxRun[i])<<endl;

    fChain[0]->Reset();
  }

  // init fChain[0] for  #input 0
  for(int run=MinRun[0]; run<=MaxRun[0]; run++){
    fChain[0]->AddFile(Form("%s/G4SimData%04d.root",path[0].c_str(),run),0,"tree");
  }
  Init(0);

  SourceE = fSourceE[0];
  SourcePos= fSourcePos[0];
}

void TreeReaderPulse::ScanPS(AGATA *agata, int nevts){
  ScanPS(agata, nevts, -1);
}

void TreeReaderPulse::ScanPS(AGATA *agata, int nevts, double Diff){
  //fChain[0]->GetEntry(0);
  //int tdet = obj[0].pdet->at(0);
  //cout<<"compare pulse shape for det "<<tdet<<endl;

  //int itype = tdet%3;

  int nentries = fChain[0]->GetEntriesFast();
  cout<<"find "<<nentries<<" events from tree"<<endl;
  irun = MinRun[0];
  ievt = 0;
  kInterrupt = 0;

#ifndef NTHREADS
  ScanPSLoop1(0, agata, nevts);

#else
  // loop trees with multi threads
  thread th[NTHREADS];
  cout<<"using "<<NTHREADS<<" threads:"<<endl;

  for(int i=0; i<NTHREADS; i++){
    th[i] = thread(&TreeReaderPulse::ScanPSLoop1, this, i, ref(agata), nevts);
  }

  for(int i=0; i<NTHREADS; i++){
    if(th[i].joinable())
      th[i].join();
  }
#endif
  
  cout<<"\r load "<<fPSs[0].size()<<"__"<<fPSs[1].size()<<"__"<<fPSs[2].size()<<" / "<<nevts<<" pulse shape"<<endl;
  
  sort(fPSs[0].begin(),fPSs[0].end(),[](const PS& lhs, const PS& rhs){return lhs.seg<rhs.seg;});
  sort(fPSs[1].begin(),fPSs[1].end(),[](const PS& lhs, const PS& rhs){return lhs.seg<rhs.seg;});
  sort(fPSs[2].begin(),fPSs[2].end(),[](const PS& lhs, const PS& rhs){return lhs.seg<rhs.seg;});

  // output
  TFile *fout = new TFile("ComparePS.root","RECREATE");

#ifdef REALPOS
  TTree *postree[3];
  for(int itype=0; itype<3; itype++)
    postree[itype] = new TTree(Form("postree%d",itype),Form("PS position tree of type%d",itype));
#endif
  TTree *anatree[3];
  for(int itype=0; itype<3; itype++)
    anatree[itype] = new TTree(Form("tree%d",itype),Form("analyzed tree of type%d",itype));


  
#ifndef NTHREADS
  for(int itype=0; itype<3; itype++)
    ScanPSLoop2(0, postree[itype], anatree[itype], agata, nevts, Diff);

#else
  thread th2[NTHREADS];
  cout<<"using "<<NTHREADS<<" threads:"<<endl;

  for(int itype=0; itype<3; ){
    for(int i=0; i<NTHREADS; i++){
      if(itype<3)
	th2[i] = thread(&TreeReaderPulse::ScanPSLoop2, this, itype, ref(postree[itype]), ref(anatree[itype]), ref(agata), nevts, Diff);
      itype++;
    }

    for(int i=0; i<NTHREADS; i++){
      if(th2[i].joinable())
	th2[i].join();
    }
  }
#endif
  

  fout->cd();
#ifdef REALPOS
  for(int itype=0; itype<3; itype++) postree[itype]->Write();
#endif
  for(int itype=0; itype<3; itype++) anatree[itype]->Write();
  fout->Close();

  return;
}

void TreeReaderPulse::ScanPSLoop1(int iChain, AGATA *agata, int nevts){
  int run;
  for(; irun<=MaxRun[0]; ){ // loop runs
    if(kInterrupt) break;

    {
#ifdef NTHREADS
      lock_guard<mutex> lock(treemtx); // lock tree read
#endif
      if(irun>MaxRun[0]) return;
      run = irun;
      irun++;
    }

    // initial tree
    fChain[iChain]->Reset();
    fChain[iChain]->AddFile(Form("%s/G4SimData%04d.root",path[0].c_str(),run),0,"tree");
    int nentrytmp = fChain[iChain]->GetEntriesFast();
    Init(iChain);

    for(int ientry=0; ientry<nentrytmp; ientry++){ // loop evts
      if(kInterrupt) break;

      {
#ifdef NTHREADS
	lock_guard<mutex> lock(scanmtx); // lock scan
#endif
	if(ievt%10000==0) cout<<"\r load "<<fPSs[0].size()<<"__"<<fPSs[1].size()<<"__"<<fPSs[2].size()<<" / "<<nevts<<" pulse shape... ievt = "<<ievt<<flush;
      }
      
      fChain[iChain]->GetEntry(ientry);
      for(int idet=0; idet<obj[iChain].pdet->size(); idet++){
	//if(obj[iChain].pdet->at(idet)!=tdet) continue;
	int itype = obj[iChain].pdet->at(idet)%3;

	int tmpnidx = -1, tmpnidxshift = 0;
#ifdef NOISE
	tmpnidx = (int)gRandom->Uniform(0,NOISE);
	tmpnidxshift = (int)gRandom->Uniform(0,NOISE/(NSig*NSegCore));
#endif
	int segidx = -1;
	PS aps = GetAPS(iChain,agata,idet,tmpnidx,tmpnidxshift,false,segidx); // aps w/ PS
	if(aps.det<0) continue;

#ifdef SINGLEHIT
	if(aps.nhits>1) continue;
#endif
	{
#ifdef NTHREADS
	  lock_guard<mutex> lock(scanmtx); // lock scan
#endif
	  if(fPSs[itype].size()<nevts) fPSs[itype].push_back(aps);
	}
      }

      ievt++;

      {
#ifdef NTHREADS
	lock_guard<mutex> lock(scanmtx); // lock scan
#endif
	if(fPSs[0].size()>=nevts && fPSs[1].size()>=nevts && fPSs[2].size()>=nevts){
	  kInterrupt = 1;
	  break;
	}
      }

      
    } // end of loop evts
    
  } // end of loop runs

  return;  
}


void TreeReaderPulse::ScanPSLoop2(int itype, TTree *postree, TTree *anatree, AGATA *agata, int nevts, double Diff){
  AGATAgeo* agatageo = agata->GetGeo();

  Int_t simseg;
  Int_t nhits1, nhits2;
  vector<float> hiteng1;
  vector<float> hiteng2;
  Float_t SimPos[3];
  Float_t Phi, Radius, Z, PhiC, ZC;
  Float_t PhiRZ1[3], PhiRZ2[3];
  Float_t dist;
  Float_t chi2;
  Float_t chi2s[3];
  Int_t nfired;
  Float_t diffphi, diffr, diffz, rdiffphi;
  Float_t Energy1, Energy2;
  Float_t pulse1[NSegCore][NSig], pulse2[NSegCore][NSig], chis[NSegCore][NSig], zero[NSegCore][NSig];
  for(int iseg=0; iseg<NSegCore; iseg++)
    for(int isig=0; isig<NSig; isig++)
      zero[iseg][isig] = 0;

#ifdef REALPOS
  postree->Branch("simseg", &simseg, "simseg/I");
  postree->Branch("nhits", &nhits1, "nhits/I");
  postree->Branch("SimPos", SimPos, "SimPos[3]/F");
#endif

  anatree->Branch("type", &itype, "type/I");
  anatree->Branch("simseg", &simseg, "simseg/I");
  anatree->Branch("SimPos", SimPos, "SimPos[3]/F");
  anatree->Branch("Phi", &Phi, "Phi/F");
  anatree->Branch("Radius", &Radius, "Radius/F");
  anatree->Branch("Z", &Z, "Z/F");
  anatree->Branch("PhiC", &PhiC, "PhiC/F"); // phi relative to seg center
  anatree->Branch("ZC", &ZC, "ZC/F"); // z relative to seg center
  anatree->Branch("chi2", &chi2, "chi2/F");
  anatree->Branch("chi2s", chi2s, "chi2s[3]/F");
  anatree->Branch("nfired", &nfired, "nfired/I");
#ifdef REALPOS
  anatree->Branch("dist", &dist, "dist/F");
  anatree->Branch("diffphi", &diffphi, "diffphi/F");
  anatree->Branch("diffr", &diffr, "diffr/F");
  anatree->Branch("diffz", &diffz, "diffz/F");
  anatree->Branch("rdiffphi", &rdiffphi, "rdiffphi/F");
#endif
  anatree->Branch("nhits1", &nhits1, "nhits1/I");
  anatree->Branch("nhits2", &nhits2, "nhits2/I");
  anatree->Branch("Energy1", &Energy1, "Energy1/F");
  anatree->Branch("Energy2", &Energy2, "Energy2/F");
  if(nevts<=1000){
    anatree->Branch("hiteng1",&hiteng1);
    anatree->Branch("hiteng2",&hiteng2);
    anatree->Branch("PhiRZ1", PhiRZ1, "PhiRZ1[3]/F");
    anatree->Branch("PhiRZ2", PhiRZ2, "PhiRZ2[3]/F");
    anatree->Branch("pulse1", pulse1, Form("pulse1[%d][%d]/F",NSegCore,NSig));
    anatree->Branch("pulse2", pulse2, Form("pulse2[%d][%d]/F",NSegCore,NSig));
    anatree->Branch("chis", chis, Form("chis[%d][%d]/F",NSegCore,NSig));
  }

  cout<<"start compare..."<<endl;
  float aspulse[NSig_comp], bspulse[NSig_comp];

  int iievt;
  for(iievt=0; iievt<fPSs[itype].size(); iievt++){
    if(iievt%1000==0) cout<<"\r type "<<itype<<": finish "<<iievt<<" / "<<fPSs[itype].size()<<" pulse..."<<flush;

    simseg = fPSs[itype][iievt].seg;
    Energy1 = fPSs[itype][iievt].energy;

    if(Energy1<PSCEMIN) continue;

    nhits1 = fPSs[itype][iievt].nhits;
    if(nevts<=1000){
      hiteng1.clear();
      for(int ii=0; ii<fPSs[itype][iievt].hiteng.size(); ii++) hiteng1.push_back(fPSs[itype][iievt].hiteng[ii]);
    }

    TMatrixD SegPos(3,1);
    SegPos = agatageo->GetLocalSegPos(itype,simseg);
    TVector3 segvec(SegPos(0,0),SegPos(1,0),0);
    float SegPhi = segvec.Phi()/TMath::Pi()*180;
    float SegR   = segvec.Mag();
    float SegZ   = SegPos(2,0);
#ifdef REALPOS
    for(int ix=0; ix<3; ix++) SimPos[ix] = fPSs[itype][iievt].detpos[ix];
    TVector3 ivec(fPSs[itype][iievt].detpos[0],fPSs[itype][iievt].detpos[1],0);
    PhiRZ1[0] = ivec.Phi()/TMath::Pi()*180;
    PhiRZ1[1] = ivec.Mag();
    PhiRZ1[2] = fPSs[itype][iievt].detpos[2];

    Phi    = PhiRZ1[0];    
    Radius = PhiRZ1[1];    
    Z      = PhiRZ1[2];
    PhiC   = Phi - SegPhi; if(PhiC>180) PhiC-=360; if(PhiC<-180) PhiC+=360;
    ZC     = Z - SegZ;
    postree->Fill();
#endif

    int fseg[NSeg_comp]; //0,1:fired seg, core; 2,3:next sectors; 4,5:next slice
    agatageo->GetNextSegs(simseg, fseg);

    for(int jevt = iievt+1; jevt<fPSs[itype].size(); jevt++){
      if(simseg!=fPSs[itype][jevt].seg) break;

      Energy2 = fPSs[itype][jevt].energy;
      nhits2 = fPSs[itype][jevt].nhits;
      if(nevts<=1000){
	hiteng2.clear();
	for(int jj=0; jj<fPSs[itype][jevt].hiteng.size(); jj++) hiteng2.push_back(fPSs[itype][jevt].hiteng[jj]);
      }

#ifdef REALPOS
      dist = 0;
      for(int ix=0; ix<3; ix++) dist += pow(fPSs[itype][iievt].detpos[ix]-fPSs[itype][jevt].detpos[ix],2);
      dist = sqrt(dist);

      TVector3 jvec(fPSs[itype][jevt].detpos[0],fPSs[itype][jevt].detpos[1],0);
      PhiRZ2[0] = jvec.Phi()/TMath::Pi()*180;
      PhiRZ2[1] = jvec.Mag();
      PhiRZ2[2] = fPSs[itype][jevt].detpos[2];

      diffphi = PhiRZ1[0] - PhiRZ2[0];
      if(diffphi>180)  diffphi-=360;
      if(diffphi<-180) diffphi+=360;
      diffphi = fabs(diffphi);
      rdiffphi = (PhiRZ1[1]+PhiRZ2[1])/2*diffphi/180*TMath::Pi();
      diffr = fabs(PhiRZ1[1] - PhiRZ2[1]);
      diffz = fabs(PhiRZ1[2] - PhiRZ2[2]);

      if(Diff>0){
	double difflimit = 0.3;
	int ndiff = 0;
	if(diffr>difflimit) ndiff++;
	if(diffz>difflimit) ndiff++;
	if(rdiffphi>difflimit) ndiff++;
	if(ndiff!=1) continue; // select data diff only in one dimension

	ndiff = 0;
	if(fabs(diffr-Diff)<0.2) ndiff++;
	if(fabs(diffz-Diff)<0.2) ndiff++;
	if(fabs(rdiffphi-Diff)<0.2) ndiff++;
	if(ndiff!=1) continue; // select data only match Diff
      }
#endif

      nfired=0;
      chi2 = 0;
      for(int ix=0; ix<3; ix++) chi2s[ix]=0;
      int uflg[NSegCore]; for(int iseg=0; iseg<NSegCore; iseg++) uflg[iseg]=0;

      for(int ix=0; ix<3; ix++){
	for(int ii=0; ii<2; ii++){
	  int iseg = fseg[2*ix+ii];
	  if(uflg[iseg]!=0) continue;
	  if(nevts<=1000){
	    copy_n(fPSs[itype][iievt].opulse[iseg], NSig, pulse1[iseg]);
	    copy_n(fPSs[itype][jevt].opulse[iseg], NSig, pulse2[iseg]);
	    copy_n(zero[iseg], NSig, chis[iseg]);
	  }

	  copy_n(fPSs[itype][iievt].apulse[2*ix+ii], NSig_comp, aspulse);
	  copy_n(fPSs[itype][jevt].apulse[2*ix+ii], NSig_comp, bspulse);

	  if(nevts<=1000){
	    for(int isig=0; isig<NSig_comp; isig++){
	      chis[iseg][isig] = pow(aspulse[isig]-bspulse[isig],2); //SQ
	      float sigm = fabs(bspulse[isig]); if(sigm<0.01) sigm=0.01;
	      chis[iseg][isig] = chis[iseg][isig]/sigm; //chi2
	    }
	  }
	  
	  float tmpchi2 = agata->Chi2seg(aspulse, bspulse);
	  chi2 += tmpchi2;
	  chi2s[ix] += tmpchi2; // sum
	  nfired++;

	  uflg[iseg]=1;
	}
      }

      // compare other segments
      for(int iseg=0; iseg<NSegCore; iseg++){
	if(uflg[iseg]!=0) continue;
	if(nevts<=1000){
	  copy_n(fPSs[itype][iievt].opulse[iseg], NSig, pulse1[iseg]);
	  copy_n(fPSs[itype][jevt].opulse[iseg], NSig, pulse2[iseg]);
	  copy_n(zero[iseg], NSig, chis[iseg]);
	}

	copy_n(fPSs[itype][iievt].opulse[iseg], NSig_comp, aspulse);
	copy_n(fPSs[itype][jevt].opulse[iseg], NSig_comp, bspulse);

	if(nevts<=1000){
	  for(int isig=0; isig<NSig_comp; isig++){
	    chis[iseg][isig] = pow(aspulse[isig]-bspulse[isig],2); //SQ
	    float sigm = fabs(bspulse[isig]); if(sigm<0.01) sigm=0.01;
	    chis[iseg][isig] = chis[iseg][isig]/sigm; //chi2
	  }
	}

	float tmpchi2 = agata->Chi2seg(aspulse, bspulse);
	//tmpchi2 = fPSs[itype][iievt].segwgt[iseg]>0? tmpchi2*fPSs[itype][iievt].segwgt[iseg] : tmpchi2*fPSs[itype][jevt].segwgt[iseg];
	//if(tmpchi2>chi2) chi2=tmpchi2; // maximum
	chi2 += tmpchi2; // sum
	nfired++;

	uflg[iseg]=1;
      }
      //if(nfired>0) chi2 = chi2/nfired;

      anatree->Fill();
    }
    
  }
  cout<<"\r type "<<itype<<": finish "<<iievt<<" / "<<fPSs[itype].size()<<" pulse..."<<endl;
  
  return;
}


void TreeReaderPulse::MakeInit(){
  cout<<"Make folders for SelfCalib..."<<endl;
  gROOT->ProcessLine(".!rm -rf ./share/*");
  gROOT->ProcessLine(".!mkdir ./share/Hits");
  gROOT->ProcessLine(".!mkdir ./share/HCs");
  //gROOT->ProcessLine(".!mkdir ./share/PSCs");

  for(int i=0; i<nConfig; i++){
    gROOT->ProcessLine(Form(".!mkdir ./share/Hits/input%d",i));
  }

  MakeNoise();
  return;
}


void TreeReaderPulse::MakeNoise(){

#ifdef NOISE
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

  ofstream fnoise("./share/NoiseBase.txt");
  for(int i=0; i<NOISE; i++) fnoise<<noise[i]<<endl;
  fnoise.close();
#endif

  return;
}


void TreeReaderPulse::LoadNoise(){

#ifdef NOISE
  if(gSystem->AccessPathName("./share/NoiseBase.txt")) MakeNoise();
  
  ifstream fnoise("./share/NoiseBase.txt");
  for(int i=0; i<NOISE; i++) fnoise>>noise[i];
  fnoise.close();
#endif
  
  return;
}


void TreeReaderPulse::Init(int i){
  if (i<0 || i>=NChain) return;

  // set branch addresses and branch pointers
  fChain[i]->SetBranchAddress("ievent",&obj[i].ievent);
  fChain[i]->SetBranchAddress("ndet",&obj[i].ndet);
  fChain[i]->SetBranchAddress("g4seg",&obj[i].g4seg);
  fChain[i]->SetBranchAddress("energy",&obj[i].energy);
#ifdef REALPOS
  fChain[i]->SetBranchAddress("posa",&obj[i].posa);
  fChain[i]->SetBranchAddress("posr",&obj[i].posr);
#endif
  fChain[i]->SetBranchAddress("pdet",&obj[i].pdet);
  fChain[i]->SetBranchAddress("ecore",&obj[i].ecore);
  fChain[i]->SetBranchAddress("inter",&obj[i].inter);

  fChain[i]->SetBranchAddress("pseg",&obj[i].pseg);
  fChain[i]->SetBranchAddress("ngrid",&obj[i].ngrid);

#ifndef ADDPS
  if(kWithPS){
    //fChain[i]->SetBranchAddress("extrpl",&obj[i].extrpl);
    fChain[i]->SetBranchAddress("core",&obj[i].core);
    fChain[i]->SetBranchAddress("spulse",&obj[i].spulse);
  }
#endif
  
  fChain[i]->SetBranchAddress("category",&obj[i].category);

  return;
}


void TreeReaderPulse::GenerateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
					int ientry, int nentries){

  // get entry-------------------------------------------------------------
  {
#ifdef NTHREADS
    lock_guard<mutex> lock(treemtx); // lock tree read
#endif
    if(ievt>=nentries) return;
    if(ientry>=nentries) return;  

    ievt++;
  }
  fChain[iChain]->GetEntry(ientry);

  // get PS in one event-----------------------------------
  vector<PS> fPS;
  vector<int> fSegIdx;
  vector<int> Nidx;
  vector<int> Nidxshift;
  for(int idet=0; idet<obj[iChain].pdet->size(); idet++){ //loop dets
    int detid = obj[iChain].pdet->at(idet);

    if(Detid>-1 && detid!=Detid) continue; // get PS for selected Detid

    int tmpnidx = -1, tmpnidxshift = 0;
#ifdef NOISE
    tmpnidx = (int)gRandom->Uniform(0,NOISE);
    tmpnidxshift = (int)gRandom->Uniform(0,NOISE/(NSig*NSegCore));
#endif
    int segidx = -1;
    PS aps = GetAPS(iChain,agata,idet,tmpnidx,tmpnidxshift,false,segidx); // aps w/ PS

    if(aps.det<0 && segidx>-1){ // multi-segment fired
#ifdef MULTISEG
      for(int idx=0; idx<segidx; idx++){
	PS aps = GetAPS(iChain,agata,idet,-1,0,true,idx); // aps w/o PS
        if(aps.det<0) continue;
#ifdef SINGLEHIT
	if(aps.nhits>1) continue;
#endif
	fPS.push_back(aps);
	fSegIdx.push_back(idx);
	Nidx.push_back(-1);
	Nidxshift.push_back(0);
      }
#endif
      
    }else{ // one segment fired
      if(aps.det<0) continue;
#ifdef SINGLEHIT
      if(aps.nhits>1) continue;
#endif
      fPS.push_back(aps);
      fSegIdx.push_back(segidx);
      Nidx.push_back(tmpnidx);
      Nidxshift.push_back(tmpnidxshift);
    }

  }//end of loop dets

  
  if(fPS.size()==0) return; // skip if selected Det not fired

  if(Detid>-1){ // one det mode
    // get other PSs in one event-----------------------------------
    for(int idet=0; idet<obj[iChain].pdet->size(); idet++){ //loop dets
      int detid = obj[iChain].pdet->at(idet);

      if(Detid>-1 && detid==Detid) continue; // get PS for other Dets

      int segidx = -1;
      PS aps = GetAPS(iChain,agata,idet,-1,0,true,segidx); // aps w/o PS

      if(aps.det<0 && segidx>-1){ // multi-segment fired
#ifdef MULTISEG
	for(int idx=0; idx<segidx; idx++){
	  PS aps = GetAPS(iChain,agata,idet,-1,0,true,idx); // aps w/o PS
	  if(aps.det<0) continue;
#ifdef SINGLEHIT
	  if(aps.nhits>1) continue;
#endif
	  fPS.push_back(aps);
	  fSegIdx.push_back(idx);
	}
#endif
	
      }else{ // one segment fired
        if(aps.det<0) continue;
#ifdef SINGLEHIT
	if(aps.nhits>1) continue;
#endif
	fPS.push_back(aps);
	fSegIdx.push_back(segidx);

      }

    }//end of loop dets
  }

  if(fPS.size()<2) return; // at least one Compton scattering

  // create Hit-----------------------------------------------
  vector<int> uflag;
  EventHits* fEvent = new EventHits(SourceE, SourcePos);
  fEvent->SetIdx(iconfig,run,ientry);

  for(int i=0; i<fPS.size(); i++){ //loop fPS
#ifdef REALPOS
    TVector3 hitpos(fPS[i].labpos[0],fPS[i].labpos[1],fPS[i].labpos[2]);
#else
    TVector3 hitpos(0,0,0);
#endif
    TVector3 initpos;
    if(Detid<0 || fPS[i].det==Detid){
      initpos = agata->GetPSpos(fPS[i].det, fPS[i].seg, &fPS[i]);
    }else{
      TMatrixD SegPos = agata->GetGeo()->GetSegPos(fPS[i].det,fPS[i].seg);
      initpos.SetXYZ(SegPos(0,0), SegPos(1,0), SegPos(2,0));
    }
    
    Hit *ahit = new Hit(fPS[i].det, fPS[i].seg, fPS[i].energy, hitpos, initpos); //keV and mm
    ahit->SetInterid(fPS[i].interid);
    if(Detid<0 || fPS[i].det==Detid){
#ifdef NOISE
      ahit->SetNoiseIdx(Nidx[i]);
      ahit->SetNoiseIdxShift(Nidxshift[i]);
#endif
    }
    
    fEvent->Add(ahit);
    uflag.push_back(1);
  }

  int iEvtHit = agata->AddEventHits(fEvent);
  vector<Hit*>* fHits = fEvent->GetfHits();
  
#ifdef CHECKTRACK
  // check track----------------------------------------------
  Tracker tracker(fHits, (Double_t)SourceE, SourcePos);
  tracker.OFTtracking();
  vector<int> atrack = tracker.GetTrack();
  for(int i=0; i<fPS.size(); i++) uflag[i] = 0;
  if(atrack.size()>1) for(int i=0; i<atrack.size(); i++) uflag[atrack[i]] = 1;
#endif
  
  // group PS-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    if(Detid>-1 && fPS[i].det!=Detid) continue; // one det mode
    if(fSegIdx[i]>-1) continue; // multi segment fired
    if(uflag[i]!=1) continue;
    vector<int> entrylist;
    double tmpminchi2 = agata->AddPStoPSC(&fPS[i], fHits->at(i), entrylist); // add to pulse shape collection
    if(tmpminchi2<minchi2) minchi2 = tmpminchi2;
  }//end of loop fPS
  
  return;

}

void TreeReaderPulse::GenerateHCsLoop(int iconfig, int iChain, AGATA *agata, int nentries){
  
  int PSCstat[10];
  time(&start);

  int istart = 0;
  int run;
  for(; irun<=MaxRun[iconfig]; ){ //loop runs
    if(kInterrupt) break;

    {
#ifdef NTHREADS
      lock_guard<mutex> lock(treemtx); // lock tree read
#endif
      if(irun>MaxRun[iconfig]) return;
      run = irun;
      irun++;
    }

    // initial tree
    fChain[iChain]->Reset();
    fChain[iChain]->AddFile(Form("%s/G4SimData%04d.root",path[iconfig].c_str(),run),0,"tree");
    int nentrytmp = fChain[iChain]->GetEntriesFast();
    Init(iChain);

    
    for(int ientry=0; ientry<nentrytmp; ientry++){ //loop evts
      if(kInterrupt) break;
    
      // output state
      if(ievt%10000==0 && kcout){
	kcout = false;      
	time(&stop);

#ifdef NTHREADS
	lock_guard<mutex> lock(treemtx); // lock tree read
#endif
	agata->GetPSCstat(PSCstat);
	double MemUsageGB = GetCurrentMemoryUsage()/GB;
	double MemTotalGB = GetTotalSystemMemory()/GB;
	double MemUsage = MemUsageGB / MemTotalGB * 100;

	if(kUpdateHCs==0){
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start))
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]<<" Mem-"<<PSCstat[2]<<".." //" File-"<<PSCstat[3]<<".."
	      <<"maxnhits-"<<PSCstat[4]
	      <<Form(" minchi2-%.3f",minchi2)<<flush;

	}else{ // update HCs
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start))
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]<<" Mem-"<<PSCstat[2]<<".." //" File-"<<PSCstat[3]<<".."
	      <<"RemoveHit-"<<cRemoveHit<<" AddHit-"<<cAddHit<<" NotMatch-"<<cNotMatch<<".."
	      <<"maxnhits-"<<PSCstat[4]<<flush;
	}
      
	if(MemUsage>MaxMemUsage){
	  cout<<endl<<"exceed memory limit. Write to PSCfiles..."<<endl;
	  kInterrupt = 1;
	  break;
	}

	time(&start);
	kcout = true;
      }

      // work on entry
      if(kUpdateHCs==0) GenerateHCsworker(iconfig, run, iChain, agata, ientry, nentries);
      else                UpdateHCsworker(iconfig, run, iChain, agata, ientry, nentries, istart);

    }//end of loop evts

  }//end of loop runs
  
  return;
}


void TreeReaderPulse::GenerateHCs(AGATA *agata){
  if(kUpdateHCs>0){
    cRemoveHit = 0;
    cAddHit    = 0;
    cNotMatch  = 0;
  }

  if(nConfig<1){
    cerr<<"cannot find input, check configure file..."<<endl;
    return;
  }

  for(int i=0; i<nConfig; i++){
    cout<<"\e[1;31m"<<"#input "<<i<<"\e[0m"<<endl;
    SourceE = fSourceE[i];
    SourcePos = fSourcePos[i];
    GenerateHCs(agata, Nevts[i], i);
  }

  
#ifdef NTHREADS
  if(kUpdateHCs==0){
    //sort EventHits
    cout<<endl<<"\r sort fEventHits..."<<flush;
    time(&start);
    agata->SortEventHits();
    time(&stop);
    cout<<Form("\r sort fEventHits..%.0fs",difftime(stop,start))<<endl;
  }
#endif
  
}


void TreeReaderPulse::GenerateHCs(AGATA *agata, int nevts){
  if(kUpdateHCs>0){
    cRemoveHit = 0;
    cAddHit    = 0;
    cNotMatch  = 0;
  }

  GenerateHCs(agata, nevts, 0);

  
#ifdef NTHREADS
  if(kUpdateHCs==0){
    //sort EventHits
    cout<<endl<<"\r sort fEventHits..."<<flush;
    time(&start);
    agata->SortEventHits();
    time(&stop);
    cout<<Form("\r sort fEventHits..%.0fs",difftime(stop,start))<<endl;
  }
#endif
  
}


void TreeReaderPulse::GenerateHCs(AGATA *agata, int nevts, int iconfig){

  if(kUpdateHCs>0){
    cout<<"\e[1;31m UpdateHCs Level "<<kUpdateHCs<<" ... \e[0m"<<endl;
  }
  
  if(!kGroupPos && !kWithPS){
    cerr<<"No Pulse Shape read from tree!!! kGroupPos="<<kGroupPos<<"; kWithPS="<<kWithPS<<endl;
    return;
  }

  // statistics in total
  fChain[0]->Reset();
  for(int run=MinRun[iconfig]; run<=MaxRun[iconfig]; run++){
    fChain[0]->AddFile(Form("%s/G4SimData%04d.root",path[iconfig].c_str(),run),0,"tree");
  }
  Init(0);

  int nentries = fChain[0]->GetEntriesFast();
  if(nevts>0) nentries = TMath::Min(nentries,nevts);
  irun = MinRun[iconfig];
  ievt = 0;

  cout<<"\e[1;33m Read "<<nentries<<" events from rootfiles "<<Form("%s/G4SimData%04d ~ %04d",path[iconfig].c_str(),MinRun[iconfig],MaxRun[iconfig])<<" ... \e[0m"<<endl;

  // Loop entries-------------------------------------------
#ifndef NTHREADS
  GenerateHCsLoop(iconfig, 0, agata, nentries);

#else
  // loop trees with multi threads
  thread th[NTHREADS];
  cout<<"using "<<NTHREADS<<" threads:"<<endl;
  
  for(int i=0; i<NTHREADS; i++){
    th[i] = thread(&TreeReaderPulse::GenerateHCsLoop, this, iconfig, i, ref(agata),nentries);
  }

  for(int i=0; i<NTHREADS; i++){
    if(th[i].joinable())
      th[i].join();
  }

#endif

  // output final statistics
  NEventHits = agata->GetEventHitsSize();
  
  time(&stop);
  int PSCstat[10];
  agata->GetPSCstat(PSCstat);
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;

  if(kUpdateHCs==0){
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start))
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]<<" Mem-"<<PSCstat[2]<<".." //" File-"<<PSCstat[3]<<".."
	<<"maxnhits-"<<PSCstat[4]
	<<Form(" minchi2-%.3f..",minchi2)
	<<"fEventHits-"<<NEventHits<<endl;

  }else{ // update HCs
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start))
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]<<" Mem-"<<PSCstat[2]<<".." //" File-"<<PSCstat[3]<<".."
	<<"RemoveHit-"<<cRemoveHit<<" AddHit-"<<cAddHit<<" NotMatch-"<<cNotMatch<<".."
	<<"maxnhits-"<<PSCstat[4]<<endl;
  }

  return;
}



void TreeReaderPulse::UpdateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
				      int ientry, int nentries, int &istart){

  if(kUpdateHCs<1) return;

  if(ievt>=nentries) return;
  
  // find Hit-----------------------------------------------
  int iEvtHit = agata->FindiEvtHit(iconfig, run, ientry, istart);
  if(iEvtHit<0 || iEvtHit>NEventHits-1){ ievt++; return;}

  vector<Hit*>* fHits = agata->FindEventHits(iEvtHit)->GetfHits();
  vector<int> uflag;
  for(int i=0; i<fHits->size(); i++){ //loop fPS
    uflag.push_back(1);
  }
  istart = iEvtHit;

#ifdef CHECKTRACK
  // check track----------------------------------------------
  Tracker tracker(fHits, (Double_t)SourceE, SourcePos);
  tracker.OFTtracking();
  vector<int> atrack = tracker.GetTrack();
  for(int i=0; i<fHits->size(); i++) uflag[i] = 0;
  if(atrack.size()>1) for(int i=0; i<atrack.size(); i++) uflag[atrack[i]] = 1;
#endif
  
  bool kFindPS = false;
  // check if need to compare hits with HCs
  for(int i=0; i<fHits->size(); i++){
    if(fHits->at(i)->GetLevel()!=kUpdateHCs) continue;

    if(Detid>-1 && fHits->at(i)->GetDet()!=Detid) continue; // selected Detid
      
    if(kNewPSC){ // if New PSC added from last step, compare again all hits with HCs
      kFindPS = true;

    }else if(uflag[i]==0 && fHits->at(i)->hasHitCollection()>0){ // remove from PSCs
      kFindPS = true;

    }else if(uflag[i]==1 && fHits->at(i)->hasHitCollection()==0){ // add to PSCs
      kFindPS = true;
    }
  }
  
  if(!kFindPS){ ievt++; return;}

  // find corresponding fPS
  vector<PS> fPS;
  {
#ifdef NTHREADS
    lock_guard<mutex> lock(treemtx); // lock tree read
#endif
    if(ievt>=nentries) return;
    if(ientry>=nentries) return;  

    ievt++;
  }
  fChain[iChain]->GetEntry(ientry);

  int ihit = 0;
  for(int idet=0; idet<obj[iChain].pdet->size(); idet++){ //loop dets
    int detid = obj[iChain].pdet->at(idet);

    if(Detid>-1 && detid!=Detid) continue; // get PS for selected Detid
    
    int tmpnidx = -1, tmpnidxshift=0;
#ifdef NOISE
    if(ihit >= fHits->size()) ihit=fHits->size()-1;
    tmpnidx = fHits->at(ihit)->GetNoiseIdx();
    tmpnidxshift = fHits->at(ihit)->GetNoiseIdxShift();
#endif
    PS aps = GetAPS(iChain,agata,idet,tmpnidx,tmpnidxshift); // aps w/ PS
    if(aps.det<0) return;
#ifdef SINGLEHIT
    if(aps.nhits>1) continue;
#endif
    fPS.push_back(aps);
    ihit++;
  }//end of loop dets
  
  // check fPS match with fHits
  cNotMatch++;
  if(Detid<0 && fHits->size()!=fPS.size()) return;
  for(int i=0; i<fPS.size(); i++){
    if(fHits->at(i)->GetDet()!=fPS[i].det || fHits->at(i)->GetSeg()!=fPS[i].seg) return;
    if(fHits->at(i)->GetE()!=fPS[i].energy) return;
#ifdef REALPOS
    TVector3 hitpos(fPS[i].labpos[0],fPS[i].labpos[1],fPS[i].labpos[2]);
    if(fHits->at(i)->GetRealPosition()!=hitpos) return;
#endif
  }
  cNotMatch--;

  // Update PSCs-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS
    if(fHits->at(i)->GetLevel()!=kUpdateHCs) continue;
    
    if(kNewPSC){ // if New PSC added from last step, compare again all hits with HCs
      if(uflag[i]!=1) continue;
      vector<int> entrylist;
      double tmpminchi2 = agata->AddPStoPSC(&fPS[i], fHits->at(i), entrylist); // add to pulse shape collection
      if(entrylist.size()>0)
	cAddHit++;

    }else if(uflag[i]==0 && fHits->at(i)->hasHitCollection()>0){ // remove from PSCs
      agata->RemovePSfromPSC(&fPS[i], fHits->at(i)); // remove a PS from pulse shape collection
      cRemoveHit++;

    }else if(uflag[i]==1 && fHits->at(i)->hasHitCollection()==0){ // add to PSCs
      vector<int> entrylist;
      double tmpminchi2 = agata->AddPStoPSC(&fPS[i], fHits->at(i), entrylist); // add to pulse shape collection
      if(entrylist.size()>0)
	cAddHit++;
    }
    
  }//end of loop fPS

  return;

}


PS TreeReaderPulse::GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift){
  return GetAPS(iChain, agata, idet, nidx, nidxshift, false);
}

PS TreeReaderPulse::GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift, bool skipPS){
  int segidx = -1;
  return GetAPS(iChain, agata, idet, nidx, nidxshift, skipPS, segidx);
}

PS TreeReaderPulse::GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift, bool skipPS, int &segidx){
  PS aps;
  aps.det = -1;
  // calc average position
  vector<int>           simseg;
  vector<int>           siminterid;
  vector<int>           simnhits;
  vector<vector<float>> simhiteng;
  vector<float>         simeng;
#ifdef REALPOS
  vector<vector<float>> simlabpos;
  vector<vector<float>> simdetpos;
#endif
#ifdef ADDPS
  TMatrixD              simspulse(NSig*NSegCore,1);
#endif
  for(int i=0; i<obj[iChain].inter->at(idet).size(); i++){
    int interid = obj[iChain].inter->at(idet)[i];

    simseg.push_back(obj[iChain].pseg->at(idet)[i]); //pseg start from 0...
    siminterid.push_back(interid);
    simnhits.push_back(1);
    vector<float> tmphiteng;
    tmphiteng.push_back(obj[iChain].energy->at(interid));
    simhiteng.push_back(tmphiteng);
    simeng.push_back(obj[iChain].energy->at(interid));
#ifdef REALPOS
    vector<float> tmplabpos;
    vector<float> tmpdetpos;
    for(int ii=0; ii<3; ii++){
      tmplabpos.push_back(obj[iChain].posa->at(interid)[ii]); // lab position
      tmpdetpos.push_back(obj[iChain].posr->at(interid)[ii]); // det position
    }
    simlabpos.push_back(tmplabpos);
    simdetpos.push_back(tmpdetpos);
#endif
#ifdef ADDPS
    if(!skipPS){
      int itype = obj[iChain].pdet->at(idet)%3;
      TMatrixD tmppos(3,1);
      for(int ix=0; ix<3; ix++) tmppos(ix,0) = tmpdetpos[ix];
      double tmpenergy = tmphiteng[0];
      int tmpseg;
      TMatrixD tmpspulse(NSig*NSegCore,1);
      int tmpngrid = apsb->GetPS(itype, tmppos, tmpenergy, tmpseg, tmpspulse);
      if(tmpseg!=obj[iChain].pseg->at(idet)[i]){
	cerr<<"segment not match!!! tmpseg = "<<tmpseg
	    <<", simseg = "<<obj[iChain].pseg->at(idet)[i]<<endl;
	return aps;
      }
      if(tmpngrid<1) return aps; // cannot get PS
      simspulse = simspulse + tmpspulse;
    }
#endif
  }

  // sum up in a segment
  
  for(int i=0; i<simeng.size(); i++){
    for(int j=i+1; j<simeng.size(); j++){
      if(simseg[i]==simseg[j]){
	siminterid[i] = simeng[i] > simeng[j] ? siminterid[i] : siminterid[j];
	simnhits[i] = simnhits[i] + simnhits[j];
	for(int jj=0; jj<simhiteng[j].size(); jj++) simhiteng[i].push_back(simhiteng[j][jj]);
	double tmpe = simeng[i]+simeng[j];
#ifdef REALPOS
	for(int ix=0; ix<3; ix++){
	  simlabpos[i][ix] = simeng[i]/tmpe*simlabpos[i][ix] + simeng[j]/tmpe*simlabpos[j][ix];
	  simdetpos[i][ix] = simeng[i]/tmpe*simdetpos[i][ix] + simeng[j]/tmpe*simdetpos[j][ix];
	}
#endif
	simeng[i] = tmpe;

	simseg.erase(simseg.begin()+j);
	siminterid.erase(siminterid.begin()+j);
	simnhits.erase(simnhits.begin()+j);
	simhiteng[j].clear();
	simhiteng.erase(simhiteng.begin()+j);
	simeng.erase(simeng.begin()+j);
#ifdef REALPOS
	simlabpos.erase(simlabpos.begin()+j);
	simdetpos.erase(simdetpos.begin()+j);
#endif
	j--;
      }
    }
  }

  if(segidx<0 && simeng.size()!=1){
    segidx = simeng.size(); // multi-segment fired
    return aps;// require only 1 seg fired in a det
  }

  int idx;
  if(segidx<0) idx = 0;
  else         idx = segidx;
  
  aps.det = obj[iChain].pdet->at(idet);
  aps.seg = simseg[idx];
  aps.interid = siminterid[idx];
  aps.nhits = simnhits[idx];
  for(int jj=0; jj<simhiteng[idx].size(); jj++) aps.hiteng.push_back(simhiteng[idx][jj]);
  aps.energy = simeng[idx];
#ifdef REALPOS
  for(int ix=0; ix<3; ix++){
    aps.labpos[ix] = simlabpos[idx][ix];
    aps.detpos[ix] = simdetpos[idx][ix];
  }
#endif

  if(skipPS || segidx>-1) return aps; // skip PS
  
  if(kWithPS){
#ifdef ADDPS
    for(int iseg=0; iseg<NSegCore; iseg++){
      for(int isig=0; isig<NSig; isig++){
	double tmpamp = simspulse(iseg*NSig+isig,0);
	aps.opulse[iseg][isig] = tmpamp/simeng[0];
      }
    }
  
#else
    // read segments
    for(int iseg=0; iseg<NSeg; iseg++){
      for(int isig=0; isig<NSig; isig++){
	double tmpamp = obj[iChain].spulse->at(idet)[iseg*NSig+isig];
	aps.opulse[iseg][isig] = tmpamp/simeng[0];
      }
    }

    // read core
    for(int isig=0; isig<NSig; isig++){
      double tmpamp = obj[iChain].core->at(idet)[isig];
      aps.opulse[NSegCore-1][isig] = tmpamp/simeng[0];
    }
#endif

#ifdef NOISE
    if(nidx>=0){
      // add noise
      for(int iseg=0; iseg<NSegCore; iseg++){
	for(int isig=0; isig<NSig; isig++){
	  nidx = nidx%NOISE;
	  aps.opulse[iseg][isig] += noise[nidx]/simeng[0];
	  nidx++;
	}
	nidx+=nidxshift;
      }
    }
#endif

    // pulse shape for comparison
    int fseg[NSeg_comp]; //0,1:fired seg, core; 2,3:next sectors; 4,5:next slice
    agata->GetGeo()->GetNextSegs(aps.seg, fseg);

    for(int iseg=0; iseg<NSeg_comp; iseg++){
      copy_n(aps.opulse[fseg[iseg]], NSig, aps.apulse[iseg]);
      aps.segwgt[iseg]=1;
    }
    
  }
  
  return aps;
}


Double_t TreeReaderPulse::GetTotalSystemMemory(){
  Long64_t pages = sysconf(_SC_PHYS_PAGES);
  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) pages)*page_size;
}

Double_t TreeReaderPulse::GetCurrentMemoryUsage(){
  Int_t tSize, resident, share;
  ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) resident) * page_size;
}

#endif // ifndef TREEREADERPULSE_CC
