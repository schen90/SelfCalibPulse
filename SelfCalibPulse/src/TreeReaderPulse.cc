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
  ClearSkipDetId();
  
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
  long long nevt;
  while(!fin.eof()){
    fin.getline(buffer,500);
    if(strncmp(buffer,"#input",6)==0){
      Config aconfig;
      int nsource = 0;
      vector<float>    tmpSE;
      vector<TVector3> tmpSPos;
      vector<int>      tmpskip;
      int tmpdetid;
      fin >> buffer >> nsource;
      for(int i=0; i<nsource; i++){
	fin >> buffer >> sE >> spos[0] >> spos[1] >> spos[2];
	tmpSE.push_back(sE);
	tmpSPos.push_back(TVector3(spos[0],spos[1],spos[2]));
      }
      fin >> buffer >> pathtmp;
      fin >> buffer >> run[0] >> run[1];
      // skip det
      fin >> buffer >> tmpdetid;
      while(tmpdetid>-1){
        tmpskip.push_back(tmpdetid);
        fin >> tmpdetid;
      }
      fin >> buffer >> nevt;

      aconfig.NSource    = nsource;
      aconfig.fSourceE   = tmpSE;
      aconfig.fSourcePos = tmpSPos;
      aconfig.path       = pathtmp;
      aconfig.MinRun     = run[0];
      aconfig.MaxRun     = run[1];
      aconfig.Nevts      = nevt;
      aconfig.skipdetid  = tmpskip;

      fConfigs.push_back(aconfig);
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
    cout<<"#input "<<i<<":  "<<fConfigs[i].NSource<<" sources"<<endl;
    for(int is=0; is<fConfigs[i].NSource; is++){
      cout<<Form("     source %.1fkeV at %.2f %.2f %.2f",
		 fConfigs[i].fSourceE[is],
		 fConfigs[i].fSourcePos[is].X(),fConfigs[i].fSourcePos[is].Y(),fConfigs[i].fSourcePos[is].Z())<<endl;
    }

    for(int run=fConfigs[i].MinRun; run<=fConfigs[i].MaxRun; run++){
      fChain[0]->AddFile(Form("%s/G4SimData%04d.root",fConfigs[i].path.c_str(),run),0,"tree");
    }

    long long nentries = fChain[0]->GetEntriesFast();
    cout<<" find \e[1m"<<nentries<<"\e[0m events from rootfiles "
	<<Form("%s/G4SimData%04d ~ %04d",fConfigs[i].path.c_str(),fConfigs[i].MinRun,fConfigs[i].MaxRun)<<endl;

    if(fConfigs[i].skipdetid.size()>0){
      cout<<"\e[1;31m Skip Det: ";
      for(int idet : fConfigs[i].skipdetid){
        int cluster = idet/3;
        string detname = Form("%02d",cluster);
        int itype   = idet%3;
        if(itype==0) detname += "A ";
        if(itype==1) detname += "B ";
        if(itype==2) detname += "C ";
        cout<<detname;
      }
      cout<<"\e[0m"<<endl;
    }

    fChain[0]->Reset();
  }

  // init fChain[0] for  #input 0
  for(int run=fConfigs[0].MinRun; run<=fConfigs[0].MaxRun; run++){
    fChain[0]->AddFile(Form("%s/G4SimData%04d.root",fConfigs[0].path.c_str(),run),0,"tree");
  }
  Init(0);

  SourceE   = fConfigs[0].fSourceE;
  SourcePos = fConfigs[0].fSourcePos;

  ClearSkipDetId();
  for(int idet : fConfigs[0].skipdetid) SkipDetId(idet);

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
    for(int isig=0; isig<BSIZE; isig++){
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

  fChain[i]->SetBranchAddress("posa",&obj[i].posa);
  fChain[i]->SetBranchAddress("posr",&obj[i].posr);

  fChain[i]->SetBranchAddress("pdet",&obj[i].pdet);
  fChain[i]->SetBranchAddress("ecore",&obj[i].ecore);
  fChain[i]->SetBranchAddress("inter",&obj[i].inter);

  fChain[i]->SetBranchAddress("pseg",&obj[i].pseg);
  fChain[i]->SetBranchAddress("ngrid",&obj[i].ngrid);

#ifndef ADDPS
  //fChain[i]->SetBranchAddress("extrpl",&obj[i].extrpl);
  fChain[i]->SetBranchAddress("core",&obj[i].core);
  fChain[i]->SetBranchAddress("spulse",&obj[i].spulse);
#endif
  
  fChain[i]->SetBranchAddress("category",&obj[i].category);

  return;
}


// loop opt=0: generate initial PSC; opt=1: findmaxdev; opt=2: dividepsc
void TreeReaderPulse::GenerateHCs(int opt, AGATA *agata){
  cnevents = 0;
  cievt = 0;
  cievthitfind = 0;
  cievthitnotfind = 0;

  if     (opt==0) cout<<"\e[1;31m Generate initial PSC ... \e[0m"<<endl;
  else if(opt==1) cout<<"\e[1;31m Find Max Deviation ... \e[0m"<<endl;
  else if(opt==2) cout<<"\e[1;31m Divide PSC ... \e[0m"<<endl;
  else if(opt==3) cout<<"\e[1;31m Find Dev Sigma ... \e[0m"<<endl;
  else if(opt==4) cout<<"\e[1;31m Remove strange PS from PSCs ... \e[0m"<<endl;

  if(opt>0){
    cNotMatch  = 0;
  }

  if(opt==4){
    cRemovePS = 0;
  }

  if(nConfig<1){
    cerr<<"cannot find input, check configure file..."<<endl;
    return;
  }

  for(int i=0; i<nConfig; i++){
    cout<<"\e[1;32m"<<" #input "<<i<<"\e[0m"<<endl;
    SourceE = fConfigs[i].fSourceE;
    SourcePos = fConfigs[i].fSourcePos;
    ClearSkipDetId();
    for(int idet : fConfigs[i].skipdetid) SkipDetId(idet);

    GenerateHCs(opt, agata, fConfigs[i].Nevts, i);
  }

  cout<<endl<<"\e[1;32m"<<" #input 0 ~ "<<nConfig-1<<"\e[0m"<<endl;
  long long PSCstat[10];
  agata->GetPSCstat(PSCstat);
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;

  if(opt==0){ // initial PSC
    cout<<"\r finish read "<<cievt<<" / "<<cnevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."
	<<"fEventHits-"<<NEventHits<<".."<<endl;

  }else if(opt==1){ // find max dev
    cout<<"\r finish read "<<cievt<<" / "<<cnevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."
	<<"MaxDev-"<<Form("%.2f",(float)MaxDev)<<".."<<endl;
    cout<<" find iEvtHit: "<<cievthitfind
        <<"   Not find iEvtHit: "<<cievthitnotfind<<endl;

  }else if(opt==2){ // divide PSC
    cout<<"\r finish read "<<cievt<<" / "<<cnevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" DivPS-"<<cDivPS<<" NotMatch-"<<cNotMatch<<".."
	<<" maxnhitsdiv-"<<maxnhitsdiv
	<<" maxndiv-"<<PSCstat[8]<<".."<<endl;
    cout<<" find iEvtHit: "<<cievthitfind
        <<"   Not find iEvtHit: "<<cievthitnotfind<<endl;

  }else if(opt==3){ // find dev sigma
    cout<<"\r finish read "<<cievt<<" / "<<cnevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."<<endl;
    cout<<" find iEvtHit: "<<cievthitfind
        <<"   Not find iEvtHit: "<<cievthitnotfind<<endl;

  }else if(opt==4){ // remove PS
    cout<<"\r finish read "<<cievt<<" / "<<cnevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" RemovePS-"<<cRemovePS<<" NotMatch-"<<cNotMatch<<".."<<endl;
    cout<<" find iEvtHit: "<<cievthitfind
        <<"   Not find iEvtHit: "<<cievthitnotfind<<endl;

  }

  
#ifdef NTHREADS
  if(opt==0){
    //sort EventHits
    cout<<endl<<"\r sort fEventHits..."<<flush;
    time(&start);
    agata->SortEventHits();
    time(&stop);
    cout<<Form("\r sort fEventHits..%.0fs",difftime(stop,start))<<endl;
  }
#endif
  
}


void TreeReaderPulse::GenerateHCs(int opt, AGATA *agata, long long nevts){
  if(opt>0){
    cNotMatch  = 0;
  }

  GenerateHCs(opt, agata, nevts, 0);

  
#ifdef NTHREADS
  if(opt==0){
    //sort EventHits
    cout<<endl<<"\r sort fEventHits..."<<flush;
    time(&start);
    agata->SortEventHits();
    time(&stop);
    cout<<Form("\r sort fEventHits..%.0fs",difftime(stop,start))<<endl;
  }
#endif
  
}


void TreeReaderPulse::GenerateHCs(int opt, AGATA *agata, long long nevts, int iconfig){

  if(opt==1){
    MaxDev = 0;
  }
  if(opt==2){
    cDivPS = 0;
    maxnhitsdiv = 0;
    agata->SetMaxNDiv(0);
  }
  
  // statistics in total
  fChain[0]->Reset();
  for(int run=fConfigs[iconfig].MinRun; run<=fConfigs[iconfig].MaxRun; run++){
    fChain[0]->AddFile(Form("%s/G4SimData%04d.root",fConfigs[iconfig].path.c_str(),run),0,"tree");
  }
  Init(0);

  long long nentries = fChain[0]->GetEntriesFast();
  if(nevts>0) nentries = TMath::Min(nentries,nevts);
  cnevents += nentries;
  irun = fConfigs[iconfig].MinRun;
  ievt = 0;

  cout<<"\e[1;33m Read "<<nentries<<" events from rootfiles "<<Form("%s/G4SimData%04d ~ %04d",fConfigs[iconfig].path.c_str(),fConfigs[iconfig].MinRun,fConfigs[iconfig].MaxRun)<<" ... \e[0m"<<endl;

  // Loop entries-------------------------------------------
#ifndef NTHREADS
  GenerateHCsLoop(opt, iconfig, 0, agata, nentries);

#else
  // loop trees with multi threads
  thread th[NTHREADS];
  cout<<"using "<<NTHREADS<<" threads:"<<endl;
  
  for(int i=0; i<NTHREADS; i++){
    th[i] = thread(&TreeReaderPulse::GenerateHCsLoop, this, opt, iconfig, i, ref(agata),nentries);
  }

  for(int i=0; i<NTHREADS; i++){
    if(th[i].joinable())
      th[i].join();
  }

#endif

  // output final statistics
  cievt += ievt;
  NEventHits = agata->GetEventHitsSize();
  
  time(&stop);
  long long PSCstat[10];
  agata->GetPSCstat(PSCstat);
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;


  if(opt==0){ // initial PSC
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."
	<<"fEventHits-"<<NEventHits<<".."<<endl;

  }else if(opt==1){ // find max dev
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."
	<<"MaxDev-"<<Form("%.2f",(float)MaxDev)<<".."<<endl;

  }else if(opt==2){ // divide PSC
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" DivPS-"<<cDivPS<<" NotMatch-"<<cNotMatch<<".."
	<<" maxnhitsdiv-"<<maxnhitsdiv
	<<" maxndiv-"<<PSCstat[8]<<".."<<endl;

  }else if(opt==3){ // find dev sigma
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."<<endl;

  }else if(opt==4){ // remove PS
    cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" RemovePS-"<<cRemovePS<<" NotMatch-"<<cNotMatch<<".."<<endl;

  }

  return;
}


// loop opt=0: generate initial PSC; opt=1: findmaxdev; opt=2: dividepsc
void TreeReaderPulse::GenerateHCsLoop(int opt, int iconfig, int iChain, AGATA *agata, long long nentries){
  
  long long PSCstat[10];
  time(&start);

  long long istart = 0;
  int run;
  for(; irun<=fConfigs[iconfig].MaxRun; ){ //loop runs
    if(kInterrupt) break;

    {
#ifdef NTHREADS
      lock_guard<mutex> lock(treemtx); // lock tree read
#endif
      if(irun>fConfigs[iconfig].MaxRun) return;
      run = irun;
      irun++;
    }

    // initial tree
    fChain[iChain]->Reset();
    fChain[iChain]->AddFile(Form("%s/G4SimData%04d.root",fConfigs[iconfig].path.c_str(),run),0,"tree");
    int nentrytmp = fChain[iChain]->GetEntriesFast();
    Init(iChain);

    
    for(int ientry=0; ientry<nentrytmp; ientry++){ //loop evts
      if(kInterrupt) break;
    
      // output state
      if(ievt%100000==0 && kcout){
	kcout = false;      
	time(&stop);

#ifdef NTHREADS
	lock_guard<mutex> lock(treemtx); // lock tree read
#endif
	agata->GetPSCstat(PSCstat);
	if(PSCstat[1]>2000000){
	  //agata->SetAddNewPSC(false);
	  cout<<endl<<"PSC-"<<PSCstat[1]<<" : remove small PSC"<<endl;
	  agata->RemoveSmallPSC(5);
	}

	double MemUsageGB = GetCurrentMemoryUsage()/GB;
	double MemTotalGB = GetTotalSystemMemory()/GB;
	double MemUsage = MemUsageGB / MemTotalGB * 100;

	if(opt==0){ // initial PSC
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" maxnhits-"<<PSCstat[4]<<".."<<flush;

	}else if(opt==1){ // find max absdev
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" maxnhits-"<<PSCstat[4]<<".."
	      <<"MaxDev-"<<Form("%.2f",(float)MaxDev)<<".."<<flush;

	}else if(opt==2){ // divide PSC
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" DivPS-"<<cDivPS<<" NotMatch-"<<cNotMatch<<".."
	      <<" maxnhitsdiv-"<<maxnhitsdiv
	      <<" maxndiv-"<<PSCstat[8]<<".."<<flush;

	}else if(opt==3){ // find dev sigma
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" maxnhits-"<<PSCstat[4]<<".."<<flush;

	}else if(opt==4){ // remove PS
	  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" RemovePS-"<<cRemovePS<<" NotMatch-"<<cNotMatch<<".."<<flush;

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
      if     (opt==0) GenerateHCsworker(iconfig, run, iChain, agata, ientry, nentries);
      else if(opt==1) FindDevworker( 0, iconfig, run, iChain, agata, ientry, nentries, istart); // find absdev
      else if(opt==2) UpdateHCsworker( 0, iconfig, run, iChain, agata, ientry, nentries, istart); // divide HCs
      else if(opt==3) FindDevworker( 1, iconfig, run, iChain, agata, ientry, nentries, istart); // find dev
      else if(opt==4) UpdateHCsworker( 1, iconfig, run, iChain, agata, ientry, nentries, istart); // remove strange PS from PSCs

    }//end of loop evts

  }//end of loop runs
  
  return;
}


void TreeReaderPulse::GenerateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
					int ientry, long long nentries){

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

    if(SkipDet[detid]) continue;
    if(Detid>-1 && detid!=Detid) continue; // get PS for selected Detid

    int tmpnidx = -1, tmpnidxshift = 0;
#ifdef NOISE
    tmpnidx = (int)gRandom->Uniform(0,NOISE);
    tmpnidxshift = (int)gRandom->Uniform(0,NOISE/(BSIZE*NCHAN));
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
	if(aps.energy<PSCEMIN) continue;
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
      if(aps.energy<PSCEMIN) continue;
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

      if(SkipDet[detid]) continue;
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
	  if(aps.energy<PSCEMIN) continue;
	  fPS.push_back(aps);
	  fSegIdx.push_back(idx);
	}
#endif
	
      }else{ // one segment fired
        if(aps.det<0) continue;
#ifdef SINGLEHIT
	if(aps.nhits>1) continue;
#endif
	if(aps.energy<PSCEMIN) continue;
	fPS.push_back(aps);
	fSegIdx.push_back(segidx);

      }

    }//end of loop dets
  }

  if(fPS.size()<2) return; // at least one Compton scattering

#ifdef DIFFTOTE
  float Etot = 0;
  for(int i=0; i<fPS.size(); i++) Etot += fPS[i].energy;
#endif

  // create Hits-----------------------------------------------
  EventHits* fEvent = new EventHits(SourceE, SourcePos);
  fEvent->SetIdx(iconfig,run,ientry);

#ifdef DIFFTOTE
  fEvent->Etot = Etot;
#endif
  
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    TVector3 hitpos(fPS[i].labpos[0],fPS[i].labpos[1],fPS[i].labpos[2]);
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
  }

  long long iEvtHit = agata->AddEventHits(fEvent);
  vector<Hit*>* fHits = fEvent->GetfHits();
  vector<int> sign = fEvent->GetSign();

  // analysis clusters in tracking
  int Nunsigned = sign.size();
  int iclust = 0;
  while(Nunsigned>1){
    vector<Hit*>* tHits = new vector<Hit*>();
    vector<int> hitid;
    // add unused hits
    for(int i=0; i<sign.size(); i++){
      if(sign[i]<0){
	tHits->push_back( fHits->at(i) );
	hitid.push_back(i);
      }
    }
    if(tHits->size()<2){ delete tHits;  break;}

    int nsource = SourceE.size();
    vector<int> atrack;
    int bestis = 0;
    double minchi2 = 1e9;
    for(int is=0; is<nsource; is++){ // try different sources
      Tracker tracker(tHits, SourceE[is], SourcePos[is]);

#ifdef DIFFTOTE
#if    DIFFTOTE < 0 // DIFFTOTE<0 : only accept TotE match events
	if( !(fabs(Etot-SourceE[is])<fabs(DIFFTOTE)) ) continue;
#else   // DIFFTOTE>0 : if TotE match, put all hits in one clust
	if( fabs(Etot-SourceE[is])<fabs(DIFFTOTE) ) tracker.SetOneClust(true);
#endif
#endif

#ifdef ONECLUST
      tracker.SetOneClust(true);
#endif

      tracker.OFTtracking();
      //tracker.Simpletracking();
      double tmpchi2 = tracker.CalcChi2();
      if(tmpchi2>0 && tmpchi2<minchi2){
	bestis = is;
	minchi2 = tmpchi2;
	atrack = tracker.GetTrack();
      }
    }
    if(atrack.size()<2){ delete tHits; break;} // cannot find good track

    fEvent->SetBestis( iclust, bestis);
    for(int i=0; i<atrack.size(); i++){
      int hid = hitid[atrack[i]];
      sign[hid] = iclust;
      fEvent->SignClust( iclust, hid);
      Nunsigned--;
    }
    iclust++;
    delete tHits;
  }

  
  // group PS-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    if(Detid>-1 && fPS[i].det!=Detid) continue; // one det mode
    if(fSegIdx[i]>-1) continue; // multi segment fired
    if(sign[i]<0) continue;

    int tmp = agata->AddPS(&fPS[i], fHits->at(i)); // add to pulse shape collection
  }//end of loop fPS
  
  return;

}


void TreeReaderPulse::FindDevworker(int opt, int iconfig, int run, int iChain, AGATA *agata,
				    int ientry, long long nentries, long long &istart){

  {
#ifdef NTHREADS
    lock_guard<mutex> lock(treemtx); // lock tree read
#endif
    if(ievt>=nentries) return;
    if(ientry>=nentries) return;  

    ievt++;
  }
  fChain[iChain]->GetEntry(ientry);

  bool kFindDet = false;
  for(int idet=0; idet<obj[iChain].pdet->size(); idet++){
    if(obj[iChain].pdet->at(idet)==Detid) kFindDet = true;
  }

  if(Detid<0) kFindDet = true; // all det selected

  if(!kFindDet) return;

  // find Hit-----------------------------------------------
  long long iEvtHit = agata->FindiEvtHit(iconfig, run, ientry, istart);
  if(iEvtHit<0 || iEvtHit>NEventHits-1){  cievthitnotfind++; return;}
  cievthitfind++;

  vector<Hit*>* fHits = agata->FindEventHits(iEvtHit)->GetfHits();
  istart = iEvtHit;

  bool kFindPS = false;
  // check if need to compare hits with HCs
  for(int i=0; i<fHits->size(); i++){

    if(Detid>-1 && fHits->at(i)->GetDet()!=Detid) continue; // selected Detid

    kFindPS = true;
    if(kFindPS) break;
  }
  
  if(!kFindPS){ return;}

  // find corresponding fPS
  int ihit = 0;
  vector<PS> fPS;
  for(int idet=0; idet<obj[iChain].pdet->size(); idet++){ //loop dets
    int detid = obj[iChain].pdet->at(idet);

    if(SkipDet[detid]) continue;
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
    if(aps.energy<PSCEMIN) continue;
    fPS.push_back(aps);
    ihit++;
  }//end of loop dets
  
  // check fPS match with fHits
  cNotMatch++;
  if(Detid<0 && fHits->size()!=fPS.size()) return;
  for(int i=0; i<fPS.size(); i++){
    if(fHits->at(i)->GetDet()!=fPS[i].det || fHits->at(i)->GetSeg()!=fPS[i].seg) return;
    if(fHits->at(i)->GetE()!=fPS[i].energy) return;

    TVector3 hitpos(fPS[i].labpos[0],fPS[i].labpos[1],fPS[i].labpos[2]);
    if(fHits->at(i)->GetRealPosition()!=hitpos) return;
  }
  cNotMatch--;

  // check PSCs-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    if(opt==0){
      float tmpdev = agata->FindMaxDev(&fPS[i], fHits->at(i)); // find max abs deviation
      if(tmpdev>MaxDev) MaxDev = tmpdev;    
    }

    if(opt==1){
      agata->FindDevSigma(&fPS[i], fHits->at(i)); // find deviation sigma
    }
    
  }//end of loop fPS

  return;

}


void TreeReaderPulse::UpdateHCsworker(int opt, int iconfig, int run, int iChain, AGATA *agata,
				      int ientry, long long nentries, long long &istart){

  {
#ifdef NTHREADS
    lock_guard<mutex> lock(treemtx); // lock tree read
#endif
    if(ievt>=nentries) return;
    if(ientry>=nentries) return;  

    ievt++;
  }
  fChain[iChain]->GetEntry(ientry);

  bool kFindDet = false;
  for(int idet=0; idet<obj[iChain].pdet->size(); idet++){
    if(obj[iChain].pdet->at(idet)==Detid) kFindDet = true;
  }

  if(Detid<0) kFindDet = true; // all det selected

  if(!kFindDet) return;
  
  // find Hit-----------------------------------------------
  long long iEvtHit = agata->FindiEvtHit(iconfig, run, ientry, istart);
  if(iEvtHit<0 || iEvtHit>NEventHits-1){  cievthitnotfind++; return;}
  cievthitfind++;

  EventHits *fEvent = agata->FindEventHits(iEvtHit);
  vector<Hit*>* fHits = fEvent->GetfHits();
  vector<int> sign = fEvent->GetSign();

#ifdef DIFFTOTE
  float Etot = fEvent->Etot;
#endif

  istart = iEvtHit;

#ifdef CHECKTRACK  // check track----------------------------------------------
  int Nunsigned = sign.size();
  for(int i=0; i<Nunsigned; i++) sign[i] = -1; // remove previous clust
  int iclust = 0;

  // analysis clusters in tracking
  while(Nunsigned>1){
    vector<Hit*>* tHits = new vector<Hit*>();
    vector<int> hitid;
    // add unsigned hits
    for(int i=0; i<sign.size(); i++){
      if(sign[i]<0){
        tHits->push_back( fHits->at(i) );
        hitid.push_back(i);
      }
    }
    if(tHits->size()<2){ delete tHits;  break;}

    int nsource = SourceE.size();
    vector<int> atrack;
    int bestis = 0;
    double minchi2 = 1e9;
    for(int is=0; is<nsource; is++){
      Tracker tracker(tHits, SourceE[is], SourcePos[is]);

#ifdef DIFFTOTE
#if    DIFFTOTE < 0 // DIFFTOTE<0 : only accept TotE match events
	if( !(fabs(Etot-SourceE[is])<fabs(DIFFTOTE)) ) continue;
#else   // DIFFTOTE>0 : if TotE match, put all hits in one clust
	if( fabs(Etot-SourceE[is])<fabs(DIFFTOTE) ) tracker.SetOneClust(true);
#endif
#endif

#ifdef ONECLUST
      tracker.SetOneClust(true);
#endif

      tracker.OFTtracking();
      //tracker.Simpletracking();
      double tmpchi2 = tracker.CalcChi2();
      if(tmpchi2>0 && tmpchi2<minchi2){
	bestis = is;
	minchi2 = tmpchi2;
	atrack = tracker.GetTrack();
      }
    }
    if(atrack.size()<2){ delete tHits; break;} // cannot find good track

    fEvent->SetBestis(iclust, bestis);
    for(int i=0; i<atrack.size(); i++){
      int hid = hitid[atrack[i]];
      sign[hid] = iclust;
      //fEvent->SignClust( iclust, hid);
      Nunsigned--;
    }
    iclust++;
    delete tHits;
  }
#endif
  
  bool kFindPS = false;
  // check if need to compare hits with HCs
  for(int i=0; i<fHits->size(); i++){

    if(Detid>-1 && fHits->at(i)->GetDet()!=Detid) continue; // selected Detid

    if(opt==0){
      vector<HitCollection*>* hcs = fHits->at(i)->GetHitCollections();
      for(HitCollection* ahc : *hcs){
	if(ahc->GetSize() > MAXHITS) kFindPS = true;
	if(kFindPS) break;
      }
    }

    if(opt==1){
      kFindPS = true;
    }
    
    if(kFindPS) break;
  }
  
  if(!kFindPS){ return;}

  // find corresponding fPS
  int ihit = 0;
  vector<PS> fPS;
  for(int idet=0; idet<obj[iChain].pdet->size(); idet++){ //loop dets
    int detid = obj[iChain].pdet->at(idet);

    if(SkipDet[detid]) continue;
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
    if(aps.energy<PSCEMIN) continue;
    fPS.push_back(aps);
    ihit++;
  }//end of loop dets
  
  // check fPS match with fHits
  cNotMatch++;
  if(Detid<0 && fHits->size()!=fPS.size()) return;
  for(int i=0; i<fPS.size(); i++){
    if(fHits->at(i)->GetDet()!=fPS[i].det || fHits->at(i)->GetSeg()!=fPS[i].seg) return;
    if(fHits->at(i)->GetE()!=fPS[i].energy) return;

    TVector3 hitpos(fPS[i].labpos[0],fPS[i].labpos[1],fPS[i].labpos[2]);
    if(fHits->at(i)->GetRealPosition()!=hitpos) return;
  }
  cNotMatch--;

  // Update PSCs-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    if(opt==0){
      if(sign[i]<0) continue;
      int tmp = agata->AddPStoDiv(&fPS[i], fHits->at(i)); // add to divided pulse shape collection
      if(tmp>maxnhitsdiv) maxnhitsdiv = tmp;
      if(tmp>0) cDivPS++;
    }

    if(opt==1){
      if(sign[i]<0) continue;
      int tmp = agata->CheckPSinPSC(&fPS[i], fHits->at(i)); // check if PS within 3 sigma of PSC
      if(tmp>0) cRemovePS++;
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

  vector<vector<float>> simlabpos;
  vector<vector<float>> simdetpos;

#ifdef ADDPS
  TMatrixD              simspulse(BSIZE*NCHAN,1);
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

    vector<float> tmplabpos;
    vector<float> tmpdetpos;
    for(int ii=0; ii<3; ii++){
      tmplabpos.push_back(obj[iChain].posa->at(interid)[ii]); // lab position
      tmpdetpos.push_back(obj[iChain].posr->at(interid)[ii]); // det position
    }
    simlabpos.push_back(tmplabpos);
    simdetpos.push_back(tmpdetpos);

#ifdef ADDPS
    if(!skipPS){
      int itype = obj[iChain].pdet->at(idet)%3;
      TMatrixD tmppos(3,1);
      for(int ix=0; ix<3; ix++) tmppos(ix,0) = tmpdetpos[ix];
      double tmpenergy = tmphiteng[0];
      int tmpseg;
      TMatrixD tmpspulse(BSIZE*NCHAN,1);
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

  for(int ix=0; ix<3; ix++){
    aps.labpos[ix] = simlabpos[idx][ix];
    aps.detpos[ix] = simdetpos[idx][ix];
  }

  if(skipPS || segidx>-1) return aps; // skip PS
  
#ifdef ADDPS
  for(int iseg=0; iseg<NCHAN; iseg++){
    for(int isig=0; isig<BSIZE; isig++){
      double tmpamp = simspulse(iseg*BSIZE+isig,0);
      aps.opulse[iseg][isig] = tmpamp/simeng[0];
    }
  }
  
#else
  // read segments
  for(int iseg=0; iseg<NSEGS; iseg++){
    for(int isig=0; isig<BSIZE; isig++){
      double tmpamp = obj[iChain].spulse->at(idet)[iseg*BSIZE+isig];
      aps.opulse[iseg][isig] = tmpamp/simeng[0];
    }
  }

  // read core
  for(int isig=0; isig<BSIZE; isig++){
    double tmpamp = obj[iChain].core->at(idet)[isig];
    aps.opulse[NCHAN-1][isig] = tmpamp/simeng[0];
  }
#endif

#ifdef NOISE
  if(nidx>=0){
    // add noise
    for(int iseg=0; iseg<NCHAN; iseg++){
      for(int isig=0; isig<BSIZE; isig++){
	nidx = nidx%NOISE;
	aps.opulse[iseg][isig] += noise[nidx]/simeng[0];
	nidx++;
      }
      nidx+=nidxshift;
    }
  }
#endif

  // partial pulse shape for comparison
  int fseg[NCOMP];
  agata->GetGeo()->GetNextSegs(aps.seg, fseg);

  for(int is=0; is<NCOMP; is++){
    int iseg = fseg[is];
    copy_n(aps.opulse[iseg], BSIZE, aps.cpulse[is]);
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
