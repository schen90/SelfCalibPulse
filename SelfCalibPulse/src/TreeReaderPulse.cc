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
      int nsource = 0;
      vector<float>    tmpSE;
      vector<TVector3> tmpSPos;
      fin >> buffer >> nsource;
      for(int i=0; i<nsource; i++){
	fin >> buffer >> sE >> spos[0] >> spos[1] >> spos[2];
	tmpSE.push_back(sE);
	tmpSPos.push_back(TVector3(spos[0],spos[1],spos[2]));
      }
      fin >> buffer >> pathtmp;
      fin >> buffer >> run[0] >> run[1];
      fin >> buffer >> nevt;

      NSource.push_back(nsource);
      fSourceE.push_back(tmpSE);
      fSourcePos.push_back(tmpSPos);
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
    cout<<"#input "<<i<<":  "<<NSource[i]<<" sources"<<endl;
    for(int is=0; is<NSource[i]; is++){
      cout<<Form("     source %.1fkeV at %.2f %.2f %.2f",
		 fSourceE[i][is],
		 fSourcePos[i][is].X(),fSourcePos[i][is].Y(),fSourcePos[i][is].Z())<<endl;
    }

    for(int run=MinRun[i]; run<=MaxRun[i]; run++){
      fChain[0]->AddFile(Form("%s/Tree_%04d.root",path[i].c_str(),run),0,"TreeMaster");
    }

    long long nentries = fChain[0]->GetEntriesFast();
    cout<<" find \e[1m"<<nentries<<"\e[0m entries from rootfiles "
	<<Form("%s/Tree_%04d ~ %04d",path[i].c_str(),MinRun[i],MaxRun[i])<<endl;

    fChain[0]->Reset();
  }

  // init fChain[0] for  #input 0
  for(int run=MinRun[0]; run<=MaxRun[0]; run++){
    fChain[0]->AddFile(Form("%s/Tree_%04d.root",path[0].c_str(),run),0,"TreeMaster");
  }
  Init(0);

  SourceE   = fSourceE[0];
  SourcePos = fSourcePos[0];
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

  return;
}


void TreeReaderPulse::Init(int i){
  if (i<0 || i>=NChain) return;

  // set branch addresses and branch pointers
  fChain[i]->SetBranchAddress("EntryID",   &obj[i].EntryID);
  fChain[i]->SetBranchAddress("SegTraces",  obj[i].SegTraces);
  fChain[i]->SetBranchAddress("CoreTraces", obj[i].CoreTraces);
  fChain[i]->SetBranchAddress("SegE",       obj[i].SegE);
  fChain[i]->SetBranchAddress("CoreE",      obj[i].CoreE);
  fChain[i]->SetBranchAddress("CoreT",      obj[i].CoreT);
  fChain[i]->SetBranchAddress("CrystalId", &obj[i].CrystalId);
  fChain[i]->SetBranchAddress("CrystalTS", &obj[i].CrystalTS);

  return;
}


// loop opt=0: generate initial PSC; opt=1: findmaxdev; opt=2: dividepsc
void TreeReaderPulse::GenerateHCs(int opt, AGATA *agata){
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
    SourceE = fSourceE[i];
    SourcePos = fSourcePos[i];
    GenerateHCs(opt, agata, Nevts[i], i);
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
  long long nevents = 0;
  for(int run=MinRun[iconfig]; run<=MaxRun[iconfig]; run++){
    fChain[0]->Reset();
    fChain[0]->AddFile(Form("%s/Tree_%04d.root",path[iconfig].c_str(),run),0,"TreeMaster");
    Init(0);
    long tmpnentries = fChain[0]->GetEntriesFast();
    fChain[0]->GetEntry(tmpnentries-1);
    nevents += obj[0].EntryID+1;
  }
  if(nevts>0) nevents = TMath::Min(nevents,nevts);

  fChain[0]->Reset();
  for(int run=MinRun[iconfig]; run<=MaxRun[iconfig]; run++){
    fChain[0]->AddFile(Form("%s/Tree_%04d.root",path[iconfig].c_str(),run),0,"TreeMaster");
  }
  Init(0);

  long long nentries = fChain[0]->GetEntriesFast();
  irun = MinRun[iconfig];
  ievt = 0;

  cout<<"\e[1;33m Read "<<nentries<<" entries ( "<<nevents<<" events ) from rootfiles "
      <<Form("%s/Tree_%04d ~ %04d",path[iconfig].c_str(),MinRun[iconfig],MaxRun[iconfig])<<" ... \e[0m"<<endl;

  // Loop entries-------------------------------------------
#ifndef NTHREADS
  GenerateHCsLoop(opt, iconfig, 0, agata, nevents);

#else
  // loop trees with multi threads
  thread th[NTHREADS];
  cout<<"using "<<NTHREADS<<" threads:"<<endl;
  
  for(int i=0; i<NTHREADS; i++){
    th[i] = thread(&TreeReaderPulse::GenerateHCsLoop, this, opt, iconfig, i, ref(agata),nevents);
  }

  for(int i=0; i<NTHREADS; i++){
    if(th[i].joinable())
      th[i].join();
  }

#endif

  // output final statistics
  NEventHits = agata->GetEventHitsSize();
  
  time(&stop);
  long long PSCstat[10];
  agata->GetPSCstat(PSCstat);
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;


  if(opt==0){ // initial PSC
    cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."
	<<"fEventHits-"<<NEventHits<<".."<<endl;
    
  }else if(opt==1){ // find max dev
    cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."
	<<"MaxDev-"<<Form("%.2f",(float)MaxDev)<<".."<<endl;

  }else if(opt==2){ // divide PSC
    cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" DivPS-"<<cDivPS<<" NotMatch-"<<cNotMatch<<".."
	<<" maxnhitsdiv-"<<maxnhitsdiv
	<<" maxndiv-"<<PSCstat[8]<<".."<<endl;

  }else if(opt==3){ // find dev sigma
    cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" maxnhits-"<<PSCstat[4]<<".."<<endl;

  }else if(opt==4){ // remove PS
    cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	<<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"PS-"<<PSCstat[0]
	<<" PSC-"<<PSCstat[1]
	<<" RemovePS-"<<cRemovePS<<" NotMatch-"<<cNotMatch<<".."<<endl;

  }

  return;
}


// loop opt=0: generate initial PSC; opt=1: findmaxdev; opt=2: dividepsc
void TreeReaderPulse::GenerateHCsLoop(int opt, int iconfig, int iChain, AGATA *agata, long long nevents){
  
  long long PSCstat[10];
  time(&start);

  long long istart = 0;
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
    fChain[iChain]->AddFile(Form("%s/Tree_%04d.root",path[iconfig].c_str(),run),0,"TreeMaster");
    int nentrytmp = fChain[iChain]->GetEntriesFast();
    Init(iChain);

    
    for(int ientry=0; ientry<nentrytmp; ){ //loop evts
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
	  cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" maxnhits-"<<PSCstat[4]<<".."<<flush;

	}else if(opt==1){ // find max absdev
	  cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" maxnhits-"<<PSCstat[4]<<".."
	      <<"MaxDev-"<<Form("%.2f",(float)MaxDev)<<".."<<flush;

	}else if(opt==2){ // divide PSC
	  cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" DivPS-"<<cDivPS<<" NotMatch-"<<cNotMatch<<".."
	      <<" maxnhitsdiv-"<<maxnhitsdiv
	      <<" maxndiv-"<<PSCstat[8]<<".."<<flush;

	}else if(opt==3){ // find dev sigma
	  cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
	      <<Form("(%.0fs/10kevts)..",difftime(stop,start)/10.)
	      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<PSCstat[0]
	      <<" PSC-"<<PSCstat[1]
	      <<" maxnhits-"<<PSCstat[4]<<".."<<flush;

	}else if(opt==4){ // remove PS
	  cout<<"\r finish read "<<ievt<<" / "<<nevents<<" evts"
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
      int multi = 0;
      if     (opt==0) multi=GenerateHCsworker(  iconfig, run, iChain, agata, ientry, nentrytmp);
      else if(opt==1) multi=FindDevworker(   0, iconfig, run, iChain, agata, ientry, nentrytmp, istart); // find absdev
      else if(opt==2) multi=UpdateHCsworker( 0, iconfig, run, iChain, agata, ientry, nentrytmp, istart); // divide HCs
      else if(opt==3) multi=FindDevworker(   1, iconfig, run, iChain, agata, ientry, nentrytmp, istart); // find dev
      else if(opt==4) multi=UpdateHCsworker( 1, iconfig, run, iChain, agata, ientry, nentrytmp, istart); // remove strange PS from PSCs

      {
#ifdef NTHREADS
	lock_guard<mutex> lock(treemtx); // lock tree read
#endif
	if(ievt>=nevents) return;
	ientry += multi;
	ievt++;
	if(ievt>=nevents) return;
      }

    }//end of loop evts

  }//end of loop runs
  
  return;
}


int TreeReaderPulse::GenerateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
				       int ientry, long long nentries){

  // get PS in one event-----------------------------------
  vector<PS> fPS;
  vector<int> fSegIdx;

  fChain[iChain]->GetEntry(ientry);
  int multi = 0;
  int EvtID = obj[iChain].EntryID;

  while( EvtID == obj[iChain].EntryID ){ // find det with same EntryID
    multi++;

    int detid = obj[iChain].CrystalId;

    // get PS for selected Detid
    if(Detid>-1 && detid!=Detid){
      if(ientry+multi >= nentries) break;
      fChain[iChain]->GetEntry(ientry + multi); // get next entry
      continue;
    }

    int segidx = -1;
    PS aps = GetAPS(iChain, false, segidx); // aps w/ PS

    if(aps.det<0 && segidx>-1){ // multi-segment fired
#ifdef MULTISEG
      for(int idx=0; idx<segidx; idx++){
	PS aps = GetAPS(iChain, true, idx); // aps w/o PS
        if(aps.det<0) continue;
	fPS.push_back(aps);
	fSegIdx.push_back(idx);
      }
#endif
      
    }else{ // one segment fired
      if(aps.det<0) continue;
      fPS.push_back(aps);
      fSegIdx.push_back(segidx);
    }

    if(ientry+multi >= nentries) break;
    fChain[iChain]->GetEntry(ientry + multi); // get next entry
  }//end of loop dets

  
  if(fPS.size()==0) return multi; // skip if selected Det not fired

  if(Detid>-1){ // one det mode
    // get other PSs in one event-----------------------------------
    for(int idet=0; idet<multi; idet++){ //loop dets
      fChain[iChain]->GetEntry(ientry + idet);

      int detid = obj[iChain].CrystalId;

      if(Detid>-1 && detid==Detid) continue; // get PS for other Dets

      int segidx = -1;
      PS aps = GetAPS(iChain, true, segidx); // aps w/o PS

      if(aps.det<0 && segidx>-1){ // multi-segment fired
#ifdef MULTISEG
	for(int idx=0; idx<segidx; idx++){
	  PS aps = GetAPS(iChain, true, idx); // aps w/o PS
	  if(aps.det<0) continue;
	  fPS.push_back(aps);
	  fSegIdx.push_back(idx);
	}
#endif
	
      }else{ // one segment fired
        if(aps.det<0) continue;
	fPS.push_back(aps);
	fSegIdx.push_back(segidx);

      }

    }//end of loop dets
  }

  if(fPS.size()<2) return multi; // at least one Compton scattering

#ifdef DIFFTOTE
  float Etot = 0;
  for(int i=0; i<fPS.size(); i++) Etot += fPS[i].energy;
  bool kskip = true;
  for(int is=0; is<SourceE.size(); is++) if(fabs(Etot-SourceE[is])<DIFFTOTE) kskip = false;
  if(kskip) return multi;
#endif

  // create Hit-----------------------------------------------
  vector<int> uflag;
  EventHits* fEvent = new EventHits(SourceE, SourcePos);
  fEvent->SetIdx(iconfig,run,ientry,EvtID);

#ifdef DIFFTOTE
  fEvent->Etot = Etot;
#endif

  for(int i=0; i<fPS.size(); i++){ //loop fPS

    TVector3 initpos;
    if(Detid<0 || fPS[i].det==Detid){
      initpos = agata->GetPSpos(fPS[i].det, fPS[i].seg, &fPS[i]);
    }else{
      TMatrixD SegPos = agata->GetGeo()->GetSegPos(fPS[i].det,fPS[i].seg);
      initpos.SetXYZ(SegPos(0,0), SegPos(1,0), SegPos(2,0));
    }
    
    Hit *ahit = new Hit(fPS[i].det, fPS[i].seg, fPS[i].energy, initpos); //keV and mm
    
    fEvent->Add(ahit);
    uflag.push_back(1);
  }

  long long iEvtHit = agata->AddEventHits(fEvent);
  vector<Hit*>* fHits = fEvent->GetfHits();
  
  int nsource = SourceE.size();
  vector<int> atrack;
  int bestis = 0;
  double minchi2 = 1e9;
  for(int is=0; is<nsource; is++){
#ifdef DIFFTOTE
    if( !(fabs(Etot-SourceE[is])<DIFFTOTE) ) continue;
#endif
    Tracker tracker(fHits, SourceE[is], SourcePos[is]);
    tracker.OFTtracking();
    //tracker.Simpletracking();
    double tmpchi2 = tracker.CalcChi2();
    if(tmpchi2>0 && tmpchi2<minchi2){
      bestis = is;
      minchi2 = tmpchi2;
      atrack = tracker.GetTrack();
    }
  }
  fEvent->SetBestis(bestis);
  
#ifdef CHECKTRACK
  // check track----------------------------------------------
  for(int i=0; i<fPS.size(); i++) uflag[i] = 0;
  if(atrack.size()>1) for(int i=0; i<atrack.size(); i++) uflag[atrack[i]] = 1;
#endif
  
  // group PS-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    if(Detid>-1 && fPS[i].det!=Detid) continue; // one det mode
    if(fSegIdx[i]>-1) continue; // multi segment fired
    if(uflag[i]!=1) continue;

    int tmp = agata->AddPS(&fPS[i], fHits->at(i)); // add to pulse shape collection
  }//end of loop fPS
  
  return multi;

}


int TreeReaderPulse::FindDevworker(int opt, int iconfig, int run, int iChain, AGATA *agata,
				   int ientry, long long nentries, long long &istart){

  // find Hit-----------------------------------------------
  long long iEvtHit = agata->FindiEvtHit(iconfig, run, ientry, istart);
  if(iEvtHit<0 || iEvtHit>NEventHits-1){ ievt++; return 1;}

  vector<Hit*>* fHits = agata->FindEventHits(iEvtHit)->GetfHits();
  istart = iEvtHit;

  bool kFindPS = false;
  // check if need to compare hits with HCs
  for(int i=0; i<fHits->size(); i++){

    if(Detid>-1 && fHits->at(i)->GetDet()!=Detid) continue; // selected Detid

    kFindPS = true;
    if(kFindPS) break;
  }
  
  if(!kFindPS){ ievt++; return fHits->size();}

  // find corresponding fPS
  vector<PS> fPS;

  fChain[iChain]->GetEntry(ientry);
  int multi = 0;
  int EvtID = obj[iChain].EntryID;

  int ihit = 0;
  while( EvtID == obj[iChain].EntryID ){ // find det with same EntryID
    multi++;

    int detid = obj[iChain].CrystalId;

    // get PS for selected Detid
    if(Detid>-1 && detid!=Detid){
      if(ientry+multi >= nentries) break;
      fChain[iChain]->GetEntry(ientry + multi); // get next entry
      continue;
    }

    int segidx = -1;
    PS aps = GetAPS(iChain, false, segidx); // aps w/ PS

    if(aps.det<0) return fHits->size();
    fPS.push_back(aps);
    ihit++;

    if(ientry+multi >= nentries) break;
    fChain[iChain]->GetEntry(ientry + multi); // get next entry
  }//end of loop dets
  
  // check fPS match with fHits
  cNotMatch++;
  if(Detid<0 && fHits->size()!=fPS.size()) return multi;
  for(int i=0; i<fPS.size(); i++){
    if(fHits->at(i)->GetDet()!=fPS[i].det || fHits->at(i)->GetSeg()!=fPS[i].seg) return multi;
    if(fHits->at(i)->GetE()!=fPS[i].energy) return multi;
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

  return multi;

}


int TreeReaderPulse::UpdateHCsworker(int opt, int iconfig, int run, int iChain, AGATA *agata,
				     int ientry, long long nentries, long long &istart){

  // find Hit-----------------------------------------------
  long long iEvtHit = agata->FindiEvtHit(iconfig, run, ientry, istart);
  if(iEvtHit<0 || iEvtHit>NEventHits-1){ ievt++; return 1;}

  vector<Hit*>* fHits = agata->FindEventHits(iEvtHit)->GetfHits();
#ifdef DIFFTOTE
  float Etot = agata->FindEventHits(iEvtHit)->Etot;
#endif
  vector<int> uflag;
  for(int i=0; i<fHits->size(); i++){ //loop fPS
    uflag.push_back(1);
  }
  istart = iEvtHit;

#ifdef CHECKTRACK
  // check track----------------------------------------------
  int nsource = SourceE.size();
  vector<int> atrack;
  int bestis = 0;
  double minchi2 = 1e9;
  for(int is=0; is<nsource; is++){
#ifdef DIFFTOTE
    if( !(fabs(Etot-SourceE[is])<DIFFTOTE) ) continue;
#endif
    Tracker tracker(fHits, SourceE[is], SourcePos[is]);
    tracker.OFTtracking();
    //tracker.Simpletracking();
    double tmpchi2 = tracker.CalcChi2();
    if(tmpchi2>0 && tmpchi2<minchi2){
      bestis = is;
      minchi2 = tmpchi2;
      atrack = tracker.GetTrack();
    }
  }
  for(int i=0; i<fHits->size(); i++) uflag[i] = 0;
  if(atrack.size()>1) for(int i=0; i<atrack.size(); i++) uflag[atrack[i]] = 1;
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
  
  if(!kFindPS){ ievt++; return fHits->size();}

  // find corresponding fPS
  vector<PS> fPS;

  fChain[iChain]->GetEntry(ientry);
  int multi = 0;
  int EvtID = obj[iChain].EntryID;

  int ihit = 0;
  while( EvtID == obj[iChain].EntryID ){ // find det with same EntryID
    multi++;

    int detid = obj[iChain].CrystalId;

    // get PS for selected Detid
    if(Detid>-1 && detid!=Detid){
      if(ientry+multi >= nentries) break;
      fChain[iChain]->GetEntry(ientry + multi); // get next entry
      continue;
    }
    
    int segidx = -1;
    PS aps = GetAPS(iChain, false, segidx); // aps w/ PS

    if(aps.det<0) return fHits->size();
    fPS.push_back(aps);
    ihit++;

    if(ientry+multi >= nentries) break;
    fChain[iChain]->GetEntry(ientry + multi); // get next entry
  }//end of loop dets
  
  // check fPS match with fHits
  cNotMatch++;
  if(Detid<0 && fHits->size()!=fPS.size()) return multi;
  for(int i=0; i<fPS.size(); i++){
    if(fHits->at(i)->GetDet()!=fPS[i].det || fHits->at(i)->GetSeg()!=fPS[i].seg) return multi;
    if(fHits->at(i)->GetE()!=fPS[i].energy) return multi;
  }
  cNotMatch--;

  // Update PSCs-------------------------------------------------
  for(int i=0; i<fPS.size(); i++){ //loop fPS

    if(opt==0){
      if(uflag[i]!=1) continue;
      int tmp = agata->AddPStoDiv(&fPS[i], fHits->at(i)); // add to divided pulse shape collection
      if(tmp>maxnhitsdiv) maxnhitsdiv = tmp;
      if(tmp>0) cDivPS++;
    }

    if(opt==1){
      if(uflag[i]!=1) continue;
      int tmp = agata->CheckPSinPSC(&fPS[i], fHits->at(i)); // check if PS within 3 sigma of PSC
      if(tmp>0) cRemovePS++;
    }
    
  }//end of loop fPS

  return multi;

}


PS TreeReaderPulse::GetAPS(int iChain, bool skipPS, int &segidx){
  PS aps;
  aps.det = -1;
  vector<int>           vseg;
  vector<float>         veng;
  for(int iseg=0; iseg<NSeg; iseg++){
    if(obj[iChain].SegE[iseg]>10){
      vseg.push_back(iseg);
      veng.push_back(obj[iChain].SegE[iseg]);
    }
  }

  if(segidx<0 && veng.size()!=1){
    segidx = veng.size(); // multi-segment fired
    return aps;// require only 1 seg fired in a det
  }
  if(segidx<0 && veng.size()==1){
    if( !(fabs(veng[0]-obj[iChain].CoreE[0])<3) ){
      segidx = veng.size(); // multi-segment hit below threshold
      return aps;// require only 1 seg fired in a det
    }
  }

  int idx;
  if(segidx<0) idx = 0;
  else         idx = segidx;

  aps.det = obj[iChain].CrystalId;
  aps.seg = vseg[idx];
  aps.energy = veng[idx];

  if(skipPS || segidx>-1) return aps; // skip PS


  // read segments
  for(int iseg=0; iseg<NSeg; iseg++){
    for(int isig=0; isig<NSig; isig++){
      double tmpamp = obj[iChain].SegTraces[iseg*NSig+isig];
      aps.opulse[iseg][isig] = tmpamp/veng[0];
    }
  }

  // read core
  for(int isig=0; isig<NSig; isig++){
    double tmpamp = obj[iChain].CoreTraces[isig];
    aps.opulse[NSegCore-1][isig] = tmpamp/veng[0];
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
