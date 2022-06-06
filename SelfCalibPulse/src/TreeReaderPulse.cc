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
	for(int idx=0; idx<segidx; idx++){
	  PS aps = GetAPS(iChain,agata,idet,-1,0,true,idx); // aps w/o PS
	  if(aps.det<0) continue;
#ifdef SINGLEHIT
	  if(aps.nhits>1) continue;
#endif
	  fPS.push_back(aps);
	  fSegIdx.push_back(idx);
	}

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
    uflag.push_back(0);
  }
  istart = iEvtHit;

  // check track----------------------------------------------
  Tracker tracker(fHits, (Double_t)SourceE, SourcePos);
  tracker.OFTtracking();
  vector<int> atrack = tracker.GetTrack();
  if(atrack.size()>1) for(int i=0; i<atrack.size(); i++) uflag[atrack[i]] = 1;

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
    PS aps = GetAPS(iChain,agata,idet,tmpnidx,tmpnidxshift);
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

    simseg.push_back(obj[iChain].pseg->at(idet)[i]-1); //pseg start from 1...
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
