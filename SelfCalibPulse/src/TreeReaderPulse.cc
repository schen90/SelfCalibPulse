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

using namespace std;

TreeReaderPulse::TreeReaderPulse(){

#ifdef ADDPS
  apsb = new PSbasis();
#endif

#ifdef NOISE
  LoadNoise();
#endif

  for(int i=0; i<NChain; i++){
    fChain[i] = new TChain();
  }

  ievt = 0;
  cPaths = 0;
  iter = 0;
  kcout = true;
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

  for(int i=0; i<nConfig; i++){
    gROOT->ProcessLine(Form(".!mkdir ./share/Hits/input%d",i));
    gROOT->ProcessLine(Form(".!mkdir ./share/Hits/input%d/tmp",i));
    gROOT->ProcessLine(Form(".!mkdir ./share/Hits/input%d/it",i));
  }

#ifdef NOISE
  MakeNoise();
#endif
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

  fChain[i]->SetBranchAddress("posa",&obj[i].posa);
  fChain[i]->SetBranchAddress("posr",&obj[i].posr);

  fChain[i]->SetBranchAddress("pdet",&obj[i].pdet);
  fChain[i]->SetBranchAddress("ecore",&obj[i].ecore);
  fChain[i]->SetBranchAddress("inter",&obj[i].inter);

  fChain[i]->SetBranchAddress("pseg",&obj[i].pseg);
  fChain[i]->SetBranchAddress("ngrid",&obj[i].ngrid);

#ifndef ADDPS
  fChain[i]->SetBranchAddress("core",&obj[i].core);
  fChain[i]->SetBranchAddress("spulse",&obj[i].spulse);
#endif

  fChain[i]->SetBranchAddress("category",&obj[i].category);

  return;
}


void TreeReaderPulse::InitEvtTree(int i){
  if (i<0 || i>=NChain) return;

  // set branch
  //EvtTree[i]->Branch("iconfig",&evtobj[i].iconfig);
  //EvtTree[i]->Branch("irun",&evtobj[i].irun);
  EvtTree[i]->Branch("ientry",&evtobj[i].ientry);
  EvtTree[i]->Branch("nhits",&evtobj[i].nhits);

  EvtTree[i]->Branch("interid",&evtobj[i].interid);
  EvtTree[i]->Branch("idet",&evtobj[i].idet);
  EvtTree[i]->Branch("iseg",&evtobj[i].iseg);
  EvtTree[i]->Branch("depE",&evtobj[i].depE);
  EvtTree[i]->Branch("labpos",evtobj[i].labpos,"labpos[3]/F");
  EvtTree[i]->Branch("calpos",evtobj[i].calpos,"calpos[3]/F");
  EvtTree[i]->Branch("dist",&evtobj[i].dist);
#ifdef NOISE
  EvtTree[i]->Branch("noiseidx",&evtobj[i].noiseidx);
  EvtTree[i]->Branch("noiseidxshift",&evtobj[i].noiseidxshift);
#endif
  EvtTree[i]->Branch("sourcepos",evtobj[i].sourcepos,"sourcepos[3]/F");
  EvtTree[i]->Branch("sourceeng",&evtobj[i].sourceeng);

  return;
}


void TreeReaderPulse::GeneratePSC(AGATA *agata){
  cout<<"\e[1;31m Generate PSC ... \e[0m"<<endl;

  cPaths = 0;
  if(nConfig<1){
    cerr<<"cannot find input, check configure file..."<<endl;
    return;
  }

  for(int i=0; i<nConfig; i++){
    cout<<"\e[1;32m"<<" #input "<<i<<"\e[0m"<<endl;
    SourceE = fSourceE[i];
    SourcePos = fSourcePos[i];
    GeneratePSC( agata, i, Nevts[i]);
  }

  iter++;
  return;
}


void TreeReaderPulse::GeneratePSC(AGATA *agata, int nevts){

  cPaths = 0;
  GeneratePSC( agata, 0, nevts);
  return;
}


void TreeReaderPulse::GeneratePSC(AGATA *agata, int iconfig, int nevts){

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
  GeneratePSCLoop( 0, agata, iconfig, nentries);

#else
  // loop trees with multi threads
  thread th[NTHREADS];
  cout<<"using "<<NTHREADS<<" threads:"<<endl;

  for(int i=0; i<NTHREADS; i++){
    th[i] = thread(&TreeReaderPulse::GeneratePSCLoop, this, i, ref(agata), iconfig, nentries);
  }

  for(int i=0; i<NTHREADS; i++){
    if(th[i].joinable())
      th[i].join();
  }

#endif

  // output final statistics
  time(&stop);
  int PSCstat[10];
  agata->GetPSCstat(PSCstat);
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;

  cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
      <<Form("(%.0fs/10kevts)..",difftime(stop,start))
      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
      <<" PSC-"<<PSCstat[0]
      <<" Path-"<<cPaths
      <<" maxnhits-"<<PSCstat[2]<<endl;
  
  return;
}


void TreeReaderPulse::GeneratePSCLoop(int iChain, AGATA *agata,
				      int iconfig, int nentries){

  int PSCstat[10];
  time(&start);

  int istart = 0;
  int run;
  for(; irun<=MaxRun[iconfig]; ){ //loop runs

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

    // write to EvtTree
    EvtFile[iChain] = new TFile(Form("./share/Hits/input%d/tmp/run%04d.root",iconfig,run),"RECREATE");
    EvtTree[iChain] = new TTree("tree","EventHits");
    InitEvtTree(iChain);
    
    for(int ientry=0; ientry<nentrytmp; ientry++){ //loop evts

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

	cout<<"\r finish read "<<ievt<<" / "<<nentries<<" evts"
	    <<Form("(%.0fs/10kevts)..",difftime(stop,start))
	    <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	    <<" PSC-"<<PSCstat[0]
	    <<" Path-"<<cPaths
	    <<" maxnhits-"<<PSCstat[2]<<flush;

	time(&start);
        kcout = true;
      }

      // process one event
      ProcessOneEvent( iChain, agata, iconfig, run, ientry, nentries);
      
    }//end of loop evts

    EvtFile[iChain]->cd();
    EvtTree[iChain]->Write();
    EvtFile[iChain]->Close();
    gROOT->ProcessLine(Form(".!mv -f ./share/Hits/input%d/tmp/run%04d.root ./share/Hits/input%d/run%04d.root",iconfig,run,iconfig,run));

    if(run==0){
      gROOT->ProcessLine(Form(".!cp -pdr ./share/Hits/input%d/run%04d.root ./share/Hits/input%d/it/run%04d_Fit%d.root",iconfig,run,iconfig,run,iter));
    }
  }//end of loop runs

  return;
}


void TreeReaderPulse::ProcessOneEvent(int iChain, AGATA *agata,
				      int iconfig, int run,
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

  if(fPS.size()<2) return; // at least one Compton scattering

  // create Hit-----------------------------------------------
  vector<int> uflag;
  EventHits* fEvent = new EventHits(SourceE, SourcePos);
  fEvent->SetIdx(iconfig,run,ientry);

  for(int i=0; i<fPS.size(); i++){ //loop fPS
    
    TVector3 hitpos(fPS[i].labpos[0],fPS[i].labpos[1],fPS[i].labpos[2]);
    TVector3 initpos = agata->GetPSpos(fPS[i].det, fPS[i].seg, &fPS[i]); // initpos from PSA

    Hit *ahit = new Hit(fPS[i].det, fPS[i].seg, fPS[i].energy, hitpos, initpos); //keV and mm
    ahit->SetInterid(fPS[i].interid);
#ifdef NOISE
    ahit->SetNoiseIdx(Nidx[i]);
    ahit->SetNoiseIdxShift(Nidxshift[i]);
#endif

    fEvent->Add(ahit);
    uflag.push_back(1);
  }

  vector<Hit*>* fHits = fEvent->GetfHits();

  // write EventHits -----------------------------------------
  evtobj[iChain].iconfig = iconfig;
  evtobj[iChain].irun = run;
  evtobj[iChain].ientry = ientry;
  evtobj[iChain].nhits = fHits->size();
  evtobj[iChain].sourcepos[0] = SourcePos.X();
  evtobj[iChain].sourcepos[1] = SourcePos.Y();
  evtobj[iChain].sourcepos[2] = SourcePos.Z();
  evtobj[iChain].sourceeng = SourceE;
  for(int i=0; i<fHits->size(); i++){
    evtobj[iChain].interid = i;
    evtobj[iChain].idet = fHits->at(i)->GetDet();
    evtobj[iChain].iseg = fHits->at(i)->GetSeg();
    evtobj[iChain].depE = fHits->at(i)->GetE();
    TVector3 LabPos = fHits->at(i)->GetRealPosition();
    TVector3 CalPos = fHits->at(i)->GetPosition();
    evtobj[iChain].labpos[0] = LabPos.X();  evtobj[iChain].calpos[0] = CalPos.X();
    evtobj[iChain].labpos[1] = LabPos.Y();  evtobj[iChain].calpos[1] = CalPos.Y();
    evtobj[iChain].labpos[2] = LabPos.Z();  evtobj[iChain].calpos[2] = CalPos.Z();
    evtobj[iChain].dist = (LabPos - CalPos).Mag();
#ifdef NOISE
    evtobj[iChain].noiseidx = fHits->at(i)->GetNoiseIdx();
    evtobj[iChain].noiseidxshift = fHits->at(i)->GetNoiseIdxShift();
#endif
    EvtTree[iChain]->Fill();
  }
  
  // tracking ----------------------------------------------
  Tracker tracker(fHits, (Double_t)SourceE, SourcePos);
  //tracker.OFTtracking();
  tracker.Simpletracking();
  vector<int> atrack = tracker.GetTrack();

  // fit paths
  if(atrack.size()>1){ // at least two hits
    double incE = SourceE;
    double depE = fHits->at(atrack[0])->GetE(); // keV

    Hit *sourcehit = new Hit(SourcePos);
    sourcehit->SetInterid(-1);

    // first scattering
    Path *apath = new Path(sourcehit,fHits->at(atrack[0]),fHits->at(atrack[1]),
			   incE, depE, incE, depE);
    int ipoint = agata->FitPath(apath);
    if(ipoint>-1){
      agata->AddPStoPSC(&fPS[atrack[0]], ipoint);
      cPaths++;
    }

    delete apath;
    delete sourcehit;
    
    // next scattering
    incE = incE - depE;
    for(int i=1; i<atrack.size()-1; i++){
      depE = fHits->at(atrack[i])->GetE(); // keV

      Path *apath = new Path(fHits->at(atrack[i-1]),fHits->at(atrack[i]),fHits->at(atrack[i+1]),
			     incE, depE, incE, depE);
      int ipoint = agata->FitPath(apath);
      if(ipoint>-1){
	agata->AddPStoPSC(&fPS[atrack[i]], ipoint);
	cPaths++;
      }

      delete apath;
      
      incE = incE - depE;
    }
  }

  delete fEvent;
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
