#ifndef AGATA_CC
#define AGATA_CC

#include <fstream>
#include <iostream>
#include <algorithm>
#include <x86intrin.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>

#include "AGATA.hh"
#include "AGATAgeo.hh"

using namespace std;

AGATA::AGATA(int detid){
  Detid = detid; // selected det

  agatageo = new AGATAgeo();
  NDets = agatageo->GetNDets();

#ifdef PSA
  ReadPSAbasis();
#endif
  
  cPStotal  = 0;
  cPSCtotal = 0;
  cPSCmem   = 0;
  cPSCfile  = 0;
  maxnhits  = 0;
  cPaths    = 0;
  cHits     = 0;
  cHCs      = 0;

  for(int idet=0; idet<MaxNDets; idet++)
    for(int iseg=0; iseg<NSeg; iseg++)
      fHCs[idet][iseg] = new vector<HitCollection*>();
  fAllHCs = new vector<HitCollection*>();
  
  fEventHits = new vector<EventHits*>();
  fPaths = new vector<Path*>();

  // parameters for HC pos optimize
  fitlimit = 5;
  // PSC number limit
  for(int idet=0; idet<MaxNDets; idet++){

    if(Detid<0) PSClimit[idet] = 2500;
    else        PSClimit[idet] = 5000;
    if(idet==0) PSClimit[idet] = 5000;
    
    for(int iseg=0; iseg<NSeg; iseg++){
      kAddNewPSC[idet][iseg] = false;
    }
  }
}

AGATA::~AGATA(){
  ClearPSCMem();
}

//------------------------------------------------
// Pulse Shape Collections
//------------------------------------------------

void AGATA::WritePSCfiles(){
  for(int idet=0; idet<NDets; idet++){
    LoadPSCfiles(idet);
    WritePSCfiles(idet);
    cout<<endl;
  }
}

void AGATA::WritePSCfiles(int detid){ // create Pulse Shape Collection files

  //ClosePSCFiles();

  vector <int> idlist; // detid list to write
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  string pscfilesname0 = "./PSCfiles/tmp/Det";
  string pscfilesname = "./PSCfiles/Det";
  //if(Detid>-1) pscfilesname = "./share/PSCs/Det";
  cout<<"\e[1;33m Create PSCfiles for "<<pscfilesname;
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;


  cPSCfile = 0;
  for(int idet : idlist){

    pscfile[idet] = new TFile(Form("%s%04d.root",pscfilesname0.c_str(),idet),"RECREATE");

    for(int iseg=0; iseg<NSeg; iseg++){
      cout<<"\r writing tree for Det "<<idet<<" Seg "<<iseg<<" ..."<<flush;
      psctree[idet][iseg] = new TTree(Form("tree%d",iseg),Form("pulse shape tree for det%04d seg%d",idet,iseg));

      InitTreeWrite(psctree[idet][iseg]);

      for(int ic=0; ic<fHCs[idet][iseg]->size(); ic++){
	det = fHCs[idet][iseg]->at(ic)->GetDet();
	seg = fHCs[idet][iseg]->at(ic)->GetSeg();
	index = fHCs[idet][iseg]->at(ic)->GetPid();
	nhits = fHCs[idet][iseg]->at(ic)->GetSize();

	Marker = fHCs[idet][iseg]->at(ic)->Marker;
	for(int ix=0; ix<3; ix++) chi2limit[ix] = fHCs[idet][iseg]->at(ic)->MaxChi2s[ix];
	
	TVector3 vtmp;

	vtmp = fHCs[idet][iseg]->at(ic)->GetInitPosition();
	calpos[0] = vtmp.X();  calpos[1] = vtmp.Y();  calpos[2] = vtmp.Z();
	TMatrixD CalPos(3,1);
	for(int ix=0; ix<3; ix++) CalPos(ix,0) = calpos[ix];
	TMatrixD CadPos = agatageo->Lab2DetPos(det,CalPos);
	for(int ix=0; ix<3; ix++) cadpos[ix] = CadPos(ix,0);

	vtmp = fHCs[idet][iseg]->at(ic)->GetPosition();
	calpos2[0] = vtmp.X();  calpos2[1] = vtmp.Y();  calpos2[2] = vtmp.Z();
	TMatrixD CalPos2(3,1);
	for(int ix=0; ix<3; ix++) CalPos2(ix,0) = calpos2[ix];
	TMatrixD CadPos2 = agatageo->Lab2DetPos(det,CalPos2);
	for(int ix=0; ix<3; ix++) cadpos2[ix] = CadPos2(ix,0);

#ifdef REALPOS
	vtmp = fHCs[idet][iseg]->at(ic)->GetRealPosition();
	labpos[0] = vtmp.X();  labpos[1] = vtmp.Y();  labpos[2] = vtmp.Z();
	TMatrixD LabPos(3,1);
	for(int ix=0; ix<3; ix++) LabPos(ix,0) = labpos[ix];
	TMatrixD DetPos = agatageo->Lab2DetPos(det,LabPos);
	for(int ix=0; ix<3; ix++) detpos[ix] = DetPos(ix,0);

	dist = 0;
	for(int ix=0; ix<3; ix++) dist += SQ(calpos[ix]-labpos[ix]);
	dist = sqrt(dist);

	dist2 = 0;
	for(int ix=0; ix<3; ix++) dist2 += SQ(calpos2[ix]-labpos[ix]);
	dist2 = sqrt(dist2);

	if(kWithPS) copy_n(fPSC[idet][iseg][ic].cpos, 3, cpos);
#endif

#ifdef WITHPS
	if(kWithPS){
	  cpulsehits = fPSC[idet][iseg][ic].cpulsehits;
	  for(int iiseg=0; iiseg<NSegCore; iiseg++){
	    copy_n(fPSC[idet][iseg][ic].spulse[iiseg], NSig, spulse[iiseg]);
	  }
	  for(int iiseg=0; iiseg<NSeg_comp; iiseg++){
	    copy_n(fPSC[idet][iseg][ic].cpulse[iiseg], NSig, cpulse[iiseg]);
	  }
	  copy_n(fPSC[idet][iseg][ic].segwgt, NSeg_comp, segwgt);
	}
#endif
	npaths = fHCs[idet][iseg]->at(ic)->GetPaths()->size();
	
	psctree[idet][iseg]->Fill();
	cPSCfile++;
      }
      pscfile[idet]->cd();
      psctree[idet][iseg]->Write();
    }//end of loop seg

    pscfile[idet]->Close();
    gROOT->ProcessLine(Form(".!mv -f %s%04d.root %s%04d.root", pscfilesname0.c_str(),idet, pscfilesname.c_str(),idet));
  }

  cout<<"\e[1;33m Write "<<cPSCfile<<" PSC to files \e[0m"<<endl;
  
  //ClearPSCMem();
  return;
}

void AGATA::LoadPSCfiles(){ // load all Pulse Shape Collection in memory
  LoadPSCfiles(-1);
}

void AGATA::LoadPSCfiles(int detid){ // load all Pulse Shape Collection in memory

  ClearPSCMem();
  //ClosePSCFiles();

  vector <int> idlist; // detid list to load
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  string pscfilesname = "PSCfiles/Det";
  //if(Detid>-1) pscfilesname = "./share/PSCs/Det";
  cout<<"\e[1;33m Load Pulse Shape Collections from "<<pscfilesname;
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;

  
  for(int idet : idlist){
    pscfile[idet] = new TFile(Form("%s%04d.root",pscfilesname.c_str(),idet));

    for(int iseg=0; iseg<NSeg; iseg++){
      psctree[idet][iseg] = (TTree*)pscfile[idet]->Get(Form("tree%d",iseg));

      InitTreeRead(psctree[idet][iseg]);

      int nentries = psctree[idet][iseg]->GetEntriesFast();

      for(int i=0; i<nentries; i++){

	if(cPSCmem%10000==0){
	  double MemUsageGB = GetCurrentMemoryUsage()/GB;
	  double MemTotalGB = GetTotalSystemMemory()/GB;
	  double MemUsage = MemUsageGB / MemTotalGB * 100;
	
	  cout<<"\r Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<"PS-"<<cPStotal
	      <<" PSC-"<<cPSCtotal<<flush;

	  if(MemUsage>MaxMemUsage){
	    cout<<endl<<"exceed memory limit..."<<endl;
	    return;
	  }
	}
	
	psctree[idet][iseg]->GetEntry(i);

	PSC apsc;
	apsc.det = det;
	apsc.seg = seg;
	apsc.index = index;
	apsc.nhits = nhits;

#ifdef REALPOS
	copy_n(labpos, 3, apsc.labpos);
	copy_n(detpos, 3, apsc.detpos);
#endif
	//copy_n(calpos, 3, apsc.calpos);
	//copy_n(cadpos, 3, apsc.cadpos);
#ifdef WITHPS
	if(kWithPS){
	  apsc.cpulsehits = cpulsehits;
	  for(int iiseg=0; iiseg<NSegCore; iiseg++){
	    copy_n(spulse[iiseg], NSig, apsc.spulse[iiseg]);
	  }
	  for(int iiseg=0; iiseg<NSeg_comp; iiseg++){
	    copy_n(cpulse[iiseg], NSig, apsc.cpulse[iiseg]);
	  }
	  copy_n(segwgt, NSeg_comp, apsc.segwgt);
	}
#endif
	
	fPSC[idet][iseg].push_back(apsc);
	cPSCtotal++;   cPSCmem++;
      }
    } //end of loop seg
    
    pscfile[idet]->Close();
  }
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;
  cout<<" Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
      <<"PS-"<<cPStotal
      <<" PSC-"<<cPSCtotal<<endl;

  cout<<"\e[1;33m Load "<<cPSCmem<<" PSC to memory \e[0m"<<endl;

  return;
}


void AGATA::ClearPSCMem(){
  for(int idet = 0; idet<NDets; idet++)
    for(int iseg=0; iseg<NSeg; iseg++){
      cPSCmem-=fPSC[idet][iseg].size();
      fPSC[idet][iseg].clear();
      fPSC[idet][iseg].shrink_to_fit();
    }
  cPSCtotal = 0;
  cPSCmem   = 0;
  cPSCfile  = 0;

  return;
}

void AGATA::GetPSCstat(int *PSCstat){
  PSCstat[0] = cPStotal;
  PSCstat[1] = cPSCtotal;
  PSCstat[2] = cPSCmem;
  PSCstat[3] = cPSCfile;
  PSCstat[4] = maxnhits;
  PSCstat[5] = cPaths;
  PSCstat[6] = cHits;
  PSCstat[7] = cHCs;
  return;
}


//------------------------------------------------
// HitCollection
//------------------------------------------------

// Write files for HitCollection
void AGATA::WriteHCfiles(){
  WriteHCfiles(-1);
}

void AGATA::WriteHCfiles(int detid){

  vector <int> idlist; // detid list to write
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  string hcfilesname = "./share/HCs/Det";
  cout<<"\e[1;33m Create "<<hcfilesname;
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;

  
  int idet, iseg;
  for(int id : idlist){
    idet = id;
    hcfilesname = (string) Form("./share/HCs/Det%04d.root",idet);

    TFile *hcfile = new TFile(hcfilesname.c_str(),"RECREATE");
    TTree *hctree[NSeg];
    for(iseg=0; iseg<NSeg; iseg++){
      cout<<"\r writing hctree for Det "<<idet<<" Seg "<<iseg<<" ..."<<flush;

      hctree[iseg] = new TTree(Form("tree%d",iseg),Form("tree for Hit Collections det%04d seg%d",idet,iseg));
      hctree[iseg]->Branch("det",&det);
      hctree[iseg]->Branch("seg",&seg);
      hctree[iseg]->Branch("index",&index); // PSCid

      hctree[iseg]->Branch("calpos",calpos,"calpos[3]/F"); // init position in lab frame
      hctree[iseg]->Branch("calpos2",calpos2,"calpos2[3]/F"); // calib position in lab frame
#ifdef REALPOS
      hctree[iseg]->Branch("labpos",labpos,"labpos[3]/F"); // real position in lab frame
#endif

      for(HitCollection* ahc : *fHCs[idet][iseg]){
	det = ahc->GetDet();
	seg = ahc->GetSeg();
	index = ahc->GetPid();
	TVector3 tmp;
	tmp = ahc->GetInitPosition();
	calpos[0] = tmp.X(); calpos[1] = tmp.Y(); calpos[2] = tmp.Z();

	tmp = ahc->GetPosition();
	calpos2[0] = tmp.X(); calpos2[1] = tmp.Y(); calpos2[2] = tmp.Z();

#ifdef REALPOS
	tmp = ahc->GetRealPosition();
	labpos[0] = tmp.X(); labpos[1] = tmp.Y(); labpos[2] = tmp.Z();
#endif
	hctree[iseg]->Fill();
      }
      hctree[iseg]->Write();
    }
    hcfile->Close();
  }
  cout<<"\r writing hctree for Det "<<idet<<" Seg "<<iseg<<" ..."<<endl;
  
  return;
}


// Load files for HitCollection
void AGATA::LoadHCfiles(){
  LoadHCfiles(-1);
}

void AGATA::LoadHCfiles(int detid){

  ClearHCMem();

  vector <int> idlist; // detid list to load
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  string hcfilesname = "./share/HCs/Det";
  cout<<"\e[1;33m Load HitCollections from "<<hcfilesname;
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;

  TChain *hctree = new TChain();
  for(int idet : idlist)
    for(int iseg=0; iseg<NSeg; iseg++)
      hctree->AddFile(Form("./share/HCs/Det%04d.root",idet),0,Form("tree%d",iseg));

  hctree->SetBranchAddress("det",&det);
  hctree->SetBranchAddress("seg",&seg);
  hctree->SetBranchAddress("index",&index); // PSCid

  hctree->SetBranchAddress("calpos",calpos); // init position in lab frame
#ifdef REALPOS
  hctree->SetBranchAddress("labpos",labpos); // real position in lab frame
#endif
  
  int Nhcs = hctree->GetEntriesFast();

  for(ihc=0; ihc<Nhcs; ihc++){
    if(ihc%10000==0) cout<<"\r load "<<ihc<<" / "<<Nhcs<<" HitCollections..."<<flush;
    hctree->GetEntry(ihc);

    if(fHCs[det][seg]->size()!=index){ cerr<<"HC index not match!!!"<<endl;}

#ifdef REALPOS
    HitCollection* ahc = new HitCollection(det, seg, index, labpos, calpos);
#else
    float dummy[3] = {0,0,0};
    HitCollection* ahc = new HitCollection(det, seg, index, dummy, calpos);    
#endif
    fHCs[det][seg]->push_back(ahc);
    ahc->SetGid(ihc);
    fAllHCs->push_back(ahc);
    cHCs++;
  }
  cout<<"\r load "<<ihc<<" / "<<Nhcs<<" HitCollections..."<<endl;

  return;
}


void AGATA::ClearHCMem(){
  for(int idet = 0; idet<NDets; idet++)
    for(int iseg=0; iseg<NSeg; iseg++){
      fHCs[idet][iseg]->clear();
      fHCs[idet][iseg]->shrink_to_fit();
    }

  for(HitCollection* ahc : *fAllHCs) delete ahc;
  fAllHCs->clear();
  fAllHCs->shrink_to_fit();
  cHCs = 0;

  return;
}



//------------------------------------------------
// Hit
//------------------------------------------------

// Write files for EvtHits after TreeReaderPulse->GenerateHCs()
void AGATA::WriteEvtHitsfiles(int detid){

  vector <int> idlist; // detid list to write
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  cout<<"\e[1;33m Create hit files for Det";
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;

  
  int iconfig = -1, irun = -1, ientry = -1;
  int tmpconfig, tmprun;

  string hfilename;
  TFile *hfile[MaxNDets];
  TTree *htree[MaxNDets];

  for(int idet : idlist){
    if(idet<0 || idet>=MaxNDets){
      cerr<<"Error: idet = "<<idet<<" !!! "<<endl;
      return;
    }
    
    hfilename = (string) Form("tmp/Det%04d.root",idet);
    hfile[idet] = new TFile(hfilename.c_str(),"RECREATE");
    htree[idet] = new TTree("tree",Form("tree for det%d",idet));
  }
  
  vector<int>    vhcid;
  float          depE;
#ifdef NOISE
  int            noiseidx;
  int            noiseidxshift;
#endif
  Float_t        sourcepos[3];
  Float_t        sourceeng;

  for(ievt=0; ievt<fEventHits->size(); ievt++){ // loop events
    EventHits *fEvent = fEventHits->at(ievt);
    fEvent->GetIdx(tmpconfig, tmprun, ientry);

    if(tmpconfig!=iconfig || tmprun!=irun){ // open new file
      for(int idet : idlist){
	hfilename = (string) Form("./share/Hits/input%d/run%04d_det%04d.root",iconfig,irun,idet);
	cout<<"\r writing htree to "<<hfilename<<" ..."<<flush;
	hfile[idet]->cd();
	htree[idet]->Write();
	hfile[idet]->Close();

	iconfig = tmpconfig;
	irun = tmprun;

	hfilename = (string) Form("./share/Hits/input%d/run%04d_det%04d.root",iconfig,irun,idet);
	hfile[idet] = new TFile(hfilename.c_str(),"RECREATE");
	htree[idet] = new TTree("tree",Form("tree for Hits in config%d run%d det%d",iconfig,irun,idet));
	htree[idet]->Branch("iconfig",&iconfig);
	htree[idet]->Branch("irun",&irun);
	htree[idet]->Branch("ientry",&ientry);
	htree[idet]->Branch("det",&det);
	htree[idet]->Branch("seg",&seg);
	htree[idet]->Branch("hcid",&vhcid);
	htree[idet]->Branch("depE",&depE);
	htree[idet]->Branch("calpos",calpos,"calpos[3]/F");
#ifdef REALPOS
	htree[idet]->Branch("labpos",labpos,"labpos[3]/F");
#endif
#ifdef NOISE
	htree[idet]->Branch("noiseidx",&noiseidx);
	htree[idet]->Branch("noiseidxshift",&noiseidxshift);
#endif
	htree[idet]->Branch("sourcepos",sourcepos,"sourcepos[3]/F");
	htree[idet]->Branch("sourceeng",&sourceeng,"sourceeng/F");
      }
    }

    TVector3 SourcePos = fEvent->GetSourcePos();
    sourcepos[0] = SourcePos.X();
    sourcepos[1] = SourcePos.Y();
    sourcepos[2] = SourcePos.Z();
    sourceeng = fEvent->GetSourceE();

    // get fHit information
    vector<Hit*>* fHits = fEvent->GetfHits();
    for(Hit *ah : *fHits){
      det = ah->GetDet();
      seg = ah->GetSeg();

      bool kfill = false;
      for(int idet : idlist) if(idet==det) kfill = true;
      if(!kfill) continue;
      
      vhcid.clear(); vhcid.shrink_to_fit();
      vector<HitCollection*>* hcs = ah->GetHitCollections();
      for(HitCollection* ahc : *hcs) vhcid.push_back(ahc->GetPid());

      depE = ah->GetE(); // keV

      TVector3 tmp;
      tmp = ah->GetPosition();
      calpos[0] = tmp.X(); calpos[1] = tmp.Y(); calpos[2] = tmp.Z();
#ifdef REALPOS
      tmp = ah->GetRealPosition();
      labpos[0] = tmp.X(); labpos[1] = tmp.Y(); labpos[2] = tmp.Z();
#endif
#ifdef NOISE
      noiseidx = ah->GetNoiseIdx();
      noiseidxshift = ah->GetNoiseIdxShift();
#endif
      htree[det]->Fill();
    }
  }

  for(int idet : idlist){
    hfilename = (string) Form("./share/Hits/input%d/run%04d_det%04d.root",iconfig,irun,idet);
    cout<<"\r writing htree to "<<hfilename<<" ..."<<flush;
    hfile[idet]->cd();
    htree[idet]->Write();
    hfile[idet]->Close();
  }
  cout<<endl;
  
  return;
}


void AGATA::Load(string configfile){
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

      MinRun.push_back(run[0]);
      MaxRun.push_back(run[1]);
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
    cout<<"#input "<<i<<": "<<Form("run %d ~ %d",MinRun[i],MaxRun[i])<<endl;
  }

}


void AGATA::CombEvtHitsfiles(){

  //------------------------------------------------
  // Hit
  //------------------------------------------------
  TFile *hfile[MaxNDets];
  TTree *htree[MaxNDets];
  TFile *hfileall;
  TTree *htreeall;
  
  struct OBJ{
    int          iconfig;
    int          irun;
    int          ientry;
    int          idet;
    int          iseg;
    vector<int> *ihcid = 0;
    float        idepE;
    float        icalpos[3];
#ifdef REALPOS
    float        ilabpos[3];
#endif
#ifdef NOISE
    int          inoiseidx;
    int          inoiseidxshift;
#endif
    float        isourcepos[3];
    float        isourceeng;
  };

  OBJ obj[MaxNDets];
  
  string hfilename0;
  bool kempty;
  int idx[MaxNDets], Nevts[MaxNDets];

  // check exist
  int iconfig, irun;  
  for(iconfig=0; iconfig<nConfig; iconfig++){ // loop configs
    hfilename0 = (string) Form("./share/Hits/input%d",iconfig);
    if(gSystem->AccessPathName(hfilename0.c_str())) return;

    for(irun=MinRun[iconfig]; irun<MaxRun[iconfig]+1; irun++){ // loop runs 

      // check exist
      hfilename0 = (string) Form("./share/Hits/input%d/run%04d",iconfig,irun);
      string hfilename = (string) Form("%s.root",hfilename0.c_str());
      if(gSystem->AccessPathName(hfilename.c_str())){
	hfileall = new TFile(hfilename.c_str(),"RECREATE");
      }else{
	continue;
      }

      
      kempty = true;
      for(int detid=0; detid<NDets; detid++){
	hfilename = (string) Form("%s_det%04d.root",hfilename0.c_str(),detid);
	if(gSystem->AccessPathName(hfilename.c_str())){
	  idx[detid] = -1;
	}else{
	  idx[detid] = 0;
	  hfile[detid] = new TFile(hfilename.c_str());
	  kempty = false;
	}
      }

      if(kempty) continue;

      // load hits from iconfig, irun
      cout<<"\r \e[1;33m Load Hits from "<<hfilename0<<"... \e[0m"<<flush;

      int MaxEntry = -1;
      for(int detid=0; detid<NDets; detid++){

	if(idx[detid]<0){ 
	  cerr<<"iconfig "<<iconfig<<" irun "<<irun<<" detid "<<detid<<" not exist!!!!"<<endl; 
	  continue;
	}

	htree[detid] = (TTree *)hfile[detid]->Get("tree");
	Nevts[detid] = htree[detid]->GetEntriesFast();
	
	htree[detid]->SetBranchAddress("iconfig",       &obj[detid].iconfig);
	htree[detid]->SetBranchAddress("irun",          &obj[detid].irun);
	htree[detid]->SetBranchAddress("ientry",        &obj[detid].ientry);
	htree[detid]->SetBranchAddress("det",           &obj[detid].idet);
	htree[detid]->SetBranchAddress("seg",           &obj[detid].iseg);
	htree[detid]->SetBranchAddress("hcid",          &obj[detid].ihcid);
	htree[detid]->SetBranchAddress("depE",          &obj[detid].idepE);
	htree[detid]->SetBranchAddress("calpos",         obj[detid].icalpos);
#ifdef REALPOS
	htree[detid]->SetBranchAddress("labpos",         obj[detid].ilabpos);
#endif
#ifdef NOISE
	htree[detid]->SetBranchAddress("noiseidx",      &obj[detid].inoiseidx);
	htree[detid]->SetBranchAddress("noiseidxshift", &obj[detid].inoiseidxshift);
#endif
	htree[detid]->SetBranchAddress("sourcepos",      obj[detid].isourcepos);
	htree[detid]->SetBranchAddress("sourceeng",     &obj[detid].isourceeng);

	htree[detid]->GetEntry(Nevts[detid]-1);
	if(obj[detid].ientry>MaxEntry) MaxEntry=obj[detid].ientry;

	htree[detid]->GetEntry(idx[detid]);
      }

      MaxEntry += 1;

      hfileall->cd();
      int oconfig, orun, oentry, odet, oseg;
      vector<int>    ovhcid;
      float          odepE;
      float          ocalpos[3];
      float          olabpos[3];
#ifdef NOISE
      int            onoiseidx;
      int            onoiseidxshift;
#endif
      Float_t        osourcepos[3];
      Float_t        osourceeng;
      htreeall = new TTree("tree",Form("tree for Hits in config%d run%d alldet",iconfig,irun));
      htreeall->Branch("iconfig",&oconfig);
      htreeall->Branch("irun",&orun);
      htreeall->Branch("ientry",&oentry);
      htreeall->Branch("det",&odet);
      htreeall->Branch("seg",&oseg);
      htreeall->Branch("hcid",&ovhcid);
      htreeall->Branch("depE",&odepE);
      htreeall->Branch("calpos",ocalpos,"calpos[3]/F");
#ifdef REALPOS
      htreeall->Branch("labpos",olabpos,"labpos[3]/F");
#endif
#ifdef NOISE
      htreeall->Branch("noiseidx",&onoiseidx);
      htreeall->Branch("noiseidxshift",&onoiseidxshift);
#endif
      htreeall->Branch("sourcepos",osourcepos,"sourcepos[3]/F");
      htreeall->Branch("sourceeng",&osourceeng,"sourceeng/F");

      for(int iety=0; iety<MaxEntry; iety++){ // loop entries

	// find hits belong to the same entry
	for(int detid=0; detid<NDets; detid++){

	  if(idx[detid]<0) continue;

	  while(obj[detid].ientry<=iety){

	    if(iety==obj[detid].ientry){
	      oconfig = obj[detid].iconfig;
	      orun = obj[detid].irun;
	      oentry = obj[detid].ientry;
	      odet = obj[detid].idet;
	      oseg = obj[detid].iseg;
	      for(int id : *obj[detid].ihcid) ovhcid.push_back(id);
	      odepE = obj[detid].idepE;
	      for(int ix=0; ix<3; ix++) ocalpos[ix] = obj[detid].icalpos[ix];
#ifdef REALPOS
	      for(int ix=0; ix<3; ix++) olabpos[ix] = obj[detid].ilabpos[ix];
#endif
#ifdef NOISE
	      onoiseidx = obj[detid].inoiseidx;
	      onoiseidxshift = obj[detid].inoiseidxshift;
#endif
	      for(int ix=0; ix<3; ix++) osourcepos[ix] = obj[detid].isourcepos[ix];
	      osourceeng = obj[detid].isourceeng;
	    
	      htreeall->Fill();
	      ovhcid.clear();
	    }

	    // move to next idx
	    idx[detid]++;
	    if(idx[detid]>=Nevts[detid]) break;
	    htree[detid]->GetEntry(idx[detid]);
	    
	  }
	  
	}

      } // end of loop entries

      for(int detid=0; detid<NDets; detid++){
	if(idx[detid]<0) continue;
	hfile[detid]->Close();
      }

      hfileall->cd();
      htreeall->Write();
      hfileall->Close();
    } // end of loop runs
    cout<<endl;
  } // end of loop configs

  return;
}


void AGATA::LoadEvtHitsfiles(int iconfig){

  //------------------------------------------------
  // Hit
  //------------------------------------------------
  TFile *hfile[MaxNDets];
  TTree *htree[MaxNDets];

  struct OBJ{
    //int          iconfig;
    //int          irun;
    int          ientry;
    int          idet;
    int          iseg;
    vector<int> *ihcid = 0;
    float        idepE;
    float        icalpos[3];
#ifdef REALPOS
    float        ilabpos[3];
#endif
#ifdef NOISE
    int          inoiseidx;
    int          inoiseidxshift;
#endif
    float        isourcepos[3];
    float        isourceeng;
  };

  OBJ obj[MaxNDets];

  string hfilename0;
  bool kempty;
  int idx[MaxNDets], Nevts[MaxNDets];

  // check exist
  hfilename0 = (string) Form("./share/Hits/input%d",iconfig);
  if(gSystem->AccessPathName(hfilename0.c_str())) return;

  int irun;  
  for(; atomrun<MaxRun[iconfig]+1;){ // loop runs 

    {
#ifdef NTHREADS
      lock_guard<mutex> lock(EvtHitsmtx);
#endif
      if(atomrun<MaxRun[iconfig]+1){
	irun = atomrun;
	atomrun++;
      }else{
	continue;
      }
    }
    
    // check exist
    hfilename0 = (string) Form("./share/Hits/input%d/run%04d",iconfig,irun);
    kempty = true;
    for(int detid=0; detid<NDets; detid++){
      string hfilename = (string) Form("%s_det%04d.root",hfilename0.c_str(),detid);
      if(gSystem->AccessPathName(hfilename.c_str())){
	idx[detid] = -1;
      }else{
	idx[detid] = 0;
	hfile[detid] = new TFile(hfilename.c_str());
	kempty = false;
      }
    }

    if(kempty) continue;

    // load hits from iconfig, irun
    cout<<"\r \e[1;33m Load Hits from "<<hfilename0<<"... \e[0m"<<flush;

    int MaxEntry = -1;
    for(int detid=0; detid<NDets; detid++){

      if(idx[detid]<0){ 
	cerr<<"iconfig "<<iconfig<<" irun "<<irun<<" detid "<<detid<<" not exist!!!!"<<endl; 
	continue;
      }

      htree[detid] = (TTree *)hfile[detid]->Get("tree");
      Nevts[detid] = htree[detid]->GetEntriesFast();
	
      //htree[detid]->SetBranchAddress("iconfig",       &obj[detid].iconfig);
      //htree[detid]->SetBranchAddress("irun",          &obj[detid].irun);
      htree[detid]->SetBranchAddress("ientry",        &obj[detid].ientry);
      htree[detid]->SetBranchAddress("det",           &obj[detid].idet);
      htree[detid]->SetBranchAddress("seg",           &obj[detid].iseg);
      htree[detid]->SetBranchAddress("hcid",          &obj[detid].ihcid);
      htree[detid]->SetBranchAddress("depE",          &obj[detid].idepE);
      htree[detid]->SetBranchAddress("calpos",         obj[detid].icalpos);
#ifdef REALPOS
      htree[detid]->SetBranchAddress("labpos",         obj[detid].ilabpos);
#endif
#ifdef NOISE
      htree[detid]->SetBranchAddress("noiseidx",      &obj[detid].inoiseidx);
      htree[detid]->SetBranchAddress("noiseidxshift", &obj[detid].inoiseidxshift);
#endif
      htree[detid]->SetBranchAddress("sourcepos",      obj[detid].isourcepos);
      htree[detid]->SetBranchAddress("sourceeng",     &obj[detid].isourceeng);

      htree[detid]->GetEntry(Nevts[detid]-1);
      if(obj[detid].ientry>MaxEntry) MaxEntry=obj[detid].ientry;

      htree[detid]->GetEntry(idx[detid]);
    }

    MaxEntry += 1;
    if(MaxEntry>0) NevtsTotal += MaxEntry;
      
    vector<int>            vdet;
    vector<int>            vseg;
    vector<vector<int>>    vhcid;
    vector<float>          vdepE;
    vector<vector<double>> vcalpos;
#ifdef REALPOS
    vector<vector<double>> vlabpos;
#endif
#ifdef NOISE
    vector<int>            vnoiseidx;
    vector<int>            vnoiseidxshift;
#endif
    Float_t                sourcepos[3];
    Float_t                sourceeng;

    for(int iety=0; iety<MaxEntry; iety++){ // loop entries
      //if(iety%10==0) cout<<"\r ientry = "<<iety<<flush;
      vdet.clear();
      vseg.clear();
      vhcid.clear();
      vdepE.clear();
      vcalpos.clear();
#ifdef REALPOS
      vlabpos.clear();
#endif
#ifdef NOISE
      vnoiseidx.clear();
      vnoiseidxshift.clear();
#endif

      // find hits belong to the same entry
      for(int detid=0; detid<NDets; detid++){

	if(idx[detid]<0) continue;

	int imove=0;
	while(obj[detid].ientry<iety){
	  imove++;
	  // move to next idx
	  idx[detid]++;
	  if(idx[detid]>=Nevts[detid]) break;
	  htree[detid]->GetEntry(idx[detid]);
	}
	if(imove>1){ cout<<"move "<<imove<<" times!!!"<<endl;}


	if(iety==obj[detid].ientry){

	  vdet.push_back(obj[detid].idet);
	  vseg.push_back(obj[detid].iseg);
	  vhcid.push_back(*obj[detid].ihcid);
	  vdepE.push_back(obj[detid].idepE);

	  vector<double> tmppos1(3);
	  for(int ix=0; ix<3; ix++) tmppos1[ix]=obj[detid].icalpos[ix];
	  vcalpos.push_back(tmppos1);
#ifdef REALPOS
	  vector<double> tmppos2(3);
	  for(int ix=0; ix<3; ix++) tmppos2[ix]=obj[detid].ilabpos[ix];
	  vlabpos.push_back(tmppos2);
#endif
#ifdef NOISE
	  vnoiseidx.push_back(obj[detid].inoiseidx);
	  vnoiseidxshift.push_back(obj[detid].inoiseidxshift);
#endif
	    
	  for(int ix=0; ix<3; ix++) sourcepos[ix]=obj[detid].isourcepos[ix];
	  sourceeng = obj[detid].isourceeng;

	}

      }


      // make EventHits
      float SourceE = sourceeng;
      TVector3 SourcePos(sourcepos[0],sourcepos[1],sourcepos[2]);
      EventHits* fEvent = new EventHits(SourceE, SourcePos);
      fEvent->SetIdx(iconfig,irun,iety);

      for(int i=0; i<vdet.size(); i++){
	int detid = vdet[i];
	int segid = vseg[i];

	float depE = vdepE[i]; // keV
	TVector3 initpos(vcalpos[i][0],vcalpos[i][1],vcalpos[i][2]); // mm
#ifdef REALPOS
	TVector3 hitpos(vlabpos[i][0],vlabpos[i][1],vlabpos[i][2]); // mm
#else
	TVector3 hitpos(0,0,0);
#endif
	Hit *ahit = new Hit(detid, segid, depE, hitpos, initpos);

#ifdef NOISE
	ahit->SetNoiseIdx(vnoiseidx[i]);
	ahit->SetNoiseIdxShift(vnoiseidxshift[i]);
#endif
	fEvent->Add(ahit);
	cHits++;

	
#ifdef NTHREADS
	lock_guard<mutex> lock(PSCmtx[detid][segid]);
#endif
	// connect with HCs
	for(int ii=0; ii<vhcid[i].size(); ii++){
	  int ipsc = vhcid[i][ii];
	  if(ipsc>fHCs[detid][segid]->size()-1){
	    cerr<<"ipsc = "<<ipsc<<" outside fHCs["<<detid<<"]["<<segid<<"] range!!!!"<<endl;
	    continue;
	  }
	    
	  HitCollection *ahc = fHCs[detid][segid]->at(ipsc);
	  ahc->AddHit(ahit);
	  ahit->AddHitCollection(ahc);
	}
      
      }

#ifdef NTHREADS
      lock_guard<mutex> lock(EvtHitsmtx);
#endif
      fEventHits->push_back(fEvent);
      ievt++;

    } // end of loop entries

    for(int detid=0; detid<NDets; detid++){
      if(idx[detid]<0) continue;
      hfile[detid]->Close();
    }
      
  } // end of loop runs

  return;
}


void AGATA::LoadEvtHitsfiles2(int iconfig){

  //------------------------------------------------
  // Hit
  //------------------------------------------------
  TFile *hfile;
  TTree *htree;

  struct OBJ{
    //int          iconfig;
    //int          irun;
    int          ientry;
    int          idet;
    int          iseg;
    vector<int> *ihcid = 0;
    float        idepE;
    float        icalpos[3];
#ifdef REALPOS
    float        ilabpos[3];
#endif
#ifdef NOISE
    int          inoiseidx;
    int          inoiseidxshift;
#endif
    float        isourcepos[3];
    float        isourceeng;
  };

  OBJ obj;

  string hfilename0;
  int idx, Nevts;

  // check exist
  hfilename0 = (string) Form("./share/Hits/input%d",iconfig);
  if(gSystem->AccessPathName(hfilename0.c_str())) return;

  int irun;  
  for(; atomrun<MaxRun[iconfig]+1;){ // loop runs 

    {
#ifdef NTHREADS2
      lock_guard<mutex> lock(EvtHitsmtx);
#endif
      if(atomrun<MaxRun[iconfig]+1){
	irun = atomrun;
	atomrun++;
      }else{
	continue;
      }
    }
    
    // check exist
    hfilename0 = (string) Form("./share/Hits/input%d/run%04d",iconfig,irun);
    string hfilename = (string) Form("%s.root",hfilename0.c_str());
    if(gSystem->AccessPathName(hfilename.c_str())){
      continue;
    }else{
      hfile = new TFile(hfilename.c_str());
    }

    // load hits from iconfig, irun
    cout<<"\r \e[1;33m Load Hits from "<<hfilename0<<"... \e[0m"<<flush;

    int MaxEntry = -1;

    htree = (TTree *)hfile->Get("tree");
    if(!htree){
      hfile->Close();
      continue;
    }
    Nevts = htree->GetEntriesFast();
	
    //htree->SetBranchAddress("iconfig",       &obj.iconfig);
    //htree->SetBranchAddress("irun",          &obj.irun);
    htree->SetBranchAddress("ientry",        &obj.ientry);
    htree->SetBranchAddress("det",           &obj.idet);
    htree->SetBranchAddress("seg",           &obj.iseg);
    htree->SetBranchAddress("hcid",          &obj.ihcid);
    htree->SetBranchAddress("depE",          &obj.idepE);
    htree->SetBranchAddress("calpos",         obj.icalpos);
#ifdef REALPOS
    htree->SetBranchAddress("labpos",         obj.ilabpos);
#endif
#ifdef NOISE
    htree->SetBranchAddress("noiseidx",      &obj.inoiseidx);
    htree->SetBranchAddress("noiseidxshift", &obj.inoiseidxshift);
#endif
    htree->SetBranchAddress("sourcepos",      obj.isourcepos);
    htree->SetBranchAddress("sourceeng",     &obj.isourceeng);

    htree->GetEntry(Nevts-1);
    MaxEntry=obj.ientry;

    MaxEntry += 1;
    if(MaxEntry>0) NevtsTotal += MaxEntry;

    vector<int>            vdet;
    vector<int>            vseg;
    vector<vector<int>>    vhcid;
    vector<float>          vdepE;
    vector<vector<double>> vcalpos;
#ifdef REALPOS
    vector<vector<double>> vlabpos;
#endif
#ifdef NOISE
    vector<int>            vnoiseidx;
    vector<int>            vnoiseidxshift;
#endif
    Float_t                sourcepos[3];
    Float_t                sourceeng;

    idx = 0;
    for(int iety=0; iety<MaxEntry; iety++){ // loop entries
      //if(iety%10==0) cout<<"\r ientry = "<<iety<<flush;
      vdet.clear();
      vseg.clear();
      vhcid.clear();
      vdepE.clear();
      vcalpos.clear();
#ifdef REALPOS
      vlabpos.clear();
#endif
#ifdef NOISE
      vnoiseidx.clear();
      vnoiseidxshift.clear();
#endif

      // find hits belong to the same entry
      for(; idx<Nevts; idx++){

	htree->GetEntry(idx);
	if(obj.ientry>iety) break;
	
	if(iety==obj.ientry){

	  vdet.push_back(obj.idet);
	  vseg.push_back(obj.iseg);
	  vhcid.push_back(*obj.ihcid);
	  vdepE.push_back(obj.idepE);

	  vector<double> tmppos1(3);
	  for(int ix=0; ix<3; ix++) tmppos1[ix]=obj.icalpos[ix];
	  vcalpos.push_back(tmppos1);
#ifdef REALPOS
	  vector<double> tmppos2(3);
	  for(int ix=0; ix<3; ix++) tmppos2[ix]=obj.ilabpos[ix];
	  vlabpos.push_back(tmppos2);
#endif
#ifdef NOISE
	  vnoiseidx.push_back(obj.inoiseidx);
	  vnoiseidxshift.push_back(obj.inoiseidxshift);
#endif
	    
	  for(int ix=0; ix<3; ix++) sourcepos[ix]=obj.isourcepos[ix];
	  sourceeng = obj.isourceeng;

	}

      }


      // make EventHits
      float SourceE = sourceeng;
      TVector3 SourcePos(sourcepos[0],sourcepos[1],sourcepos[2]);
      EventHits* fEvent = new EventHits(SourceE, SourcePos);
      fEvent->SetIdx(iconfig,irun,iety);

      for(int i=0; i<vdet.size(); i++){
	int detid = vdet[i];
	int segid = vseg[i];

	float depE = vdepE[i]; // keV
	TVector3 initpos(vcalpos[i][0],vcalpos[i][1],vcalpos[i][2]); // mm
#ifdef REALPOS
	TVector3 hitpos(vlabpos[i][0],vlabpos[i][1],vlabpos[i][2]); // mm
#else
	TVector3 hitpos(0,0,0);
#endif
	Hit *ahit = new Hit(detid, segid, depE, hitpos, initpos);

#ifdef NOISE
	ahit->SetNoiseIdx(vnoiseidx[i]);
	ahit->SetNoiseIdxShift(vnoiseidxshift[i]);
#endif
	fEvent->Add(ahit);
	cHits++;

	
#ifdef NTHREADS2
	lock_guard<mutex> lock(PSCmtx[detid][segid]);
#endif
	// connect with HCs
	for(int ii=0; ii<vhcid[i].size(); ii++){
	  int ipsc = vhcid[i][ii];
	  if(ipsc>fHCs[detid][segid]->size()-1){
	    cerr<<"ipsc = "<<ipsc<<" outside fHCs["<<detid<<"]["<<segid<<"] range!!!!"<<endl;
	    continue;
	  }
	    
	  HitCollection *ahc = fHCs[detid][segid]->at(ipsc);
	  ahc->AddHit(ahit);
	  ahit->AddHitCollection(ahc);
	}
      
      }

#ifdef NTHREADS2
      lock_guard<mutex> lock(EvtHitsmtx);
#endif
      fEventHits->push_back(fEvent);
      ievt++;

    } // end of loop entries

    hfile->Close();
      
  } // end of loop runs

  return;
}


void AGATA::LoadEvtHitsconfigs(){

  if(Detid>-1) return; // only work for load all detectors
  if(cHCs==0){
    cerr<<"Need to load HCs first!!!"<<endl;
    return;
  }
  
  ClearEvtHitsMem();

#ifdef NTHREADS2
  cout<<"\e[1;33m using "<<NTHREADS2<<" threads:";
#endif
  cout<<"\e[1;33m Load EventHits... \e[0m"<<endl;
  
  int iconfig;
  NevtsTotal = 0;
  ievt = 0;

  for(iconfig=0; iconfig<nConfig; iconfig++){ // loop input config
    // check exist
    string hfilename0 = (string) Form("./share/Hits/input%d",iconfig);
    if(gSystem->AccessPathName(hfilename0.c_str())) continue;

    atomrun = MinRun[iconfig];
    
#ifndef NTHREADS2
    LoadEvtHitsfiles2(iconfig);
#else
    thread th[NTHREADS2];
  
    for(int i=0; i<NTHREADS2; i++){
      th[i] = thread(&AGATA::LoadEvtHitsfiles2, this, iconfig);
    }

    for(int i=0; i<NTHREADS2; i++){
      if(th[i].joinable())
	th[i].join();
    }
#endif
    
    hfilename0 = (string) Form("./share/Hits/input%d/run%04d",iconfig,atomrun-1);
    cout<<"\r \e[1;33m Load Hits from "<<hfilename0<<"... \e[0m"<<endl;
  } // end of loop input iconfig
  
  cout<<" load "<<ievt<<" / "<<NevtsTotal<<" EventHits..."<<endl;
  
  return;
}


void AGATA::ClearEvtHitsMem(){

  for(EventHits* fEvent : *fEventHits) delete fEvent;
  fEventHits->clear();
  fEventHits->shrink_to_fit();
  cHits = 0;
  
  for(Path* apath : *fPaths) delete apath;
  fPaths->clear();
  fPaths->shrink_to_fit();
  cPaths = 0;

  return;
}



void AGATA::SetAddNewPSC(bool val){
  for(int detid=0; detid<MaxNDets; detid++)
    for(int segid=0; segid<NSeg; segid++)
      kAddNewPSC[detid][segid] = val;
}

bool AGATA::IsNewPSC(){ // if new PSC added
  for(int detid=0; detid<MaxNDets; detid++)
    for(int segid=0; segid<NSeg; segid++)
      if(kAddNewPSC[detid][segid]) return true;
  return false;
}

void AGATA::SortPSC(){
  for(int idet=0; idet<NDets; idet++){
    for(int iseg=0; iseg<NSeg; iseg++){
      // increasing ientry
      sort(fPSC[idet][iseg].begin(),fPSC[idet][iseg].end(),
	   [](const PSC& lhs, const PSC& rhs){
	     return lhs.index<rhs.index;});
    }
  }  
  return;
}

//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// group pulse shape
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
Double_t AGATA::AddPStoPSC(PS *aps, Hit *ahit, vector<int> &entrylist){
  int detid = aps->det;
  int segid = aps->seg;

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

#ifdef REALPOS
  TMatrixD pos(3,1); for(int ix=0; ix<3; ix++) pos(ix,0) = aps->labpos[ix];
#endif

  float maxchi2;
  float maxchi2s[3];
  float dymaxchi2 = MaxChi2;
  float dymaxchi2s[3] = {MaxChi2s[0], MaxChi2s[1], MaxChi2s[2]};
  Double_t minchi2 = 1e9;

  int level = ahit->GetLevel();
  int npsc;
  if(level==1) npsc = fPSC[detid][segid].size(); // initial group fHCs
  else if(level==2) npsc = divHCs[detid][segid].size(); // divide fHCs
  else return minchi2; // skip hit

  int ipsc;
  for(int i=0; i<npsc; i++){
    if(level==1) ipsc=i;
    else if(level==2) ipsc=divHCs[detid][segid][i];
    else continue;

    if(ipsc<0 || ipsc>fHCs[detid][segid]->size()-1) continue;

    if(ahit->InHitCollection(fHCs[detid][segid]->at(ipsc))) continue; // skip if already in HC
    
    // dynamic chi2 range from HC
    if(kGroupPos) maxchi2 = fHCs[detid][segid]->at(ipsc)->MaxChi2s[0];
    else{
      for(int ix=0; ix<3; ix++) maxchi2s[ix] = fHCs[detid][segid]->at(ipsc)->MaxChi2s[ix];
    }
      
    float chi2;
    float chi2s[3];
    if(kGroupPos) chi2 = Dist(aps,&fPSC[detid][segid][ipsc]); // group according to real position
    else          chi2 = Chi2Fast3D(aps,&fPSC[detid][segid][ipsc],MaxChi2s,chi2s,true);
    //else          chi2 = Chi2Fast(aps,&fPSC[detid][segid][ipsc],MaxChi2,true);

    if(chi2<minchi2) minchi2 = chi2;

    if((kGroupPos && chi2<MaxChi2) ||
       (!kGroupPos && chi2s[0]<MaxChi2s[0] && chi2s[1]<MaxChi2s[1] && chi2s[2]<MaxChi2s[2])){ // within initial chi2 range

      // get min dynamic chi2 range for new PSC
      if(kGroupPos){
	if(maxchi2<dymaxchi2) dymaxchi2 = maxchi2;
      }else{
	for(int ix=0; ix<3; ix++)
	  if(maxchi2s[ix]<dymaxchi2s[ix]) dymaxchi2s[ix] = maxchi2s[ix];
      }
      
      // add to PSC
      if((kGroupPos && chi2<maxchi2) ||
	 (!kGroupPos && chi2s[0]<maxchi2s[0] && chi2s[1]<maxchi2s[1] && chi2s[2]<maxchi2s[2])){ // within dynamic chi2 range

	entrylist.push_back(fPSC[detid][segid][ipsc].index);

	double nhitstmp = (double)fPSC[detid][segid][ipsc].nhits;
#ifdef REALPOS
	// calc average position
	TMatrixD LabPos(3,1);
	for(int it=0; it<3; it++) LabPos(it,0) = fPSC[detid][segid][ipsc].labpos[it];
	LabPos = 1./(nhitstmp+1)*(nhitstmp*LabPos + pos);
	TMatrixD DetPos = agatageo->Lab2DetPos(detid,LabPos);
      
	for(int it=0; it<3; it++){
	  fPSC[detid][segid][ipsc].labpos[it] = LabPos(it,0);
	  fPSC[detid][segid][ipsc].detpos[it] = DetPos(it,0);
	}
#endif
      
#ifdef WITHPS
	if(kWithPS){
	  // calc average pulse
	  for(int iseg=0; iseg<NSegCore; iseg++)
	    for(int isig=0; isig<NSig; isig++)
	      fPSC[detid][segid][ipsc].spulse[iseg][isig] =
		1./(nhitstmp+1)*(fPSC[detid][segid][ipsc].spulse[iseg][isig]*nhitstmp +
				 aps->opulse[iseg][isig]);
	}
#endif
      
	fPSC[detid][segid][ipsc].nhits++;
	if(fPSC[detid][segid][ipsc].nhits>maxnhits) maxnhits = fPSC[detid][segid][ipsc].nhits;

	// add to HitCollection
	fHCs[detid][segid]->at(ipsc)->AddHit(ahit);
#ifdef REALPOS
	fHCs[detid][segid]->at(ipsc)->SetRealPosition(fPSC[detid][segid][ipsc].labpos);
#endif
	ahit->AddHitCollection(fHCs[detid][segid]->at(ipsc));

      }
    }
  }//end of loop fPSC

  
  //***************************************************//
  // new PSC
  //***************************************************//
  if(ahit->hasHitCollection()==0 &&
     aps->energy>PSCEMIN &&
     (fPSC[detid][segid].size()<PSClimit[detid]
      || freeHCs[detid][segid].size()>0
      || PSClimit[detid]<0)){

    float dummycalpos[3];
    float dummylabpos[3];

    TVector3 tmppos = ahit->GetPosition(); // initial pos from ahit
    TMatrixD CalPos(3,1);
    CalPos(0,0)=tmppos.X();  CalPos(1,0)=tmppos.Y();  CalPos(2,0)=tmppos.Z();
    //TMatrixD CadPos = agatageo->Lab2DetPos(detid,CalPos);   
#ifdef REALPOS
    TMatrixD DetPos = agatageo->Lab2DetPos(detid,pos);
#endif

    for(int ix=0; ix<3; ix++){
      dummycalpos[ix] = CalPos(ix,0);
#ifdef REALPOS
      dummylabpos[ix] = pos(ix,0);
#else
      dummylabpos[ix] = 0;
#endif
    }

    int idx;
    if(freeHCs[detid][segid].size()>0){ // use old free PSC
      idx = freeHCs[detid][segid][0];
      freeHCs[detid][segid].erase(freeHCs[detid][segid].begin());

      fHCs[detid][segid]->at(idx)->SetRealPosition(dummylabpos);
      fHCs[detid][segid]->at(idx)->SetPosition(dummycalpos);
      fHCs[detid][segid]->at(idx)->SetInitPosition(dummycalpos);

    }else{ // add new PSC
      PSC newpsc;
      newpsc.det = detid;
      newpsc.seg = segid;
      newpsc.index = fPSC[detid][segid].size();

      fPSC[detid][segid].push_back(newpsc);
      idx = newpsc.index;

      // new HitCollection
      HitCollection* newhc = new HitCollection(detid, segid, idx, dummylabpos, dummycalpos);
      fHCs[detid][segid]->push_back(newhc);

#ifdef NTHREADS
      lock_guard<mutex> lock2(AllHCmtx); // lock fAllHCs
#endif
      newhc->SetGid(fAllHCs->size());
      fAllHCs->push_back(newhc);
    }
    cPSCtotal++;   cPSCmem++;    cHCs++;
    
    PSC *apsc = &fPSC[detid][segid][idx];
    apsc->nhits = 1;
    HitCollection* ahc = fHCs[detid][segid]->at(idx);
    if(level==2){
      ahc->Marker = 3;
      divHCs[detid][segid].push_back(idx);
    }
    
    for(int ix=0; ix<3; ix++){
      //apsc->calpos[ix] = CalPos(ix,0);
      //apsc->cadpos[ix] = CadPos(ix,0);
#ifdef REALPOS
      apsc->labpos[ix] = pos(ix,0);
      apsc->detpos[ix] = DetPos(ix,0);
      apsc->cpos[ix] = apsc->detpos[ix];
#endif
    }

#ifdef WITHPS
    if(kWithPS){
      apsc->cpulsehits = aps->nhits;
      for(int iseg=0; iseg<NSeg_comp; iseg++){
	copy_n(aps->apulse[iseg], NSig, apsc->cpulse[iseg]);
      }
      copy_n(aps->segwgt, NSeg_comp, apsc->segwgt);
      
      for(int iseg=0; iseg<NSegCore; iseg++){
	copy_n(aps->opulse[iseg], NSig, apsc->spulse[iseg]);
      }
    }
#endif
    
    if(kGroupPos){
      ahc->MaxChi2s[0] = dymaxchi2;
    }else{
      for(int ix=0; ix<3; ix++) ahc->MaxChi2s[ix] = dymaxchi2s[ix];
    }
    ahc->AddHit(ahit);
    ahit->AddHitCollection(ahc);

    kAddNewPSC[detid][segid] = true;
    
    entrylist.push_back(idx);
  }

  // count if previously has no HC and now has HCs
  if(ahit->hasHitCollection()>0 && ahit->hasHitCollection()==entrylist.size()){
    cPStotal++;
    cHits++;
  }
  
  return minchi2;
}


void AGATA::RemovePSfromPSC(PS *aps, Hit *ahit){ // remove ahit from all PSCs
  int detid = aps->det;
  int segid = aps->seg;

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  vector<HitCollection*>* hcs = ahit->GetHitCollections();

#ifdef REALPOS
  TMatrixD pos(3,1); for(int ix=0; ix<3; ix++) pos(ix,0) = aps->labpos[ix];
#endif
  
  for(HitCollection* ahc : *hcs){

    int ipsc = ahc->GetPid();

    double nhitstmp = (double)fPSC[detid][segid][ipsc].nhits;
#ifdef REALPOS
    // calc average position
    TMatrixD LabPos(3,1);
    for(int it=0; it<3; it++) LabPos(it,0) = fPSC[detid][segid][ipsc].labpos[it];
#endif

    if(nhitstmp<1){
      cerr<<"something wrong, det = "<<detid<<" seg = "<<segid<<" ipsc = "<<ipsc
	  <<" fPSC.nhits = "<<fPSC[detid][segid][ipsc].nhits
	  <<" fHC.nhits = "<<fHCs[detid][segid]->at(ipsc)->GetSize()<<endl;
      continue;
    }
    double ftmp = 1;
    if(nhitstmp>1) ftmp = 1./(nhitstmp-1);

#ifdef REALPOS
    LabPos = ftmp*(nhitstmp*LabPos - pos);
    TMatrixD DetPos = agatageo->Lab2DetPos(detid,LabPos);
      
    for(int it=0; it<3; it++){
      fPSC[detid][segid][ipsc].labpos[it] = LabPos(it,0);
      fPSC[detid][segid][ipsc].detpos[it] = DetPos(it,0);
    }
#endif
    
#ifdef WITHPS
    if(kWithPS){
      // calc average pulse
      for(int iseg=0; iseg<NSegCore; iseg++)
	for(int isig=0; isig<NSig; isig++)
	  fPSC[detid][segid][ipsc].spulse[iseg][isig] =
	    ftmp*(fPSC[detid][segid][ipsc].spulse[iseg][isig]*nhitstmp -
		  aps->opulse[iseg][isig]);
    }
#endif
      
    fPSC[detid][segid][ipsc].nhits--;
    if(fPSC[detid][segid][ipsc].nhits==0){ cPSCtotal--; cPSCmem--;}

    // remove from HitCollection
    fHCs[detid][segid]->at(ipsc)->RemoveHit(ahit);
#ifdef REALPOS
    fHCs[detid][segid]->at(ipsc)->SetRealPosition(fPSC[detid][segid][ipsc].labpos);
#endif
    if(fHCs[detid][segid]->at(ipsc)->GetSize()==0){
      freeHCs[detid][segid].push_back(ipsc);
      cHCs--;
    }

  }//end of loop hcs

  ahit->ClearHitCollection();

  ahit->SetLevel(0);
  cPStotal--;
  cHits--;

  return;
}


void AGATA::RemovePSC(HitCollection *ahc){ // empty ahc from fHCs
  int detid = ahc->GetDet();
  int segid = ahc->GetSeg();
  int idx = ahc->GetPid();

  vector<Hit*>* fhits = ahc->GetHits();
  if(fhits->size()==0) return;
  
  for(Hit* ah : *fhits){
    ah->RemoveHitCollection(ahc);
    if(ah->hasHitCollection()==0){
      ah->SetLevel(0);
      cPStotal--;
      cHits--;
    }
  }
  ahc->Clear();
  ahc->Marker = 0;
  cPSCtotal--; cPSCmem--;
   
  fPSC[detid][segid][idx].nhits = 0;
#ifdef REALPOS
  for(int ix=0; ix<3; ix++){
    fPSC[detid][segid][idx].labpos[ix] = 0;
    fPSC[detid][segid][idx].detpos[ix] = 0;
    fPSC[detid][segid][idx].cpos[ix] = 0;
  }
#endif
#ifdef WITHPS
  fPSC[detid][segid][idx].cpulsehits = 0;
  for(int iseg=0; iseg<NSeg_comp; iseg++){
    fPSC[detid][segid][idx].segwgt[iseg] = 0;
    for(int isig=0; isig<NSig; isig++)
      fPSC[detid][segid][idx].cpulse[iseg][isig] = 0;
  }
  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++)
      fPSC[detid][segid][idx].spulse[iseg][isig] = 0;
  }
#endif

  freeHCs[detid][segid].push_back(idx);
  cHCs--;
  
  return;
}


void AGATA::DividePSC(HitCollection *ahc, double factor){ // divide ahc by reduced chi2s 
  for(int ix=0; ix<3; ix++) ahc->MaxChi2s[ix] = ahc->MaxChi2s[ix]*factor;
  ahc->MaxChi2s[1] = ahc->MaxChi2s[1]*factor;
  
  int detid = ahc->GetDet();
  int segid = ahc->GetSeg();
  int idx = ahc->GetPid();

  vector<Hit*>* fhits = ahc->GetHits();
  if(fhits->size()==0) return;
  
  for(Hit* ah : *fhits){
    ah->SetLevel(2); // check for divide
    ah->RemoveHitCollection(ahc);
    if(ah->hasHitCollection()==0){
      cPStotal--;
      cHits--;
    }
  }
  ahc->Clear();
   
  fPSC[detid][segid][idx].nhits = 0;
#ifdef REALPOS
  for(int ix=0; ix<3; ix++){
    fPSC[detid][segid][idx].labpos[ix] = 0;
    fPSC[detid][segid][idx].detpos[ix] = 0;
  }
#endif
#ifdef WITHPS
  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++)
      fPSC[detid][segid][idx].spulse[iseg][isig] = 0;
  }
#endif

  ahc->Marker = 2;
  divHCs[detid][segid].push_back(idx);
  
  return;
}


void AGATA::CheckPSCs(int minhits, int maxhits){

  // remove HCs with nhits<minhits
  if(minhits>0){
    int counter = 0;
    for(HitCollection* ahc : *fAllHCs){
      if(ahc->GetSize()<minhits){
	RemovePSC(ahc);
	counter++;
	if(counter%10000==0)
	  cout<<"\r Remove "<<counter<<" HCs with nhits<"<<minhits<<flush;
      }
    }
    cout<<"\r Remove "<<counter<<" HCs with nhits<"<<minhits<<endl;
  }

  // divide HCs with nhits>2*maxhits
  if(maxhits>0){
    maxnhits=0;
    for(int detid=0; detid<MaxNDets; detid++)
      for(int segid=0; segid<NSeg; segid++)
	divHCs[detid][segid].clear();
    
    int counter = 0;
    for(HitCollection* ahc : *fAllHCs){
      if(ahc->GetSize()>2*maxhits){
	double factor = 1.*maxhits/ahc->GetSize();
	factor = pow(factor, 1./4);
	DividePSC(ahc, factor);
	counter++;
	if(counter%10000==0)
	  cout<<"\r Divide "<<counter<<" HCs with nhits>"<<2*maxhits<<flush;
      }else if(ahc->GetSize()>maxnhits){
	maxnhits = ahc->GetSize();
      }
    }
    cout<<"\r Divide "<<counter<<" HCs with nhits>"<<2*maxhits<<endl;
  }
  
  return;
}


void AGATA::ClearHitLevelMarker(int val){
  int counter = 0;
  for(EventHits* fEvent : *fEventHits){
    vector<Hit*>* fhits = fEvent->GetfHits();
    for(Hit* ah : *fhits){
      if(ah->GetLevel()==val){
	ah->SetLevel(1);
	counter++;
      }
    }
  }
  cout<<endl<<" Clear Level Marker "<<val<<" for "<<counter<<" hits"<<endl;
  
  return;
}


Float_t AGATA::Dist(PS *aps, PSC *apsc){
  float dist = 0;
#ifdef REALPOS
  for(int ix=0; ix<3; ix++) dist += SQ(aps->detpos[ix]-apsc->cpos[ix]);
  dist = sqrt(dist);
#endif
  return dist;
}

Float_t AGATA::Chi2Fast3D(PS *aps, PSC *apsc, float *Chi2Limits, float *chi2s, bool kFast){
  float asegpulse[NSig_comp], bsegpulse[NSig_comp];
#ifdef WITHPS
  // compare fired seg, core and neighbor segment
  for(int i=0; i<3; i++) chi2s[i] = 0;

  for(int i=0; i<3; i++){
    for(int j=0; j<2; j++){
      int iseg = 2*i+j;
      if(aps->segwgt[iseg]>0 || apsc->segwgt[iseg]>0){
	copy_n( aps->apulse[iseg], NSig_comp, asegpulse);
	copy_n(apsc->cpulse[iseg], NSig_comp, bsegpulse);
	float tmpchi2 = Chi2seg(asegpulse, bsegpulse);

	chi2s[i] += tmpchi2; // sum
	if(kFast) if(chi2s[i]>Chi2Limits[i]) return chi2s[0]; // interrupt if exceed Chi2Limits
      }
    }
  }
#endif

  return chi2s[0];
}

Float_t AGATA::Chi2Fast(PS *aps, PSC *apsc, float Chi2Limit, bool kFast){
  int nfired = 0;
  float chi2 = 0;
  float asegpulse[NSig_comp], bsegpulse[NSig_comp];
#ifdef WITHPS
  // first compare fired seg, core and neighbor segment, then others
  for(int iseg=0; iseg<NSeg_comp; iseg++){
    if(aps->segwgt[iseg]>0 || apsc->segwgt[iseg]>0){
      copy_n( aps->apulse[iseg], NSig_comp, asegpulse);
      copy_n(apsc->cpulse[iseg], NSig_comp, bsegpulse);
      float tmpchi2 = Chi2seg(asegpulse, bsegpulse);
      //tmpchi2 = aps->segwgt[iseg]>0? tmpchi2*aps->segwgt[iseg] : tmpchi2*apsc->segwgt[iseg];
      //if(tmpchi2>chi2) chi2 = tmpchi2; // maximum
      chi2 += tmpchi2; // sum
      if(kFast) if(chi2>Chi2Limit) return chi2; // interrupt if exceed Chi2Limit
      nfired++;
    }
  }
#endif
  //if(nfired>0) chi2 = chi2/nfired;

  return chi2;
}


Float_t AGATA::Chi2seg(const float *apulse, const float *bpulse){

#ifdef SSE_M256

  const __m256 masks = _mm256_set1_ps(-0.0f);
  const __m256 misig = _mm256_set1_ps(0.01); // min sigma
  
  __m256* realtrace = (__m256*)apulse;
  __m256* basetrace = (__m256*)bpulse;

  __m256 diff = _mm256_setzero_ps();
  __m256 sigm = _mm256_setzero_ps();
  __m256 chis = _mm256_setzero_ps();

  for(int nn=0; nn<LOOP_SSE8_seg; nn++){
    diff = _mm256_sub_ps(realtrace[nn], basetrace[nn]);
#if   CHI2 == CHI2_SQ
    chis = _mm256_add_ps(chis, _mm256_mul_ps(diff,diff)); // square
#elif CHI2 == CHI2_CHI2
    sigm = _mm256_max_ps(_mm256_andnot_ps(masks, basetrace[nn]), misig);       // chi2
    chis = _mm256_add_ps(chis, _mm256_div_ps(_mm256_mul_ps(diff,diff), sigm)); // chi2
#elif CHI2 == CHI2_CHI2_2
    sigm = _mm256_max_ps(_mm256_andnot_ps(masks, basetrace[nn]), sigm); // chi2_2
    chis = _mm256_add_ps(chis, _mm256_mul_ps(diff,diff));               // chi2_2
#elif CHI2 == CHI2_ABS
    chis = _mm256_add_ps(chis, _mm256_andnot_ps(masks, diff)); // ABS value
#elif CHI2 == CHI2_FABS
    sigm = _mm256_max_ps(_mm256_andnot_ps(masks, basetrace[nn]), misig);            // fABS
    chis = _mm256_add_ps(chis, _mm256_div_ps(_mm256_andnot_ps(masks, diff), sigm)); // fABS
#elif CHI2 == CHI2_SQRT
    chis = _mm256_add_ps(chis, _mm256_sqrt_ps(_mm256_andnot_ps(masks, diff))); // sqrt
#elif CHI2 == CHI2_2SQRT
    chis = _mm256_add_ps(chis, _mm256_sqrt_ps(_mm256_sqrt_ps(_mm256_andnot_ps(masks, diff)))); // 2sqrt
#else
    chis = _mm256_add_ps(chis, _mm256_mul_ps(diff,diff)); // square
#endif
  }

#if   CHI2 == CHI2_CHI2_2
  sigm = _mm256_max_ps(sigm, misig); // chi2_2
  chis = _mm256_div_ps(chis, sigm);  // chi2_2
#endif
  
  __m256 tmp = _mm256_permute2f128_ps(chis, chis, 1);
  chis = _mm256_add_ps(chis, tmp);
  chis = _mm256_hadd_ps(chis, chis);
  chis = _mm256_hadd_ps(chis, chis);

  float chi2 = _mm256_cvtss_f32(chis);

#else // use m128

  const __m128 masks = _mm_set1_ps(-0.0f);
  const __m128 zeros = _mm_setzero_ps();
  const __m128 misig = _mm_set1_ps(0.01); // min sigma

  __m128* realtrace = (__m128*)apulse;
  __m128* basetrace = (__m128*)bpulse;

  __m128 diff = _mm_setzero_ps();
  __m128 sigm = _mm_setzero_ps();
  __m128 chis = _mm_setzero_ps();

  for(int nn=0; nn<LOOP_SSE4_seg; nn++){
    diff = _mm_sub_ps(realtrace[nn], basetrace[nn]);
#if   CHI2 == CHI2_SQ
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff)); // square
#elif CHI2 == CHI2_CHI2
    //sigm = _mm_andnot_ps(masks, basetrace[nn]);
    sigm = _mm_max_ps(_mm_andnot_ps(masks, basetrace[nn]), misig);    // chi2
    chis = _mm_add_ps(chis, _mm_div_ps(_mm_mul_ps(diff,diff), sigm)); // chi2
#elif CHI2 == CHI2_CHI2_2
    sigm = _mm_max_ps(_mm_andnot_ps(masks, basetrace[nn]), sigm); // chi2_2
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff));               // chi2_2
#elif CHI2 == CHI2_ABS
    chis = _mm_add_ps(chis, _mm_andnot_ps(masks, diff)); // ABS value
#elif CHI2 == CHI2_FABS
    sigm = _mm_max_ps(_mm_andnot_ps(masks, basetrace[nn]), misig);         // fABS
    chis = _mm_add_ps(chis, _mm_div_ps(_mm_andnot_ps(masks, diff), sigm)); // fABS
#elif CHI2 == CHI2_SQRT
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_andnot_ps(masks, diff))); // sqrt
#elif CHI2 == CHI2_2SQRT
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_sqrt_ps(_mm_andnot_ps(masks, diff)))); // 2sqrt
#else
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff)); // square
#endif
  }

#if   CHI2 == CHI2_CHI2_2
  sigm = _mm_max_ps(sigm, misig); // chi2_2
  chis = _mm_div_ps(chis, sigm);  // chi2_2
#endif
  
  chis = _mm_hadd_ps(_mm_hadd_ps(chis, zeros), zeros);

  float chi2 = _mm_cvtss_f32(chis);
#endif // end of ifdef SSE_M256
  
  return chi2;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// EventHits
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
void AGATA::SortEventHits(){
  // sort EventHits increasing iconfig, irun, ientry
  sort( fEventHits->begin(), fEventHits->end(),
	[](EventHits* lEvt, EventHits* rEvt){
	  int lconfig, lrun, lentry;
	  int rconfig, rrun, rentry;
	  lEvt->GetIdx(lconfig,lrun,lentry);
	  rEvt->GetIdx(rconfig,rrun,rentry);

	  return ((lconfig<rconfig) ||
		  (lconfig==rconfig && lrun<rrun) ||
		  (lconfig==rconfig && lrun==rrun && lentry<rentry));});
}


int AGATA::AddEventHits(EventHits* fEvent){
#ifdef NTHREADS
  lock_guard<mutex> lock(EvtHitsmtx);
#endif
  fEventHits->push_back(fEvent);
  return fEventHits->size()-1;
}


int AGATA::FindiEvtHit(int iconfig, int irun, int ientry, int &istart){
  int iEvtHit = -1;
  int tmpconf, tmprun, tmpetry;
  for(int i=istart; i<fEventHits->size(); i++){
    fEventHits->at(i)->GetIdx(tmpconf, tmprun, tmpetry);

    if(tmpconf<iconfig) continue; else if(tmpconf>iconfig) break;
    if(tmprun<irun)     continue; else if(tmprun>irun)     break;
    if(tmpetry<ientry)  continue; else if(tmpetry>ientry)  break;

    istart = iEvtHit = i;
    break;
  }
  return iEvtHit;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// tracking
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
void AGATA::TrackingLoop(){
  int Nevts = fEventHits->size();
  time(&start);

  for(; ievt<Nevts;){ // loop events

    // get Hits for ievt
    vector<Hit*>* fHits;
    TVector3 sourcepos; // mm
    Double_t EGamma; // keV
    {
#ifdef NTHREADS2
      lock_guard<mutex> lock(EvtHitsmtx);
#endif

      if(ievt%10000==0){
	time(&stop);
	double MemUsageGB = GetCurrentMemoryUsage()/GB;
	double MemTotalGB = GetTotalSystemMemory()/GB;
	double MemUsage = MemUsageGB / MemTotalGB * 100;
	
	cout<<"\r finish read "<<ievt<<" / "<<Nevts<<" evts"
	    <<Form("(%.0fs/kevts)..",difftime(stop,start))
	    <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	    <<"Hits-"<<cHits
	    <<" HCs-"<<cHCs<<".."
	    <<"Paths-"<<cPaths<<flush;

	if(MemUsage>MaxMemUsage){
	  cout<<endl<<"exceed memory limit. Write to PSCfiles..."<<endl;
	  break;
	}
	
	time(&start);
      }

      if(ievt>=Nevts) continue;
      fHits = fEventHits->at(ievt)->GetfHits();
      sourcepos = fEventHits->at(ievt)->GetSourcePos();
      EGamma = fEventHits->at(ievt)->GetSourceE();
      ievt++;
    }

    // calc calpos from HCs
    for(Hit* ahit : *fHits) ahit->CalcAveHCsPosition();

    // tracking
    Tracker tracker(fHits, EGamma, sourcepos);
    tracker.OFTtracking();
    //tracker.Simpletracking();
    vector<int> atrack = tracker.GetTrack();

    // make paths
    if(atrack.size()>1){ // at least two hits
      double incE = EGamma;
      double depE = fHits->at(atrack[0])->GetE(); // keV

      Hit *sourcehit = new Hit(sourcepos);
#ifdef NTHREADS2
      lock_guard<mutex> lock(Pathsmtx);
#endif
      if(fHits->at(atrack[0])->hasHitCollection()>0 ||
	 fHits->at(atrack[1])->hasHitCollection()>0){ // change to at least one HC
	Path *apath = new Path(sourcehit,fHits->at(atrack[0]),fHits->at(atrack[1]),
			       incE, depE, incE, depE);
	fPaths->push_back(apath);
	cPaths++;
      }
      incE = incE - depE;
      for(int i=1; i<atrack.size()-1; i++){ //<--- was atrack.size()-2 ????
	depE = fHits->at(atrack[i])->GetE(); // keV
	if(fHits->at(atrack[i-1])->hasHitCollection()>0 ||
	   fHits->at(atrack[i])->hasHitCollection()>0 ||
	   fHits->at(atrack[i+1])->hasHitCollection()>0){ // change to at least one HC
	  Path *apath = new Path(fHits->at(atrack[i-1]),fHits->at(atrack[i]),fHits->at(atrack[i+1]),
				 incE, depE, incE, depE);
	  fPaths->push_back(apath);
	  cPaths++;
	}
	incE = incE - depE;
      }
    }
    
  }

  return;
}


void AGATA::Tracking(){
  ievt = 0;

  cout<<Form("Mem %.2fGB...",GetCurrentMemoryUsage()/GB);
  cout<<"clear paths in HCs...";
  for(HitCollection* hc : *fAllHCs){
    hc->GetPaths()->clear();
    hc->GetPaths()->shrink_to_fit();
    hc->GetPathHits()->clear();
    hc->GetPathHits()->shrink_to_fit();
  }
  cout<<"clear fPaths...";
  for(Path* apath : *fPaths) delete apath;
  fPaths->clear();
  fPaths->shrink_to_fit();
  cout<<Form("Mem %.2fGB",GetCurrentMemoryUsage()/GB)<<endl;

  cPaths = fPaths->size();
  
#ifndef NTHREADS2
  TrackingLoop();

#else
  thread th[NTHREADS2];
  cout<<"using "<<NTHREADS2<<" threads:"<<endl;
  for(int i=0; i<NTHREADS2; i++){
    th[i] = thread(&AGATA::TrackingLoop, this);
  }

  for(int i=0; i<NTHREADS2; i++){
    if(th[i].joinable())
      th[i].join();
  }
  
#endif
  time(&stop);
  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  double MemUsage = MemUsageGB / MemTotalGB * 100;
  cout<<"\r finish read "<<ievt<<" / "<<fEventHits->size()<<" evts"
      <<Form("(%.0fs/kevts)..",difftime(stop,start))
      <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
      <<"Hits-"<<cHits
      <<" HCs-"<<cHCs<<".."
      <<"Paths-"<<cPaths<<endl;
  
  return;
}


void AGATA::RegisterPathswithHCs(){
  cout<<"\r clear paths in HCs..."<<flush;
  for(HitCollection* hc : *fAllHCs){ hc->GetPaths()->clear();  hc->GetPathHits()->clear();}

  cout<<"\r register paths with HCs ..."<<flush;
  for(Path* apath : *fPaths){
    apath->RegisterWithHCs();
  }
  cout<<"\r register paths with HCs finished"<<endl;

  return;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// optimize HCs position
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
thread_local HitCollection *aHC;
thread_local AGATAgeo *ageo;

Double_t WrappedEstimator(const double *par){
  float fitpos[3];  for(int i=0; i<3; i++) fitpos[i] = par[i];
  
  aHC->SetFitPosition(fitpos); // give position to the HC
  vector<Path*>* paths = aHC->GetPaths();
  vector<Hit*>* phits = aHC->GetPathHits();
  
  Double_t result = 0;
  int Ncounts = 0;
  for(int i=0; i<paths->size(); i++){
    //phits->at(i)->CalcAveHCsPosition(); // update hits position in paths
    double change = paths->at(i)->CalcChi2(phits->at(i));
    result += change;
    Ncounts++;
  }

  if(Ncounts>0) result = result/Ncounts;

  if(result==0){ cout<<endl<<"WrappedEstimator result = 0!!!"<<endl;}

  // check bounds
  bool OutOfBounds = !(ageo->CheckBounds(aHC->GetDet(), aHC->GetSeg(), par));
  if(OutOfBounds) result = result*1000;
  
  return result;
}


void MinuitFCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  f = WrappedEstimator(par);
  return;
}


void AGATA::ExecFitLoop(int it){
  int Nhcs = fAllHCs->size();
  ageo = agatageo;

  //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad"); // static method, not work with multi threads

#ifdef MINUIT2
  ROOT::Minuit2::Minuit2Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer();
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000); // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(0);

  // create function wrapper for minimizer
  ROOT::Math::Functor f(&WrappedEstimator, 3); // minize function, 3 parameters x,y,z
  min->SetFunction(f);
  min->SetMaxFunctionCalls(200);

#else // use TMinuit
  TMinuit* min = new TMinuit(3); // 3 parameters x,y,z
  min->SetMaxIterations(10000);
  min->SetPrintLevel(-1);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 200; // max function calls
  arglist[1] = 0.001; // tolerance  

  min->SetFCN(MinuitFCN);
#endif
  
  for(;ihc<Nhcs;){ // loop HC of a segment

    // get HitCollection aHC
    {
#ifdef NTHREADS2
      lock_guard<mutex> lock(AllHCmtx); // lock fAllHCs
#endif

      if(ihc%10000==0){
	time(&stop);
	double MemUsageGB = GetCurrentMemoryUsage()/GB;
	double MemTotalGB = GetTotalSystemMemory()/GB;
	double MemUsage = MemUsageGB / MemTotalGB * 100;
	
	cout<<"\r Fit "<<it<<" : "
	    <<ihc<<" / "<<Nhcs<<" HCs"
	    <<Form("(%.0fs/kevts)..",difftime(stop,start))
	    <<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	    <<"Hits-"<<cHits
	    <<" HCs-"<<cHCs<<".."
	    <<"Paths-"<<cPaths<<flush;

	if(MemUsage>MaxMemUsage){
	  cout<<endl<<"exceed memory limit. Write to PSCfiles..."<<endl;
	  break;
	}
	
	time(&start);
      }

      if(ihc>=Nhcs) continue;
      aHC = fAllHCs->at(ihc);
      ihc++;
    }

    // optimize aHC position
    if(aHC->GetPaths()->size()<4) continue;
    if(kPSA && aHC->GetPaths()->size()<50) continue; // if PSA used for initial hit pos, then only optimize for HC with npaths>50
    
    // Set the free variables to be minimized!
    TVector3 initValue = aHC->GetPosition();
    TMatrixD SegPos = ageo->GetSegPos(aHC->GetDet(), aHC->GetSeg());
    TVector3 segpos(SegPos(0,0), SegPos(1,0), SegPos(2,0));
    TVector3 tmpvec = initValue - segpos;
    // avoid init value at boundary
    if(tmpvec.Mag()>2) tmpvec.SetMag(tmpvec.Mag()-1);
    else               tmpvec.SetMag(tmpvec.Mag()/2);
    initValue = segpos + tmpvec;

    Float_t fitpos[3];
    //fitpos[0]=initValue.X(); fitpos[1]=initValue.Y(); fitpos[2]=initValue.Z();
    
#ifdef MINUIT2
    min->SetLimitedVariable(0, "x", initValue.X(), 0.01,
			    initValue.X() - fitlimit, initValue.X() + fitlimit);
    min->SetLimitedVariable(1, "y", initValue.Y(), 0.01,
			    initValue.Y() - fitlimit, initValue.Y() + fitlimit);
    min->SetLimitedVariable(2, "z", initValue.Z(), 0.01,
			    initValue.Z() - fitlimit, initValue.Z() + fitlimit);

    // do the minimization
    //cout<<endl<<"pathsize "<<aHC->GetPaths()->size()<<" call minimize..."<<endl;
    min->Minimize();
    
    // get fitted parameters
    const double *xs = min->X();
    fitpos[0]=xs[0]; fitpos[1]=xs[1]; fitpos[2]=xs[2];

#else // use TMinuit
    min->mnparm(0, "x", initValue.X(), 0.01,
		initValue.X() - fitlimit, initValue.X() + fitlimit, ierflg);
    min->mnparm(1, "y", initValue.Y(), 0.01,
		initValue.Y() - fitlimit, initValue.Y() + fitlimit, ierflg);
    min->mnparm(2, "z", initValue.Z(), 0.01,
		initValue.Z() - fitlimit, initValue.Z() + fitlimit, ierflg);
    
    // do the minimization
    min->mnexcm("MIGRAD", arglist, 2, ierflg); // minimize with MIGRAD
    
    // get fitted parameters
    Double_t val, err;
    for(int ix=0; ix<3; ix++){
      min->GetParameter(ix,val,err);
      fitpos[ix] = val;
    }
#endif

    aHC->SetFitPosition(fitpos);
    aHC->SetPosition(fitpos);
    vector<Hit*>* phits = aHC->GetPathHits();
    for(Hit* ah : *phits) ah->CalcAveHCsPosition(); // update hits pos linked with aHC
  }

  return;
}


void AGATA::ExecFit(int repeat){

#ifdef MINUIT2
  cout<<"Minuit2 ";
#else
  cout<<"TMinuit ";
#endif
  
#ifdef NTHREADS2
  cout<<"using "<<NTHREADS2<<" threads:"<<endl;
#endif

  cout<<"Fit repeat "<<repeat<<" times.."<<endl;
  for(int it=0; it<repeat; it++){// loop repeat
    ihc = 0;

    time_t fitsta, fitsto;
    time(&fitsta);
    
#ifndef NTHREADS2
    ExecFitLoop(it);

#else
    thread th[NTHREADS2];
    for(int i=0; i<NTHREADS2; i++){
      th[i] = thread(&AGATA::ExecFitLoop, this, it);
    }

    for(int i=0; i<NTHREADS2; i++){
      if(th[i].joinable())
	th[i].join();
    }
#endif

    time(&fitsto);
    time(&stop);
    double MemUsageGB = GetCurrentMemoryUsage()/GB;
    double MemTotalGB = GetTotalSystemMemory()/GB;
    double MemUsage = MemUsageGB / MemTotalGB * 100;
    cout<<"\r Fit "<<it<<" : "
	<<ihc<<" / "<<fAllHCs->size()<<" HCs"
	<<Form("(%.0fs/kevts)..",difftime(stop,start))
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	<<"Hits-"<<cHits
	<<" HCs-"<<cHCs<<".."
	<<"Paths-"<<cPaths<<".."
	<<Form("%.1fsec",difftime(fitsto,fitsta))<<endl;

    
    // save fit result for det0000
    string PSCPath = "PSCfiles";
    if(kWithPS)
      if(kGroupPos) PSCPath = "PSCfiles_GP";
      else          PSCPath = "PSCfiles";
    else            PSCPath = "noPS";

    WritePSCfiles(0);

    gROOT->ProcessLine(Form(".!cp -fpdr %s/Det0000.root %s/iter/it/Det0000_fit%d.root",PSCPath.c_str(),PSCPath.c_str(),it));
    
  }// end of loop repeat
  

}


//**************************************************//
// PSA to assign initial pos
//**************************************************//
void AGATA::ReadPSAbasis(){

  cout<<"\e[1;31m Read basis for PSA: \e[0m"<<endl;
  string dbfile[3] = {"PSBase/pulseA.root",
		      "PSBase/pulseB.root",
		      "PSBase/pulseC.root"};

  for(int itype=0; itype<NType; itype++){

    if(Detid>-1 && itype!=Detid%3) continue;
    
    TFile *fdb = new TFile(dbfile[itype].c_str());
    if(!fdb->IsOpen()){
      cerr<<"cannot find dbfile "<<dbfile[itype]<<endl;
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

    int ncounter = 0;
    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);

      dbsegi = dbsegi-1; // start from 0

      PSAbasis apsabasis;
      for(int ix=0; ix<3; ix++) apsabasis.pos[ix] = dbposi[ix];
      
      bool skip = true;
      // pulse shape for comparison
      int fseg[NSeg_comp]; //0,1:fired seg, core; 2,3:next sectors; 4,5:next slice
      agatageo->GetNextSegs(dbsegi, fseg);

      for(int i=0; i<NSeg_comp; i++){
	int iseg = fseg[i];
	for(int isig=0; isig<NSig; isig++){
	  float tmpamp = 0;
	  if(iseg<NSeg) tmpamp = dbspulsei[iseg*NSig+isig];
	  else          tmpamp = dbcorei[isig];

	  if(tmpamp>0.5) skip=false;
	  apsabasis.spulse[i][isig] = tmpamp;
	}
      }

      if(skip) continue;

      ncounter++;
      fPSAbasis[itype][dbsegi].push_back(apsabasis);
    }

    cout<<"load "<<ncounter<<" points from "<<dbfile[itype]<<endl;
    fdb->Close();
  }

  return;
}


TVector3 AGATA::GetPSpos(int detid, int segid, PS *aps){

  int itype = detid%3;
  int ipos = -1;

#ifdef PSA
#ifdef WITHPS
  if(kPSA){

    int npoint = fPSAbasis[itype][segid].size();
    float minchi2 = 1e9;
    float asegpulse[NSig_comp], bsegpulse[NSig_comp];
    // compare fired seg, core and neighbor segment
    for(int ipoint=0; ipoint<npoint; ipoint++){

      float chi2=0;
      for(int iseg=0; iseg<NSeg_comp; iseg++){
	copy_n( aps->apulse[iseg], NSig_comp, asegpulse);
	copy_n(fPSAbasis[itype][segid][ipoint].spulse[iseg], NSig_comp, bsegpulse);
	float tmpchi2 = Chi2seg(asegpulse, bsegpulse);
	chi2 += tmpchi2; // sum

	if(chi2>minchi2) break; // interrupt if exceed minchi2
      }

      if(chi2<minchi2){
	minchi2 = chi2;
	ipos = ipoint;
      }

    }

  }
#endif
#endif

  TVector3 pos;
  if(ipos>0){
    // initial pos from PSA
    TMatrixD DetPos(3,1);
    for(int ix=0; ix<3; ix++) DetPos(ix,0) = fPSAbasis[itype][segid][ipos].pos[ix];
    TMatrixD LabPos = agatageo->Det2LabPos(detid,DetPos);
    pos.SetXYZ(LabPos(0,0), LabPos(1,0), LabPos(2,0));
  }else{
    // initial pos at segment center
    TMatrixD SegPos = agatageo->GetSegPos(detid,segid);
    pos.SetXYZ(SegPos(0,0), SegPos(1,0), SegPos(2,0));
  }
  
  return pos;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// tree and file
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
void AGATA::InitTreeWrite(TTree *tree){
  tree->Branch("det",&det);
  tree->Branch("seg",&seg);
  tree->Branch("index",&index);
  tree->Branch("nhits",&nhits);
  tree->Branch("Marker",&Marker);
  tree->Branch("chi2limit",chi2limit,"chi2limit[3]/F");
  tree->Branch("calpos",calpos,"calpos[3]/F");
  tree->Branch("cadpos",cadpos,"cadpos[3]/F");
  tree->Branch("calpos2",calpos2,"calpos2[3]/F");
  tree->Branch("cadpos2",cadpos2,"cadpos2[3]/F");
#ifdef REALPOS
  tree->Branch("labpos",labpos,"labpos[3]/F");
  tree->Branch("detpos",detpos,"detpos[3]/F");
  tree->Branch("dist",&dist);
  tree->Branch("dist2",&dist2);
#endif
  if(kWithPS){
    tree->Branch("cpos",cpos,"cpos[3]/F");
    tree->Branch("cpulsehits",&cpulsehits);
    tree->Branch("cpulse",cpulse,Form("cpulse[%d][%d]/F",NSeg_comp,NSig));  
    tree->Branch("segwgt",segwgt,Form("segwgt[%d]/F",NSeg_comp));  
    tree->Branch("spulse",spulse,Form("spulse[%d][%d]/F",NSegCore,NSig));  
  }
  tree->Branch("npaths",&npaths);
}

void AGATA::InitTreeRead(TTree *tree){
  tree->SetBranchAddress("det",&det);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("index",&index);
  tree->SetBranchAddress("nhits",&nhits);
  tree->SetBranchAddress("chi2limit",chi2limit);
  tree->SetBranchAddress("calpos",calpos);
  tree->SetBranchAddress("cadpos",cadpos);
  tree->SetBranchAddress("calpos2",calpos2);
  tree->SetBranchAddress("cadpos2",cadpos2);
#ifdef REALPOS
  tree->SetBranchAddress("labpos",labpos);
  tree->SetBranchAddress("detpos",detpos);
#endif
  if(kWithPS){
    tree->SetBranchAddress("cpos",cpos);
    tree->SetBranchAddress("cpulse",cpulse);
    tree->SetBranchAddress("cpulsehits",&cpulsehits);
    tree->SetBranchAddress("segwgt",segwgt);
    tree->SetBranchAddress("spulse",spulse);
  }
  tree->SetBranchAddress("npaths",&npaths);
}

void AGATA::ClosePSCFiles(){
  for(int idet=0; idet<NDets; idet++){
    if(!pscfile[idet]) continue;
    if(pscfile[idet]->IsOpen()) pscfile[idet]->Close();
  }
}



//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// check memory
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
Double_t AGATA::GetTotalSystemMemory(){
  Long64_t pages = sysconf(_SC_PHYS_PAGES);
  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) pages)*page_size;
}

Double_t AGATA::GetCurrentMemoryUsage(){
  Int_t tSize, resident, share;
  ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) resident) * page_size;
}

#endif // #ifndef AGATA_CC
