#ifndef AGATA_CC
#define AGATA_CC

#include <numeric>
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
  maxndiv   = 0;

  for(int idet=0; idet<MaxNDets; idet++){
    SkipDet[idet] = false;
    for(int iseg=0; iseg<NSeg; iseg++){
      fPSC[idet][iseg] = new vector<PSC*>();
      fHCs[idet][iseg] = new vector<HitCollection*>();
    }
  }
  fAllHCs = new vector<HitCollection*>();
  
  fEventHits = new vector<EventHits*>();
  fPaths = new vector<Path*>();

  // parameters for HC pos optimize
  fitlimit = 5;
  // PSC number limit
  for(int idet=0; idet<MaxNDets; idet++){
    if(SkipDet[idet]) continue;

    if(Detid<0) PSClimit[idet] = 2500;
    else        PSClimit[idet] = 5000;
    if(idet==0) PSClimit[idet] = 5000;
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
    if(SkipDet[idet]) continue;
    LoadPSCfiles(idet);
    WritePSCfiles(idet);
    cout<<endl;
  }
}

void AGATA::WritePSCfiles(int detid){ // create Pulse Shape Collection files

  //ClosePSCFiles();

  vector <int> idlist; // detid list to write
  for(int idet=0; idet<NDets; idet++){
    if(SkipDet[idet]) continue;
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

	if(nhits<1) continue;
	
	if( index != fPSC[idet][iseg]->at(ic)->index){
	  cout<<endl<<"det "<<idet<<" seg "<<iseg<<" ["<<ic<<"] index not match :"
	      <<index<<" != "<<fPSC[idet][iseg]->at(ic)->index<<endl;
	  continue;
	}

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

	for(int iiseg=0; iiseg<NSegCore; iiseg++){
	  copy_n(fPSC[idet][iseg]->at(ic)->spulse[iiseg], NSig, spulse[iiseg]);
	}
	copy_n(fPSC[idet][iseg]->at(ic)->devsigma, NSeg_comp, devsigma);

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
    if(SkipDet[idet]) continue;
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

	PSC* apsc = new PSC(det, seg);
	apsc->index = index;
	apsc->nhits = nhits;

	copy_n(labpos, 3, apsc->labpos);
	copy_n(detpos, 3, apsc->detpos);

	for(int iiseg=0; iiseg<NSegCore; iiseg++){
	  copy_n(spulse[iiseg], NSig, apsc->spulse[iiseg]);
	}
	copy_n(devsigma, NSeg_comp, apsc->devsigma);
	
	fPSC[idet][iseg]->push_back(apsc);
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
  for(int idet = 0; idet<NDets; idet++){
    if(SkipDet[idet]) continue;
    for(int iseg=0; iseg<NSeg; iseg++){
      cPSCmem-=fPSC[idet][iseg]->size();
      fPSC[idet][iseg]->clear();
      fPSC[idet][iseg]->shrink_to_fit();
    }
  }
  cPSCtotal = 0;
  cPSCmem   = 0;
  cPSCfile  = 0;

  return;
}

void AGATA::GetPSCstat(long long *PSCstat){
  PSCstat[0] = cPStotal;
  PSCstat[1] = cPSCtotal;
  PSCstat[2] = cPSCmem;
  PSCstat[3] = cPSCfile;
  PSCstat[4] = maxnhits;
  PSCstat[5] = cPaths;
  PSCstat[6] = cHits;
  PSCstat[7] = cHCs;
  PSCstat[8] = maxndiv;
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
    if(SkipDet[idet]) continue;
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
      hctree[iseg]->Branch("labpos",labpos,"labpos[3]/F"); // real position in lab frame

      for(HitCollection* ahc : *fHCs[idet][iseg]){

	if(ahc->GetSize()<1) continue;

	det = ahc->GetDet();
	seg = ahc->GetSeg();
	index = ahc->GetPid();
	TVector3 tmp;
	tmp = ahc->GetInitPosition();
	calpos[0] = tmp.X(); calpos[1] = tmp.Y(); calpos[2] = tmp.Z();

	tmp = ahc->GetPosition();
	calpos2[0] = tmp.X(); calpos2[1] = tmp.Y(); calpos2[2] = tmp.Z();

	tmp = ahc->GetRealPosition();
	labpos[0] = tmp.X(); labpos[1] = tmp.Y(); labpos[2] = tmp.Z();

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
    if(SkipDet[idet]) continue;
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
  hctree->SetBranchAddress("labpos",labpos); // real position in lab frame
  
  int Nhcs = hctree->GetEntriesFast();

  for(ihc=0; ihc<Nhcs; ihc++){
    if(ihc%10000==0) cout<<"\r load "<<ihc<<" / "<<Nhcs<<" HitCollections..."<<flush;
    hctree->GetEntry(ihc);

    //if(fHCs[det][seg]->size()!=index){ cerr<<"HC index not match!!!"<<endl;}

    HitCollection* ahc = new HitCollection(det, seg, index, labpos, calpos);

    fHCs[det][seg]->push_back(ahc);
    HCMap[det][seg].push_back(index);
    ahc->SetGid(ihc);
    fAllHCs->push_back(ahc);
    cHCs++;
  }

  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  double MemTotalGB = GetTotalSystemMemory()/GB;
  cout<<"\r load "<<ihc<<" / "<<Nhcs<<" HitCollections..."<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."<<endl;

  for(int idet : idlist){
    for(int iseg=0; iseg<NSeg; iseg++){
      HCstat[idet][iseg][0] = fHCs[idet][iseg]->size();
      HCstat[idet][iseg][1] = fHCs[idet][iseg]->back()->GetPid();

      cout<<"\r fHCs["<<idet<<"]["<<iseg<<"]->"
	  <<"at("<<HCstat[idet][iseg][0]-1<<")->Pid = "<<HCstat[idet][iseg][1]<<flush;

    }
  }
  cout<<endl;

  return;
}


void AGATA::ClearHCMem(){
  for(int idet = 0; idet<NDets; idet++){
    if(SkipDet[idet]) continue;
    for(int iseg=0; iseg<NSeg; iseg++){
      fHCs[idet][iseg]->clear();
      HCMap[idet][iseg].clear();
      fHCs[idet][iseg]->shrink_to_fit();
      HCstat[idet][iseg][0] = 0;
      HCstat[idet][iseg][1] = -1;
    }
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
    if(SkipDet[idet]) continue;
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
  
  int            interid; // interaction id in a event
  vector<int>    vhcid;
  float          depE;
#ifdef NOISE
  int            noiseidx;
  int            noiseidxshift;
#endif
  int            nsource;
  vector<float>  sourcex;
  vector<float>  sourcey;
  vector<float>  sourcez;
  vector<float>  sourceeng;
  int            iclust;
  int            bestis;
#ifdef DIFFTOTE
  float          Etot;
#endif

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
	htree[idet]->Branch("interid",&interid);
	htree[idet]->Branch("hcid",&vhcid);
	htree[idet]->Branch("depE",&depE);
	htree[idet]->Branch("calpos",calpos,"calpos[3]/F");
	htree[idet]->Branch("labpos",labpos,"labpos[3]/F");
#ifdef NOISE
	htree[idet]->Branch("noiseidx",&noiseidx);
	htree[idet]->Branch("noiseidxshift",&noiseidxshift);
#endif
	htree[idet]->Branch("nsource",&nsource);
	htree[idet]->Branch("sourcex",&sourcex);
	htree[idet]->Branch("sourcey",&sourcey);
	htree[idet]->Branch("sourcez",&sourcez);
	htree[idet]->Branch("sourceeng",&sourceeng);
	htree[idet]->Branch("iclust",&iclust);
	htree[idet]->Branch("bestis",&bestis);
#ifdef DIFFTOTE
	htree[idet]->Branch("Etot",&Etot);
#endif
      }
    }

    sourcex.clear(); sourcex.shrink_to_fit();
    sourcey.clear(); sourcey.shrink_to_fit();
    sourcez.clear(); sourcez.shrink_to_fit();
    sourceeng.clear(); sourceeng.shrink_to_fit();

    nsource = fEvent->GetNSource();
    for(int is=0; is<nsource; is++){
      TVector3 SourcePos = fEvent->GetSourcePos(is);
      sourcex.push_back(SourcePos.X());
      sourcey.push_back(SourcePos.Y());
      sourcez.push_back(SourcePos.Z());
      float tmpeng = fEvent->GetSourceE(is);
      sourceeng.push_back(tmpeng);
    }
    bestis = -1;

#ifdef DIFFTOTE
    Etot = fEvent->Etot;
#endif
    
    // get fHit information
    vector<Hit*>* fHits = fEvent->GetfHits();
    //for(Hit *ah : *fHits){
    for(int i=0; i<fHits->size(); i++){
      Hit *ah = fHits->at(i);

      det = ah->GetDet();
      seg = ah->GetSeg();
      interid = ah->GetInterid();
      
      bool kfill = false;
      for(int idet : idlist) if(idet==det) kfill = true;
      if(!kfill) continue;

      iclust = fEvent->GetClust(i);
      if(iclust>-1) bestis = fEvent->GetBestis(iclust);
      
      vhcid.clear(); vhcid.shrink_to_fit();
      vector<HitCollection*>* hcs = ah->GetHitCollections();
      for(HitCollection* ahc : *hcs) vhcid.push_back(ahc->GetPid());

      depE = ah->GetE(); // keV

      TVector3 tmp;
      tmp = ah->GetPosition();
      calpos[0] = tmp.X(); calpos[1] = tmp.Y(); calpos[2] = tmp.Z();

      tmp = ah->GetRealPosition();
      labpos[0] = tmp.X(); labpos[1] = tmp.Y(); labpos[2] = tmp.Z();

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
  long long nevt;
  while(!fin.eof()){
    fin.getline(buffer,500);
    if(strncmp(buffer,"#input",6)==0){
      int nsource = 0;
      fin >> buffer >> nsource;
      for(int i=0; i<nsource; i++){
	fin >> buffer >> sE >> spos[0] >> spos[1] >> spos[2];
      }
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
    int          iinterid;
    vector<int> *ihcid = 0;
    float        idepE;
    float        icalpos[3];
    float        ilabpos[3];

#ifdef NOISE
    int          inoiseidx;
    int          inoiseidxshift;
#endif
    int            insource;
    vector<float> *isourcex = 0;
    vector<float> *isourcey = 0;
    vector<float> *isourcez = 0;
    vector<float> *isourceeng = 0;
    int            iiclust;
    int            ibestis;
#ifdef DIFFTOTE
    float          iEtot;
#endif
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
	if(SkipDet[detid]) continue;
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
	if(SkipDet[detid]) continue;

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
	htree[detid]->SetBranchAddress("interid",       &obj[detid].iinterid);
	htree[detid]->SetBranchAddress("hcid",          &obj[detid].ihcid);
	htree[detid]->SetBranchAddress("depE",          &obj[detid].idepE);
	htree[detid]->SetBranchAddress("calpos",         obj[detid].icalpos);
	htree[detid]->SetBranchAddress("labpos",         obj[detid].ilabpos);
#ifdef NOISE
	htree[detid]->SetBranchAddress("noiseidx",      &obj[detid].inoiseidx);
	htree[detid]->SetBranchAddress("noiseidxshift", &obj[detid].inoiseidxshift);
#endif
	htree[detid]->SetBranchAddress("nsource",       &obj[detid].insource);
	htree[detid]->SetBranchAddress("sourcex",       &obj[detid].isourcex);
	htree[detid]->SetBranchAddress("sourcey",       &obj[detid].isourcey);
	htree[detid]->SetBranchAddress("sourcez",       &obj[detid].isourcez);
	htree[detid]->SetBranchAddress("sourceeng",     &obj[detid].isourceeng);
	htree[detid]->SetBranchAddress("iclust",        &obj[detid].iiclust);
	htree[detid]->SetBranchAddress("bestis",        &obj[detid].ibestis);
#ifdef DIFFTOTE
	htree[detid]->SetBranchAddress("Etot",          &obj[detid].iEtot);
#endif

	htree[detid]->GetEntry(Nevts[detid]-1);
	if(obj[detid].ientry>MaxEntry) MaxEntry=obj[detid].ientry;

	htree[detid]->GetEntry(idx[detid]);
      }

      MaxEntry += 1;

      hfileall->cd();
      int oconfig, orun, oentry, odet, oseg, ointerid;
      vector<int>    ovhcid;
      float          odepE;
      float          ocalpos[3];
      float          olabpos[3];
#ifdef NOISE
      int            onoiseidx;
      int            onoiseidxshift;
#endif
      int            onsource;
      vector<float>  osourcex;
      vector<float>  osourcey;
      vector<float>  osourcez;
      vector<float>  osourceeng;
      int            oiclust;
      int            obestis;
#ifdef DIFFTOTE
      float          oEtot;
#endif
      htreeall = new TTree("tree",Form("tree for Hits in config%d run%d alldet",iconfig,irun));
      htreeall->Branch("iconfig",&oconfig);
      htreeall->Branch("irun",&orun);
      htreeall->Branch("ientry",&oentry);
      htreeall->Branch("det",&odet);
      htreeall->Branch("seg",&oseg);
      htreeall->Branch("interid",&ointerid);
      htreeall->Branch("hcid",&ovhcid);
      htreeall->Branch("depE",&odepE);
      htreeall->Branch("calpos",ocalpos,"calpos[3]/F");
      htreeall->Branch("labpos",olabpos,"labpos[3]/F");
#ifdef NOISE
      htreeall->Branch("noiseidx",&onoiseidx);
      htreeall->Branch("noiseidxshift",&onoiseidxshift);
#endif
      htreeall->Branch("nsource",&onsource);
      htreeall->Branch("sourcex",&osourcex);
      htreeall->Branch("sourcey",&osourcey);
      htreeall->Branch("sourcez",&osourcez);
      htreeall->Branch("sourceeng",&osourceeng);
      htreeall->Branch("iclust",&oiclust);
      htreeall->Branch("bestis",&obestis);
#ifdef DIFFTOTE
      htreeall->Branch("Etot",&oEtot);
#endif

      for(int iety=0; iety<MaxEntry; iety++){ // loop entries

	// find hits belong to the same entry
	for(int detid=0; detid<NDets; detid++){
	  if(SkipDet[detid]) continue;

	  if(idx[detid]<0) continue;

	  while(obj[detid].ientry<=iety){

	    if(iety==obj[detid].ientry){
	      oconfig = obj[detid].iconfig;
	      orun = obj[detid].irun;
	      oentry = obj[detid].ientry;
	      odet = obj[detid].idet;
	      oseg = obj[detid].iseg;
	      ointerid = obj[detid].iinterid;
	      for(int id : *obj[detid].ihcid) ovhcid.push_back(id);
	      odepE = obj[detid].idepE;
	      for(int ix=0; ix<3; ix++) ocalpos[ix] = obj[detid].icalpos[ix];
	      for(int ix=0; ix<3; ix++) olabpos[ix] = obj[detid].ilabpos[ix];
#ifdef NOISE
	      onoiseidx = obj[detid].inoiseidx;
	      onoiseidxshift = obj[detid].inoiseidxshift;
#endif
	      onsource = obj[detid].insource;
	      for(float x : *obj[detid].isourcex) osourcex.push_back(x);
	      for(float y : *obj[detid].isourcey) osourcey.push_back(y);
	      for(float z : *obj[detid].isourcez) osourcez.push_back(z);
	      for(float eng : *obj[detid].isourceeng) osourceeng.push_back(eng);
	      oiclust = obj[detid].iiclust;
	      obestis = obj[detid].ibestis;
#ifdef DIFFTOTE
	      oEtot = obj[detid].iEtot;
#endif
	      
	      htreeall->Fill();
	      ovhcid.clear();
	      osourcex.clear();
	      osourcey.clear();
	      osourcez.clear();
	      osourceeng.clear();
	    }

	    // move to next idx
	    idx[detid]++;
	    if(idx[detid]>=Nevts[detid]) break;
	    htree[detid]->GetEntry(idx[detid]);
	    
	  }
	  
	}

      } // end of loop entries

      for(int detid=0; detid<NDets; detid++){
	if(SkipDet[detid]) continue;
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
    int          iinterid;
    vector<int> *ihcid = 0;
    float        idepE;
    float        icalpos[3];
    float        ilabpos[3];
#ifdef NOISE
    int          inoiseidx;
    int          inoiseidxshift;
#endif
    int            insource;
    vector<float> *isourcex = 0;
    vector<float> *isourcey = 0;
    vector<float> *isourcez = 0;
    vector<float> *isourceeng = 0;
    int            iiclust;
    int            ibestis;
#ifdef DIFFTOTE
    float          iEtot;
#endif
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
    double MemUsageGB = GetCurrentMemoryUsage()/GB;
    double MemTotalGB = GetTotalSystemMemory()/GB;
    cout<<"\r \e[1;33m Load Hits from "<<hfilename0<<"... "
	<<"Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."<<"\e[0m"<<flush;

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
    htree->SetBranchAddress("interid",       &obj.iinterid);
    htree->SetBranchAddress("hcid",          &obj.ihcid);
    htree->SetBranchAddress("depE",          &obj.idepE);
    htree->SetBranchAddress("calpos",         obj.icalpos);
    htree->SetBranchAddress("labpos",         obj.ilabpos);
#ifdef NOISE
    htree->SetBranchAddress("noiseidx",      &obj.inoiseidx);
    htree->SetBranchAddress("noiseidxshift", &obj.inoiseidxshift);
#endif
    htree->SetBranchAddress("nsource",       &obj.insource);
    htree->SetBranchAddress("sourcex",       &obj.isourcex);
    htree->SetBranchAddress("sourcey",       &obj.isourcey);
    htree->SetBranchAddress("sourcez",       &obj.isourcez);
    htree->SetBranchAddress("sourceeng",     &obj.isourceeng);
    htree->SetBranchAddress("iclust",        &obj.iiclust);
    htree->SetBranchAddress("bestis",        &obj.ibestis);
#ifdef DIFFTOTE
    htree->SetBranchAddress("Etot",          &obj.iEtot);
#endif

    htree->GetEntry(Nevts-1);
    MaxEntry=obj.ientry;

    MaxEntry += 1;
    if(MaxEntry>0) NevtsTotal += MaxEntry;

    vector<int>            vdet;
    vector<int>            vseg;
    vector<int>            vinterid;
    vector<vector<int>>    vhcid;
    vector<float>          vdepE;
    vector<vector<double>> vcalpos;
    vector<vector<double>> vlabpos;
#ifdef NOISE
    vector<int>            vnoiseidx;
    vector<int>            vnoiseidxshift;
#endif
    vector<int>            iclust;
    vector<int>            bestis;

    int                    nsource;
    vector<float>          sourcex;
    vector<float>          sourcey;
    vector<float>          sourcez;
    vector<float>          sourceeng;
#ifdef DIFFTOTE
    float                  Etot;
#endif

    idx = 0;
    for(int iety=0; iety<MaxEntry; iety++){ // loop entries
      //if(iety%10==0) cout<<"\r ientry = "<<iety<<flush;
      vdet.clear();
      vseg.clear();
      vinterid.clear();
      vhcid.clear();
      vdepE.clear();
      vcalpos.clear();
      vlabpos.clear();
#ifdef NOISE
      vnoiseidx.clear();
      vnoiseidxshift.clear();
#endif
      iclust.clear();
      bestis.clear();

      nsource = 0;
      sourcex.clear();
      sourcey.clear();
      sourcez.clear();
      sourceeng.clear();

      // find hits belong to the same entry
      for(; idx<Nevts; idx++){

	htree->GetEntry(idx);
	if(obj.ientry>iety) break;
	
	if(iety==obj.ientry){

	  vdet.push_back(obj.idet);
	  vseg.push_back(obj.iseg);
	  vinterid.push_back(obj.iinterid);
	  vhcid.push_back(*obj.ihcid);
	  vdepE.push_back(obj.idepE);

	  vector<double> tmppos1(3);
	  for(int ix=0; ix<3; ix++) tmppos1[ix]=obj.icalpos[ix];
	  vcalpos.push_back(tmppos1);

	  vector<double> tmppos2(3);
	  for(int ix=0; ix<3; ix++) tmppos2[ix]=obj.ilabpos[ix];
	  vlabpos.push_back(tmppos2);

#ifdef NOISE
	  vnoiseidx.push_back(obj.inoiseidx);
	  vnoiseidxshift.push_back(obj.inoiseidxshift);
#endif

	  iclust.push_back(obj.iiclust);
	  bestis.push_back(obj.ibestis);

	  if(nsource==0){
	    nsource = obj.insource;
	    for(float x : *obj.isourcex) sourcex.push_back(x);
	    for(float y : *obj.isourcey) sourcey.push_back(y);
	    for(float z : *obj.isourcez) sourcez.push_back(z);
	    for(float eng : *obj.isourceeng) sourceeng.push_back(eng);
#ifdef DIFFTOTE
	    Etot = obj.iEtot;
#endif
	  }
	}

      }

      if(vdet.size()<2) continue;

      // make EventHits
      vector<float> SourceE;
      vector<TVector3> SourcePos;
      for(int is=0; is<nsource; is++){
	SourceE.push_back(sourceeng[is]);
	SourcePos.push_back(TVector3(sourcex[is],sourcey[is],sourcez[is]));
      }
      EventHits* fEvent = new EventHits(SourceE, SourcePos);
      fEvent->SetIdx(iconfig,irun,iety);
#ifdef DIFFTOTE
      fEvent->Etot = Etot;
#endif
      
      for(int i=0; i<vdet.size(); i++){
	int detid = vdet[i];
	int segid = vseg[i];
	int interid = vinterid[i];

	float depE = vdepE[i]; // keV
	TVector3 initpos(vcalpos[i][0],vcalpos[i][1],vcalpos[i][2]); // mm
	TVector3 hitpos(vlabpos[i][0],vlabpos[i][1],vlabpos[i][2]); // mm

	Hit *ahit = new Hit(detid, segid, depE, hitpos, initpos);
	ahit->SetInterid(interid);
	
#ifdef NOISE
	ahit->SetNoiseIdx(vnoiseidx[i]);
	ahit->SetNoiseIdxShift(vnoiseidxshift[i]);
#endif
	fEvent->Add(ahit);
	fEvent->SignClust(iclust[i],i);
	if(iclust[i]>-1) fEvent->SetBestis(iclust[i],bestis[i]);

	cHits++;

	
#ifdef NTHREADS2
	lock_guard<mutex> lock(PSCmtx[detid][segid]);
#endif
	// connect with HCs
	for(int ii=0; ii<vhcid[i].size(); ii++){
	  int pscid = vhcid[i][ii];
	  if( pscid > HCstat[detid][segid][1] ){
	    cerr<<"pscid = "<<pscid<<" outside fHCs["<<detid<<"]["<<segid<<"] range!!!!"<<endl;
	    continue;
	  }

	  int ipsc = FindHC(detid, segid, pscid);
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
  cPathsN[0] = cPathsN[1] = cPathsN[2] = cPathsN[3] = 0;

  return;
}


void AGATA::SortPSC(){
  for(int idet=0; idet<NDets; idet++){
    if(SkipDet[idet]) continue;
    for(int iseg=0; iseg<NSeg; iseg++){
      // increasing ientry
      sort(fPSC[idet][iseg]->begin(),fPSC[idet][iseg]->end(),
	   [](const PSC* lhs, const PSC* rhs){
	     return lhs->index<rhs->index;});
    }
  }  
  return;
}


int AGATA::FindHC(int detid, int segid, int pscid){
  int ipsc = -1;
  int npsc = HCstat[detid][segid][0];

  int itmp = min(pscid, npsc-1);
  for(int i=0; i<npsc; i++, itmp--){

    if( unlikely( itmp < 0 ) ) itmp+=npsc;

    int tmpid = HCMap[detid][segid][itmp];
    if     ( tmpid >  pscid ) continue;
    else if( tmpid == pscid){ ipsc=itmp; break;}
    else if( tmpid <  pscid) break;
  }
  
  return ipsc;
}

//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// group pulse shape
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
int AGATA::InitPSCandHC(int detid, int segid){

  int ipsc = -1;
  if(freeHCs[detid][segid].size()>0){ // use old HCs
    ipsc = freeHCs[detid][segid][0];
    freeHCs[detid][segid].erase(freeHCs[detid][segid].begin());

  }else if(kAddNewPSC){ // create new HCs

    ipsc = fPSC[detid][segid]->size();

    // initial PSC for detid,segid
    PSC *newpsc = new PSC(detid, segid);
    newpsc->index = ipsc;
    newpsc->nhits = 0;
    newpsc->dividx = vector<int>(0);
    for(int ix=0; ix<3; ix++){
      newpsc->labpos[ix] = 0;
      newpsc->detpos[ix] = 0;
    }
    int fseg[NSeg_comp];
    agatageo->GetNextSegs(segid, fseg);
    for(int iseg=0; iseg<NSeg_comp; iseg++){
      newpsc->divzone[iseg] = -1;
      newpsc->segcmp[iseg] = fseg[iseg];
      newpsc->devabscut[iseg][0] = -1;
      newpsc->devsigma[iseg] = -1;
      newpsc->devabs[iseg] = vector<float>(0);
      newpsc->dev[iseg] = vector<float>(0);
    }
    for(int iseg=0; iseg<NSegCore; iseg++){
      for(int isig=0; isig<NSig; isig++){
	newpsc->spulse[iseg][isig] = 0;
      }
    }

    fPSC[detid][segid]->push_back(newpsc);

    // new HitCollection
    TMatrixD SegPos = agatageo->GetSegPos(detid,segid);
    float dummycalpos[3];
    float dummylabpos[3];
    for(int ix=0; ix<3; ix++){
      dummycalpos[ix] = SegPos(ix,0);
      dummylabpos[ix] = 0;
    }
    HitCollection* newhc = new HitCollection(detid, segid, ipsc, dummylabpos, dummycalpos);
    fHCs[detid][segid]->push_back(newhc);

#ifdef NTHREADS
    lock_guard<mutex> lock2(AllHCmtx); // lock fAllHCs
#endif
    newhc->SetGid(fAllHCs->size());
    fAllHCs->push_back(newhc);
  }
  
  if(ipsc>-1){
    cPSCtotal++;   cPSCmem++;    cHCs++;
  }
  
  return ipsc;
}


void AGATA::FindInitZone(PS *aps, vector<int> *initzone){
  int segid = aps->seg;
  int fseg[NSeg_comp];
  agatageo->GetNextSegs(segid, fseg);

  float apulse[NSig_comp], bpulse[NSig_comp];
  float asum, bsum;
  float apulseabs[NSig_comp], bpulseabs[NSig_comp];
  float asumabs, bsumabs;
  float ratio;

  vector<int> zone[3];

  // 0:iseg and 1:core
  copy_n( aps->opulse[fseg[0]], NSig_comp, apulse);
  copy_n( aps->opulse[fseg[1]], NSig_comp, bpulse);

  asum = std::accumulate( apulse, apulse+NSig_comp, 0.);
  bsum = std::accumulate( bpulse, bpulse+NSig_comp, 0.);

  ratio = asum/bsum;
  if( ratio<0.8 ) zone[0].push_back(0); // inner
  if( ratio>1.2 ) zone[0].push_back(1); // outer
  if( fabs(ratio-1)<=0.2 ) zone[0].push_back(2); // center

  // neighbor sector 2 & 3
  copy_n( aps->opulse[fseg[2]], NSig_comp, apulse);
  copy_n( aps->opulse[fseg[3]], NSig_comp, bpulse);

  for(int i=0; i<NSig_comp; i++){
    apulseabs[i] = fabs(apulse[i]);
    bpulseabs[i] = fabs(bpulse[i]);
  }

  asum = std::accumulate( apulse, apulse+NSig_comp, 0.);
  bsum = std::accumulate( bpulse, bpulse+NSig_comp, 0.);

  asumabs = std::accumulate( apulseabs, apulseabs+NSig_comp, 0.);
  bsumabs = std::accumulate( bpulseabs, bpulseabs+NSig_comp, 0.);

  ratio = asumabs/bsumabs;
  if( ratio<0.8 ) zone[1].push_back(0); // left
  if( ratio>1.2 ) zone[1].push_back(1); // right
  if( fabs(ratio-1)<=0.2) zone[1].push_back(2); // center

  if( asum<-0.5*asumabs ) zone[2].push_back(0); // negative
  if( asum> 0.5*asumabs ) zone[2].push_back(1); // positive
  if( fabs(asum)<=0.5*asumabs) zone[2].push_back(2); // bipolar

  // neighbor slice 4 & 5 <--- not used


  // init zone
  for(int i0=0; i0<zone[0].size(); i0++){
    for(int i1=0; i1<zone[1].size(); i1++){
      for(int i2=0; i2<zone[2].size(); i2++){
	int tmpzone = zone[0][i0] + zone[1][i1]*3 + zone[2][i2]*3*3;
	initzone->push_back(tmpzone);
      }
    }
  } 


  return;
}


void AGATA::FindInitPSC(Hit *ahit, vector<int> *initpsc){

  TVector3 HitPos = ahit->GetPosition(); // hit pos in lab frame from PSA

  int detid = ahit->GetDetectorID();
  int segid = ahit->GetSegmentID();

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  for(HitCollection* ahc : *fHCs[detid][segid]){
    TVector3 HCPos = ahc->GetInitPosition();
    float dist = (HitPos-HCPos).Mag();
    if(dist<RADIUS0) initpsc->push_back(ahc->GetPid());
  }

  return;
}


int AGATA::AddPS(PS *aps, Hit *ahit){
  vector<int> *initpsc = new vector<int>();

  bool findpos = false;
#ifdef PSA
  if(kPSA) findpos = true;
#endif

  if(findpos){
    FindInitPSC(ahit, initpsc);
  }else{
    FindInitZone(aps, initpsc); // find all init zone the PS belong to
  }

  int detid = aps->det;
  int segid = aps->seg;

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  if(findpos){
    if(initpsc->size()<1){
      int tmp = InitPSCandHC(detid, segid);
      TVector3 LabPos = ahit->GetPosition();
      Float_t hitpos[3];
      hitpos[0] = LabPos.X();  hitpos[1] = LabPos.Y();  hitpos[2] = LabPos.Z();
      fHCs[detid][segid]->at(tmp)->SetInitPosition(hitpos);
      initpsc->push_back(tmp);
    }
  }else{
    while(fPSC[detid][segid]->size()<27){
      int tmp = InitPSCandHC(detid, segid);
    }
  }

  cPStotal++;
  int tmp;
  for(int i=0; i<initpsc->size(); i++){
    tmp = AddPStoPSC(aps, ahit, initpsc->at(i));
  }
  return tmp;
}


int AGATA::AddPStoPSC(PS *aps, Hit *ahit, int ipsc){
  
  int detid = aps->det;
  int segid = aps->seg;

  int npsc = fPSC[detid][segid]->size();
  if( ipsc > npsc-1 ){
    if(ipsc==0){
      int tmp = InitPSCandHC(detid, segid);
      if(tmp<0) return tmp;
    }else{
      cout<<Form("cannot find ipsc = %d ; fPSC[%d][%d].size() = %d",ipsc,detid,segid,(int)fPSC[detid][segid]->size())<<endl;
      return 0;
    }
  }

  double nhitstmp = (double)fPSC[detid][segid]->at(ipsc)->nhits;

  // calc average position
  TMatrixD pos(3,1); for(int ix=0; ix<3; ix++) pos(ix,0) = aps->labpos[ix];
  TMatrixD LabPos(3,1);
  for(int ix=0; ix<3; ix++) LabPos(ix,0) = fPSC[detid][segid]->at(ipsc)->labpos[ix];
  LabPos = 1./(nhitstmp+1)*(nhitstmp*LabPos + pos);
  TMatrixD DetPos = agatageo->Lab2DetPos(detid,LabPos);

  for(int ix=0; ix<3; ix++){
    fPSC[detid][segid]->at(ipsc)->labpos[ix] = LabPos(ix,0);
    fPSC[detid][segid]->at(ipsc)->detpos[ix] = DetPos(ix,0);
  }

  // calc average pulse
  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++){
      fPSC[detid][segid]->at(ipsc)->spulse[iseg][isig] =
	1./(nhitstmp+1)*(fPSC[detid][segid]->at(ipsc)->spulse[iseg][isig]*nhitstmp +
			 aps->opulse[iseg][isig]);
    }
  }

  fPSC[detid][segid]->at(ipsc)->nhits++;

  // add to HitCollection
  fHCs[detid][segid]->at(ipsc)->AddHit(ahit);
  fHCs[detid][segid]->at(ipsc)->SetRealPosition(fPSC[detid][segid]->at(ipsc)->labpos);
  ahit->AddHitCollection(fHCs[detid][segid]->at(ipsc));

  int tmpnhits = fPSC[detid][segid]->at(ipsc)->nhits;
  if(tmpnhits>maxnhits) maxnhits=tmpnhits;
  
  return tmpnhits;
}


float AGATA::FindMaxDev(PS *aps, Hit *ahit){
  int detid = aps->det;
  int segid = aps->seg;
  float MaxDev = 0;
  
#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  vector<HitCollection*>* hcs = ahit->GetHitCollections();
  for(HitCollection* ahc : *hcs){
    int ipsc = ahc->GetPid();
    PSC *apsc = fPSC[detid][segid]->at(ipsc);

    for(int is=0; is<NSeg_comp; is++){

      if( DivDir>-1 && ((int)(is/2))!=DivDir) continue;

      int iseg = apsc->segcmp[is];

      float asegpulse[NSig_comp], bsegpulse[NSig_comp];
      copy_n( aps->opulse[iseg], NSig_comp, asegpulse);
      copy_n(apsc->spulse[iseg], NSig_comp, bsegpulse);

      float dev[4]; // 0: abs dev; 1: dev; 2: dev sum sum; 3: empty
      Devseg(asegpulse, bsegpulse, dev);
      apsc->devabs[is].push_back(dev[0]);

      if(dev[0]>MaxDev) MaxDev = dev[0];

    }
  }

  return MaxDev;
}


void AGATA::FindDevCut(){
  for(int idet=0; idet<MaxNDets; idet++){
    if(SkipDet[idet]) continue;
    for(int iseg=0; iseg<NSeg; iseg++){

      for(PSC *apsc : *fPSC[idet][iseg]){

	for(int is=0; is<NSeg_comp; is++){

	  if( DivDir>-1 && ((int)(is/2))!=DivDir) continue;

	  if(apsc->devabscut[is][0]<0){
	    sort(apsc->devabs[is].begin(), apsc->devabs[is].end(), [](const float& l, const float& r){return l<r;});
	    int nsize = apsc->devabs[is].size();
	    if(nsize<1) continue;

	    if(apsc->devabs[is][0] > (0.2*apsc->devabs[is][nsize-1])){
	      // empty center zone
	      apsc->devabscut[is][0] = 0.99*apsc->devabs[is][0];
	      apsc->devabscut[is][1] = 0.99*apsc->devabs[is][0];

	    }else{
	      int ncut;
	      ncut = (int)(0.5*nsize); // 0.5
	      apsc->devabscut[is][0] = apsc->devabs[is][ncut];

	      ncut = (int)(0.5*nsize); // 0.5
	      apsc->devabscut[is][1] = apsc->devabs[is][ncut];
	    }
	  }

	  apsc->devabs[is].clear();
	}
      }

    }
  }

  return;  
}


void AGATA::FindDevSigma(PS *aps, Hit *ahit){
  int detid = aps->det;
  int segid = aps->seg;
  
#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  vector<HitCollection*>* hcs = ahit->GetHitCollections();
  for(HitCollection* ahc : *hcs){
    int ipsc = ahc->GetPid();
    PSC *apsc = fPSC[detid][segid]->at(ipsc);

    for(int is=0; is<NSeg_comp; is++){
      int iseg = apsc->segcmp[is];

      float asegpulse[NSig_comp], bsegpulse[NSig_comp];
      copy_n( aps->opulse[iseg], NSig_comp, asegpulse);
      copy_n(apsc->spulse[iseg], NSig_comp, bsegpulse);

      float dev[4]; // 0: abs dev; 1: dev; 2: dev sum sum; 3: empty
      Devseg(asegpulse, bsegpulse, dev);
      apsc->dev[is].push_back(dev[1]);

    }
  }

  return;
}


void AGATA::CalcDevSigma(){
  for(int idet=0; idet<MaxNDets; idet++){
    if(SkipDet[idet]) continue;
    for(int iseg=0; iseg<NSeg; iseg++){

      for(PSC *apsc : *fPSC[idet][iseg]){

	for(int is=0; is<NSeg_comp; is++){
	  int nsize = apsc->dev[is].size();
	  if(nsize<1) continue;

	  double sum = accumulate(apsc->dev[is].begin(), apsc->dev[is].end(), 0.0);
	  double mean = sum / nsize;

	  vector<double> diff(nsize);
	  transform(apsc->dev[is].begin(), apsc->dev[is].end(), diff.begin(),
		    bind2nd(minus<double>(), mean));
	  double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	  double stdev = sqrt(sq_sum/nsize);
	  apsc->devsigma[is] = stdev;

	  apsc->dev[is].clear();
	}
      }

    }
  }

  return;  
}


void AGATA::FindDivZone(PS *aps, PSC *apsc, vector<vector<int>> *divzone){

  for(int is=0; is<NSeg_comp; is++){

    vector<int> tmp;

    if( DivDir<0 || ((int)(is/2))==DivDir){

      // divide
      int iseg = apsc->segcmp[is];
      float asegpulse[NSig_comp], bsegpulse[NSig_comp];
      copy_n( aps->opulse[iseg], NSig_comp, asegpulse);
      copy_n(apsc->spulse[iseg], NSig_comp, bsegpulse);

      float dev[4]; // 0: abs dev; 1: dev; 2: dev sum sum; 3: empty
      Devseg(asegpulse, bsegpulse, dev);

      if(apsc->devabscut[is][0]<0){
	cout<<"devabscut[0]<0 ; something wrong!!!"<<endl;
	return;
      }

      if(dev[0] <= apsc->devabscut[is][0]) tmp.push_back(0); // center zone
      if(dev[0] >  apsc->devabscut[is][1]){
	if( dev[1] >=  0.6 * dev[0]) tmp.push_back(1);    // above zone
	if( dev[1] <= -0.6 * dev[0]) tmp.push_back(2);    // below zone
	// bipolar zone
	if( dev[1]<0.7*dev[0] && dev[1]>-0.7*dev[0]){
	  //tmp.push_back(3);                               // bipolar zone
	  if(dev[2] >  0) tmp.push_back(3);               // bipolar zone +-
	  else            tmp.push_back(4);               // bipolar zone -+
	}
      }

    }else{ // not divide
      tmp.push_back(0);
    }
    
    divzone->push_back(tmp);
  }

  return;
}


int AGATA::AddPStoDiv(PS *aps, Hit *ahit){
  int detid = aps->det;
  int segid = aps->seg;

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  int tmpnhits = 0;
  vector<HitCollection*>* hcs = ahit->GetHitCollections();
  int nhcs = hcs->size();
  for(int iHC=0; iHC<nhcs; iHC++){ // loop HitCollections
    HitCollection* ahc = hcs->at(iHC);
    int ipsc = ahc->GetPid();
    PSC *apsc = fPSC[detid][segid]->at(ipsc);

    if(apsc->nhits > MAXHITS){ // if the PSC is large size, then divide it
      vector<vector<int>> *divzone = new vector<vector<int>>();
      FindDivZone(aps, apsc, divzone); // find all divided zone the PS belong to
      int tmpidx[NSeg_comp];
      int ndiv = 1;
      for(int is=0; is<NSeg_comp; is++){
	ndiv *= divzone->at(is).size();
	tmpidx[is] = 0;
      }

      int tmpzone[NSeg_comp];
      for(int idiv=0; idiv<ndiv; idiv++){ // loop all divided zone the PS belong to
	for(int is=0; is<NSeg_comp; is++) tmpzone[is] = divzone->at(is)[tmpidx[is]];

	int ndaughter = apsc->dividx.size();
	int idaughter;
	int jpsc;
	for(idaughter=0; idaughter<ndaughter; idaughter++){
	  jpsc = apsc->dividx[idaughter];
	  bool kfound = true;
	  for(int is=0; is<NSeg_comp; is++){
	    if(tmpzone[is] != fPSC[detid][segid]->at(jpsc)->divzone[is]){
	      kfound = false;
	      break;
	    }
	  }

	  if(kfound) break;
	}

	if(idaughter==ndaughter){ // cannot find the daughter PSC
	  jpsc = InitPSCandHC(detid, segid); // create new daughter PSC
	  if(jpsc<0) continue;
	  for(int is=0; is<NSeg_comp; is++) fPSC[detid][segid]->at(jpsc)->divzone[is] = tmpzone[is];
	  apsc->dividx.push_back(jpsc);
	  if(idaughter>maxndiv) maxndiv = idaughter;
	}

	// Add PS to PSC and HC
	int tmp = AddPStoPSC(aps, ahit, jpsc);
	if(tmp>tmpnhits) tmpnhits=tmp;

	// next zone
	tmpidx[0]++;
	for(int is=0; is<NSeg_comp-1; is++){
	  if(tmpidx[is]==divzone->at(is).size()){
	    tmpidx[is] = 0;
	    tmpidx[is+1]++;
	  }
	}
	if(tmpidx[NSeg_comp-1]==divzone->at(NSeg_comp-1).size()){
	  if(idiv!=ndiv-1) cout<<"something wrong!!!"<<endl;
	  break;
	}

      } // end of loop divided zone
    } // end of if the PSC is large size
  } // end of loop HitCollections
  
  return tmpnhits;
}


int AGATA::CheckPSinPSC(PS *aps, Hit *ahit){
  int detid = aps->det;
  int segid = aps->seg;

  int nremove = 0;
  vector<HitCollection*>* hcs = ahit->GetHitCollections();
  for(HitCollection* ahc : *hcs){ // loop HitCollections
    int ipsc = ahc->GetPid();
    PSC *apsc = fPSC[detid][segid]->at(ipsc);

    bool kreject = false;
    for(int is=0; is<NSeg_comp; is++){
      int iseg = apsc->segcmp[is];

      float asegpulse[NSig_comp], bsegpulse[NSig_comp];
      copy_n( aps->opulse[iseg], NSig_comp, asegpulse);
      copy_n(apsc->cpulse[is],   NSig_comp, bsegpulse);

      float dev[4]; // 0: abs dev; 1: dev; 2: dev sum sum; 3: empty
      Devseg(asegpulse, bsegpulse, dev);

      if( fabs(dev[1]) > nSigma*apsc->devsigma[is] ) kreject = true;
      
      if(kreject) break;
    }

    if(kreject){
      RemovePSfromPSC(aps, ahit, ipsc);
      nremove++;
    }
    
  } // end of loop HitCollections
  
  return nremove;
}


void AGATA::MakeCPulse(){
#ifdef NTHREADS
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++)
      PSCmtx[detid][segid].lock(); // lock PSC for one segment
  }
#endif

  cout<<"\e[1;31m Make cpulse for all PSCs ... \e[0m"<<endl;
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++){
      for(PSC* apsc : *fPSC[detid][segid]){

	for(int is=0; is<NSeg_comp; is++){
	  int iseg = apsc->segcmp[is];
	  copy_n(apsc->spulse[iseg], NSig, apsc->cpulse[is]);
	}

      }
    }
  }

#ifdef NTHREADS
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++)
      PSCmtx[detid][segid].unlock(); // unlock PSC for one segment
  }
#endif
  
  return;
}


void AGATA::RemoveMotherPSC(){
#ifdef NTHREADS
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++)
      PSCmtx[detid][segid].lock(); // lock PSC for one segment
  }
#endif

  cout<<"\e[1;31m Remove Mother PSC ... \e[0m"<<endl;
  int counter = 0;
  for(HitCollection* ahc : *fAllHCs){
    int detid = ahc->GetDet();
    int segid = ahc->GetSeg();
    int ipsc  = ahc->GetPid();
    if(fPSC[detid][segid]->at(ipsc)->dividx.size()>0){
      RemovePSC(ahc);
      counter++;
      if(counter%1000==0)
	cout<<"\r Remove "<<counter<<" mother PSC"<<flush;
    }
  }
  cout<<"\r Remove "<<counter<<" mother PSC"<<endl;

#ifdef NTHREADS
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++)
      PSCmtx[detid][segid].unlock(); // unlock PSC for one segment
  }
#endif

  return;
}

void AGATA::RemoveSmallPSC(int minhits){
#ifdef NTHREADS
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++)
      PSCmtx[detid][segid].lock(); // lock PSC for one segment
  }
#endif

  cout<<"\e[1;31m Remove Small (<"<<minhits<<") PSC  ... \e[0m"<<endl;
  int counter = 0;
  for(HitCollection* ahc : *fAllHCs){
    int detid = ahc->GetDet();
    int segid = ahc->GetSeg();
    int ipsc  = ahc->GetPid();
    int tmpnhits = fPSC[detid][segid]->at(ipsc)->nhits;
    if(tmpnhits>0 && tmpnhits<minhits){
      RemovePSC(ahc);
      counter++;
      if(counter%1000==0)
	cout<<"\r Remove "<<counter<<" small PSC"<<flush;
    }
  }
  cout<<"\r Remove "<<counter<<" small PSC"<<endl;

#ifdef NTHREADS
  for(int detid=0; detid<MaxNDets; detid++){
    if(SkipDet[detid]) continue;
    for(int segid=0; segid<NSeg; segid++)
      PSCmtx[detid][segid].unlock(); // unlock PSC for one segment
  }
#endif

  return;
}



void AGATA::RemovePSfromPSC(PS *aps, Hit *ahit, int jpsc){ // remove ahit from PSC(s)
  int detid = aps->det;
  int segid = aps->seg;

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif

  vector<HitCollection*>* hcs = ahit->GetHitCollections();

  TMatrixD pos(3,1); for(int ix=0; ix<3; ix++) pos(ix,0) = aps->labpos[ix];
  
  for(HitCollection* ahc : *hcs){

    int ipsc = ahc->GetPid();

    if(jpsc>-1 && ipsc!=jpsc) continue;

    double nhitstmp = (double)fPSC[detid][segid]->at(ipsc)->nhits;

    // calc average position
    TMatrixD LabPos(3,1);
    for(int it=0; it<3; it++) LabPos(it,0) = fPSC[detid][segid]->at(ipsc)->labpos[it];

    if(nhitstmp<1){
      cerr<<"something wrong, det = "<<detid<<" seg = "<<segid<<" ipsc = "<<ipsc
	  <<" fPSC.nhits = "<<fPSC[detid][segid]->at(ipsc)->nhits
	  <<" fHC.nhits = "<<fHCs[detid][segid]->at(ipsc)->GetSize()<<endl;
      continue;
    }
    double ftmp = 1;
    if(nhitstmp>1) ftmp = 1./(nhitstmp-1);

    LabPos = ftmp*(nhitstmp*LabPos - pos);
    TMatrixD DetPos = agatageo->Lab2DetPos(detid,LabPos);
      
    for(int it=0; it<3; it++){
      fPSC[detid][segid]->at(ipsc)->labpos[it] = LabPos(it,0);
      fPSC[detid][segid]->at(ipsc)->detpos[it] = DetPos(it,0);
    }
    
    // calc average pulse
    for(int iseg=0; iseg<NSegCore; iseg++)
      for(int isig=0; isig<NSig; isig++)
	fPSC[detid][segid]->at(ipsc)->spulse[iseg][isig] =
	  ftmp*(fPSC[detid][segid]->at(ipsc)->spulse[iseg][isig]*nhitstmp -
		aps->opulse[iseg][isig]);
      
    fPSC[detid][segid]->at(ipsc)->nhits--;
    if(fPSC[detid][segid]->at(ipsc)->nhits==0){ cPSCtotal--; cPSCmem--;}

    // remove from HitCollection
    fHCs[detid][segid]->at(ipsc)->RemoveHit(ahit);
    fHCs[detid][segid]->at(ipsc)->SetRealPosition(fPSC[detid][segid]->at(ipsc)->labpos);

    if(fHCs[detid][segid]->at(ipsc)->GetSize()==0){
      freeHCs[detid][segid].push_back(ipsc);
      cHCs--;
    }

    if(jpsc>-1){
      ahit->RemoveHitCollection(ahc);
    }
  }//end of loop hcs

  if(jpsc<0){
    ahit->ClearHitCollection();
    cPStotal--;
    cHits--;

  }else{
    if(ahit->hasHitCollection()<1){
      cPStotal--;
      cHits--;
    }
  }

  return;
}


void AGATA::RemovePSC(HitCollection *ahc){ // empty ahc from fHCs
  int detid = ahc->GetDet();
  int segid = ahc->GetSeg();
  int ipsc  = ahc->GetPid();

  vector<Hit*>* fhits = ahc->GetHits();
  if(fhits->size()==0) return;
  
  for(Hit* ah : *fhits){
    ah->RemoveHitCollection(ahc);
    if(ah->hasHitCollection()==0){
      cPStotal--;
      cHits--;
    }
  }
  ahc->Clear();
  cPSCtotal--; cPSCmem--;
   
  fPSC[detid][segid]->at(ipsc)->nhits = 0;

  for(int ix=0; ix<3; ix++){
    fPSC[detid][segid]->at(ipsc)->labpos[ix] = 0;
    fPSC[detid][segid]->at(ipsc)->detpos[ix] = 0;
  }

  for(int is=0; is<NSeg_comp; is++){
    fPSC[detid][segid]->at(ipsc)->divzone[is] = -1;
    fPSC[detid][segid]->at(ipsc)->devabscut[is][0] = -1;
    fPSC[detid][segid]->at(ipsc)->devsigma[is] = -1;
    fPSC[detid][segid]->at(ipsc)->devabs[is].clear();
    fPSC[detid][segid]->at(ipsc)->dev[is].clear();
  }
  fPSC[detid][segid]->at(ipsc)->dividx.clear();
  
  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++)
      fPSC[detid][segid]->at(ipsc)->spulse[iseg][isig] = 0;
  }

  freeHCs[detid][segid].push_back(ipsc);
  cHCs--;
  
  return;
}


void AGATA::Devseg(const float *apulse, const float *bpulse, float *dev){

  //float apulsediff[NSig_comp], bpulsediff[NSig_comp];
  //std::adjacent_difference( std::begin(apulse), std::end(apulse), std::begin(apulsediff));
  //std::adjacent_difference( std::begin(bpulse), std::end(bpulse), std::begin(bpulsediff));

  float apulsesum[NSig_comp], bpulsesum[NSig_comp];
  std::partial_sum( apulse, apulse+NSig_comp, apulsesum);
  std::partial_sum( bpulse, bpulse+NSig_comp, bpulsesum);

  float asum, bsum;
  asum = std::accumulate( apulsesum, apulsesum+NSig_comp, 0.);
  bsum = std::accumulate( bpulsesum, bpulsesum+NSig_comp, 0.);

  dev[2] = asum - bsum;

#ifdef SSE_M256

  const __m256 masks = _mm256_set1_ps(-0.0f);
  
  __m256* realtrace = (__m256*)apulse;
  __m256* basetrace = (__m256*)bpulse;

  __m256 diff = _mm256_setzero_ps();
  __m256 devori = _mm256_setzero_ps();
  __m256 devabs = _mm256_setzero_ps();
  /*
  __m256* realdiff = (__m256*)apulsediff;
  __m256* basediff = (__m256*)bpulsediff;

  __m256 diffdiff = _mm256_setzero_ps();
  __m256 diffdevori = _mm256_setzero_ps();
  __m256 diffdevabs = _mm256_setzero_ps();
  */

  for(int nn=0; nn<LOOP_SSE8_seg; nn++){
    diff = _mm256_sub_ps(realtrace[nn], basetrace[nn]);
    devori = _mm256_add_ps(devori, diff); // original value
    devabs = _mm256_add_ps(devabs, _mm256_andnot_ps(masks, diff)); // ABS value
    /*
    diffdiff = _mm256_sub_ps(realdiff[nn], basediff[nn]);
    diffdevori = _mm256_add_ps(diffdevori, diffdiff); // original value
    diffdevabs = _mm256_add_ps(diffdevabs, _mm256_andnot_ps(masks, diffdiff)); // ABS value
    */
  }

  __m256 tmp0 = _mm256_permute2f128_ps(devabs, devabs, 1);
  devabs = _mm256_add_ps(devabs, tmp0);
  devabs = _mm256_hadd_ps(devabs, devabs);
  devabs = _mm256_hadd_ps(devabs, devabs);

  __m256 tmp1 = _mm256_permute2f128_ps(devori, devori, 1);
  devori = _mm256_add_ps(devori, tmp1);
  devori = _mm256_hadd_ps(devori, devori);
  devori = _mm256_hadd_ps(devori, devori);

  dev[0] = _mm256_cvtss_f32(devabs);
  dev[1] = _mm256_cvtss_f32(devori);
  /*
  __m256 tmp2 = _mm256_permute2f128_ps(diffdevabs, diffdevabs, 1);
  diffdevabs = _mm256_add_ps(diffdevabs, tmp2);
  diffdevabs = _mm256_hadd_ps(diffdevabs, diffdevabs);
  diffdevabs = _mm256_hadd_ps(diffdevabs, diffdevabs);

  __m256 tmp3 = _mm256_permute2f128_ps(diffdevori, diffdevori, 1);
  diffdevori = _mm256_add_ps(diffdevori, tmp3);
  diffdevori = _mm256_hadd_ps(diffdevori, diffdevori);
  diffdevori = _mm256_hadd_ps(diffdevori, diffdevori);

  dev[2] = _mm256_cvtss_f32(devabs3);
  dev[3] = _mm256_cvtss_f32(devori3);
  */

#else // use m128

  const __m128 masks = _mm_set1_ps(-0.0f);
  const __m128 zeros = _mm_setzero_ps();

  __m128* realtrace = (__m128*)apulse;
  __m128* basetrace = (__m128*)bpulse;

  __m128 diff = _mm_setzero_ps();
  __m128 devori = _mm_setzero_ps();
  __m128 devabs = _mm_setzero_ps();
  /*
  __m128* realdiff = (__m128*)apulsediff;
  __m128* basediff = (__m128*)bpulsediff;

  __m128 diffdiff = _mm_setzero_ps();
  __m128 diffdevori = _mm_setzero_ps();
  __m128 diffdevabs = _mm_setzero_ps();
  */

  for(int nn=0; nn<LOOP_SSE4_seg; nn++){
    diff = _mm_sub_ps(realtrace[nn], basetrace[nn]);
    devori = _mm_add_ps(devori, diff); // original value
    devabs = _mm_add_ps(devabs, _mm_andnot_ps(masks, diff)); // ABS value
    /*
    diffdiff = _mm_sub_ps(realdiff[nn], basediff[nn]);
    diffdevori = _mm_add_ps(diffdevori, diffdiff); // original value
    diffdevabs = _mm_add_ps(diffdevabs, _mm_andnot_ps(masks, diffdiff)); // ABS value
    */
  }

  devabs = _mm_hadd_ps(_mm_hadd_ps(devabs, zeros), zeros);
  devori = _mm_hadd_ps(_mm_hadd_ps(devori, zeros), zeros);

  dev[0] = _mm_cvtss_f32(devabs);
  dev[1] = _mm_cvtss_f32(devori);
  /*
  diffdevabs = _mm_hadd_ps(_mm_hadd_ps(diffdevabs, zeros), zeros);
  diffdevori = _mm_hadd_ps(_mm_hadd_ps(diffdevori, zeros), zeros);

  dev[2] = _mm_cvtss_f32(diffdevabs);
  dev[3] = _mm_cvtss_f32(diffdevori);
  */
#endif

  return;
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


void AGATA::CheckPSCstat(long long *PSCstat){
  maxnhits = 0;
  int nPSC = 0;
  int nEmpty = 0;

  for(HitCollection *ahc : *fAllHCs){
    int tmpnhits = ahc->GetSize();
    if(tmpnhits>maxnhits) maxnhits = tmpnhits;
    if(tmpnhits>0)  nPSC++;
    if(tmpnhits==0) nEmpty++;
  }

  PSCstat[0] = maxnhits;
  PSCstat[1] = nPSC;
  PSCstat[2] = nEmpty;

  return;
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


long long AGATA::AddEventHits(EventHits* fEvent){
#ifdef NTHREADS
  lock_guard<mutex> lock(EvtHitsmtx);
#endif
  fEventHits->push_back(fEvent);
  return fEventHits->size()-1;
}


long long AGATA::FindiEvtHit(int iconfig, int irun, int ientry, long long &istart){
  long long iEvtHit = -1;
  int tmpconf, tmprun, tmpetry;
  for(long long i=istart; i<fEventHits->size(); i++){
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
  long long Nevts = fEventHits->size();
  time(&start);

  for(; ievt<Nevts;){ // loop events

    // get Hits for ievt
    EventHits *fEvent;
    vector<Hit*>* fHits;
    vector<TVector3> sourcepos; // mm
    vector<float> EGamma; // keV
#ifdef DIFFTOTE
    float Etot;
#endif
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
      fEvent = fEventHits->at(ievt);
      fHits = fEvent->GetfHits();
#ifdef DIFFTOTE
      Etot = fEvent->Etot;
#endif
      sourcepos.clear();  EGamma.clear();
      int nsource = fEvent->GetNSource();
      for(int is=0; is<nsource; is++){
	sourcepos.push_back(fEvent->GetSourcePos(is));
	EGamma.push_back(fEvent->GetSourceE(is));
      }
      ievt++;
    }

    // calc calpos from HCs
    for(Hit* ahit : *fHits) ahit->CalcAveHCsPosition(0);

    // tracking
    vector<int> sign(fHits->size(),-1);
    int Nunsigned = sign.size();
    int iclust = 0;

    // analysis clusters
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

      int nsource = EGamma.size();
      vector<int> atrack;
      int bestis = 0;
      double minchi2 = 1e9;
      for(int is=0; is<nsource; is++){
	Tracker tracker(tHits, EGamma[is], sourcepos[is]);

#ifdef DIFFTOTE
#if    DIFFTOTE < 0 // DIFFTOTE<0 : only accept TotE match events
	if( !(fabs(Etot-EGamma[is])<fabs(DIFFTOTE)) ) continue;
#else   // DIFFTOTE>0 : if TotE match, put all hits in one clust
	if( fabs(Etot-EGamma[is])<fabs(DIFFTOTE) ) tracker.SetOneClust(true);
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
	fEvent->SignClust( iclust, hid);
	Nunsigned--;
      }


#ifdef TRACKINGTREE
      if(atrack.size()>1){
	Tracker tracker(tHits, EGamma[bestis], sourcepos[bestis]);
#ifdef ONECLUST
	tracker.SetOneClust(true);
#endif
	tracker.OFTtracking();
	//tracker.Simpletracking();

	lock_guard<mutex> lock(Trtreemtx);
	Trnhits = atrack.size();
	TrSource = true;
	TrSourceE = EGamma[bestis];
	TrSourcePos[0] = sourcepos[bestis].X();  TrSourcePos[1] = sourcepos[bestis].Y();  TrSourcePos[2] = sourcepos[bestis].Z();
	TrCorrect = tracker.CheckOrder();
	TrFOM1 = tracker.GetFOM1();
	TrFOM2 = tracker.GetFOM2();
	Trtree->Fill();
      }
#endif
    
      // make paths
      int ngoodhit = 0;
      double incE = EGamma[bestis];
      double depE = tHits->at(atrack[0])->GetE(); // keV

      Hit *sourcehit = new Hit(sourcepos[bestis]);
      sourcehit->SetInterid(-1);
#ifdef NTHREADS2
      lock_guard<mutex> lock(Pathsmtx);
#endif
      ngoodhit = 1; // sourcehit is always good hit
      //if(tHits->at(atrack[0])->hasHitCollection()>0) ngoodhit++;
      //if(tHits->at(atrack[1])->hasHitCollection()>0) ngoodhit++;
      if(tHits->at(atrack[0])->hasgoodHCs(0)>0) ngoodhit++;
      if(tHits->at(atrack[1])->hasgoodHCs(0)>0) ngoodhit++;
      cPathsN[ngoodhit]++;

      //if(ngoodhit>1){ // at least two good hit
      if(ngoodhit==3){ // three good hit
	Path *apath = new Path(sourcehit,tHits->at(atrack[0]),tHits->at(atrack[1]),
			       incE, depE, incE, depE);
	fPaths->push_back(apath);
	cPaths++;
      }
      incE = incE - depE;

      for(int i=1; i<atrack.size()-1; i++){ //<--- was atrack.size()-2 ????
	depE = tHits->at(atrack[i])->GetE(); // keV

	ngoodhit = 0;
	//if(tHits->at(atrack[i-1])->hasHitCollection()>0) ngoodhit++;
	//if(tHits->at(atrack[i]  )->hasHitCollection()>0) ngoodhit++;
	//if(tHits->at(atrack[i+1])->hasHitCollection()>0) ngoodhit++;
	if(tHits->at(atrack[i-1])->hasgoodHCs(0)>0) ngoodhit++;
	if(tHits->at(atrack[i]  )->hasgoodHCs(0)>0) ngoodhit++;
	if(tHits->at(atrack[i+1])->hasgoodHCs(0)>0) ngoodhit++;
	cPathsN[ngoodhit]++;
	
	//if(ngoodhit>1){ // at least two good hit
	if(ngoodhit==3){ // three good hit
	  Path *apath = new Path(tHits->at(atrack[i-1]),tHits->at(atrack[i]),tHits->at(atrack[i+1]),
				 incE, depE, incE, depE);
	  fPaths->push_back(apath);
	  cPaths++;
	}
	incE = incE - depE;
      }

      iclust++;
      delete tHits;
    }// end of loop clusters in one event 

  }// end of loop events

  return;
}


void AGATA::Tracking(int iter){
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

#ifdef TRACKINGTREE
  Trfile = new TFile(Form("share/Trtree%d.root",iter),"RECREATE");
  Trtree = new TTree("tree","tracking results");
  Trtree->Branch("nhits",&Trnhits);
  Trtree->Branch("sourceE",&TrSourceE);
  Trtree->Branch("source",&TrSource);
  Trtree->Branch("sourcepos",TrSourcePos,"sourcepos[3]/F");
  Trtree->Branch("correct",&TrCorrect);
  Trtree->Branch("FOM1",&TrFOM1);
  Trtree->Branch("FOM2",&TrFOM2);
#endif
  
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

  cout<<" .. Paths0-"<<cPathsN[0]<<" .. Paths1-"<<cPathsN[1]
      <<" .. Paths2-"<<cPathsN[2]<<" .. Paths3-"<<cPathsN[3]<<" .."<<endl;
  
#ifdef TRACKINGTREE
  Trfile->cd();
  Trtree->Write();
  Trfile->Close();
#endif
  
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
    //phits->at(i)->CalcAveHCsPosition(0); // update hits position in paths
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
    //if(kPSA && aHC->GetPaths()->size()<50) continue; // if PSA used for initial hit pos, then only optimize for HC with npaths>50
    
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
    for(Hit* ah : *phits) ah->CalcAveHCsPosition(0); // update hits pos linked with aHC
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
    Float_t dbposi[3];
    Float_t dbspulsei[NSig*NSegCore];

    dbtree->SetBranchAddress("seg",&dbsegi);
    dbtree->SetBranchAddress("pos",dbposi);
    dbtree->SetBranchAddress("spulse",dbspulsei);
    int npoint = dbtree->GetEntriesFast();

    int ncounter = 0;
    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);

      //dbsegi = dbsegi-1; // start from 0

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
	  else          tmpamp = dbspulsei[NSeg*NSig+isig];

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
  if(kPSA){

    int fseg[NSeg_comp];
    agatageo->GetNextSegs(segid, fseg);
    int npoint = fPSAbasis[itype][segid].size();
    float minchi2 = 1e9;
    float asegpulse[NSig_comp], bsegpulse[NSig_comp];
    // compare fired seg, core and neighbor segment
    for(int ipoint=0; ipoint<npoint; ipoint++){

      float chi2=0;
      for(int is=0; is<NSeg_comp; is++){
	int iseg = fseg[is];
        copy_n( aps->opulse[iseg], NSig_comp, asegpulse);
        copy_n(fPSAbasis[itype][segid][ipoint].spulse[is], NSig_comp, bsegpulse);
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

  TVector3 pos;
  if(ipos>-1){
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
  tree->Branch("calpos",calpos,"calpos[3]/F");
  tree->Branch("cadpos",cadpos,"cadpos[3]/F");
  tree->Branch("calpos2",calpos2,"calpos2[3]/F");
  tree->Branch("cadpos2",cadpos2,"cadpos2[3]/F");

  tree->Branch("labpos",labpos,"labpos[3]/F");
  tree->Branch("detpos",detpos,"detpos[3]/F");
  tree->Branch("dist",&dist);
  tree->Branch("dist2",&dist2);

  tree->Branch("spulse",spulse,Form("spulse[%d][%d]/F",NSegCore,NSig));  
  tree->Branch("devsigma",devsigma,Form("devsigma[%d]/F",NSeg_comp));  

  tree->Branch("npaths",&npaths);
}

void AGATA::InitTreeRead(TTree *tree){
  tree->SetBranchAddress("det",&det);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("index",&index);
  tree->SetBranchAddress("nhits",&nhits);
  tree->SetBranchAddress("calpos",calpos);
  tree->SetBranchAddress("cadpos",cadpos);
  tree->SetBranchAddress("calpos2",calpos2);
  tree->SetBranchAddress("cadpos2",cadpos2);

  tree->SetBranchAddress("labpos",labpos);
  tree->SetBranchAddress("detpos",detpos);

  tree->SetBranchAddress("spulse",spulse);
  tree->SetBranchAddress("devsigma",devsigma);

  tree->SetBranchAddress("npaths",&npaths);
}

void AGATA::ClosePSCFiles(){
  for(int idet=0; idet<NDets; idet++){
    if(SkipDet[idet]) continue;
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


