#ifndef TREEREADERPULSE_H
#define TREEREADERPULSE_H

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TVector3.h>
#include <stdlib.h>
#include <stdio.h>
#include <atomic>
#include <chrono>
#include <thread>

#include "Global.hh"
#include "AGATA.hh"
#include "Hit.hh"
#include "EventHits.hh"
#include <vector>
#include <unistd.h>


#ifdef NTHREADS
#define NChain NTHREADS
#else
#define NChain 1
#endif

using namespace std;

class TreeReaderPulse {

public:
  TreeReaderPulse(int detid);
  virtual ~TreeReaderPulse();

  virtual void Load(string configfile);
  virtual void MakeInit();
  
  virtual void Init(int iChain);

  virtual void GenerateHCs(int opt, AGATA *agata);
  virtual void GenerateHCs(int opt, AGATA *agata, long long nevts);
  virtual void GenerateHCs(int opt, AGATA *agata, long long nevts, int iconfig);

  virtual void GenerateHCsLoop(int opt, int iconfig, int iChain, AGATA *agata, long long nentries);

  virtual int GenerateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
				int ientry, long long nentries);

  virtual int FindDevworker(int opt, int iconfig, int run, int iChain, AGATA *agata,
			    int ientry, long long nentries, long long &istart);

  virtual int UpdateHCsworker(int opt, int iconfig, int run, int iChain, AGATA *agata,
			      int ientry, long long nentries, long long &istart);

  Double_t GetTotalSystemMemory();
  Double_t GetCurrentMemoryUsage();
  void SetMaxMemUsage(double value){ MaxMemUsage = value;}

  int GetRemovePSNumber(){ return cRemovePS;}

  void ClearSkipDetId(){ for(bool &val : SkipDet){ val=false;}}
  void SkipDetId(int val){ SkipDet[val]=true;}
  
private:
  TChain* fChain[NChain]; //!pointer to the analyzed tree

  int Detid = -1;
  bool SkipDet[MaxNDets];

  long long NEventHits; // size of fEventHits in AGATA.hh

  atomic_int cDivPS;
  atomic_int cRemovePS;
  atomic_int cNotMatch;

  int nConfig = 0;
  vector<Config> fConfigs;

  Double_t MaxMemUsage = 50; // max memory usage %

  vector<float> SourceE;
  vector<TVector3> SourcePos;
  
  struct OBJ{
    Int_t     EntryID;
    Float_t   SegTraces[DSIZE*NSEGS];
    Float_t   CoreTraces[DSIZE*2];
    Float_t   SegE[NSEGS];
    Float_t   CoreE[2];
    Float_t   CoreT[2];
    Int_t     CrystalId;
    ULong64_t CrystalTS;
  };

  OBJ obj[NChain];


  PS GetAPS(int iChain, bool skipPS, int &segidx);
  atomic_int irun;
  atomic<long long> ievt;
  atomic<bool> kcout;
  atomic<float> MaxDev;
  atomic_int maxnhitsdiv;
  int kInterrupt;
  mutex treemtx; // tree lock for threads read
  time_t start, stop;

  atomic<long long> cnevents; // counter for total events from all input files
  atomic<long long> cievt; // counter for total analyzied events 
  atomic_int cievthitfind; // counter for ievthit find
  atomic_int cievthitnotfind; // counter for ievthit not find

};


#endif // #ifndef TREEREADERPULSE_H
