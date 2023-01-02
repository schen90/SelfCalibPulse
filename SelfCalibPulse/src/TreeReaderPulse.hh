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
  
private:
  TChain* fChain[NChain]; //!pointer to the analyzed tree

  int Detid = -1;
  long long NEventHits; // size of fEventHits in AGATA.hh

  atomic_int cDivPS;
  atomic_int cRemovePS;
  atomic_int cNotMatch;

  int nConfig = 0;
  vector<int>              NSource;
  vector<vector<float>>    fSourceE;
  vector<vector<TVector3>> fSourcePos;
  vector<string> path;
  vector<int> MinRun;
  vector<int> MaxRun;
  vector<long long> Nevts;

  Double_t MaxMemUsage = 50; // max memory usage %

  vector<float> SourceE;
  vector<TVector3> SourcePos;
  
  struct OBJ{
    Int_t     EntryID;
    Float_t   SegTraces[2160];
    Float_t   CoreTraces[120];
    Float_t   SegE[36];
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

  // ScanPS
  vector<PS> fPSs[3];
  mutex scanmtx;

};


#endif // #ifndef TREEREADERPULSE_H
