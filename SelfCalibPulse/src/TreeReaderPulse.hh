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

#ifdef ADDPS
#include "PSbasis.hh"
#endif

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
  virtual void MakeNoise();
  virtual void LoadNoise();
  
  virtual void Init(int iChain);

  virtual void GenerateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
				 int ientry, int nentries);
  virtual void GenerateHCsLoop(int iconfig, int iChain, AGATA *agata, int nentries);
  virtual void GenerateHCs(AGATA *agata);
  virtual void GenerateHCs(AGATA *agata, int nevts);
  virtual void GenerateHCs(AGATA *agata, int nevts, int iconfig);

  void SetNewPSC(bool val){ kNewPSC = val;}
  void SetUpdateHCs(int val){ kUpdateHCs = val;}
  virtual void UpdateHCsworker(int iconfig, int run, int iChain, AGATA *agata,
			       int ientry, int nentries, int &istart);

  Double_t GetTotalSystemMemory();
  Double_t GetCurrentMemoryUsage();
  void SetMaxMemUsage(double value){ MaxMemUsage = value;}

  void SetGroupPos(bool val){ kGroupPos = val;}
  
private:
  TChain* fChain[NChain]; //!pointer to the analyzed tree

  int Detid = -1;
  int NEventHits; // size of fEventHits in AGATA.hh
#ifdef NOISE
  float noise[NOISE]; // base for noise
#endif

  bool kNewPSC;
  int kUpdateHCs = 0;
  atomic_int cRemoveHit;
  atomic_int cAddHit;
  atomic_int cNotMatch;

  int nConfig = 0;
  vector<float> fSourceE;
  vector<TVector3> fSourcePos;
  vector<string> path;
  vector<int> MinRun;
  vector<int> MaxRun;
  vector<int> Nevts;

  bool kWithPS = true;
  bool kGroupPos = false;
  
  Double_t MaxMemUsage = 50; // max memory usage %

  float SourceE;
  TVector3 SourcePos;
  
  struct OBJ{
    // Declaration of leaf types
    int                      ievent;     // event id
    vector<int>             *ndet = 0;   // detector id
    vector<int>             *g4seg = 0;  // segment id
    vector<float>           *energy = 0; // deposit energy
#ifdef REALPOS
    vector<vector<float>>   *posa = 0;   // absolute/global interaction position vector<float(3)>
    vector<vector<float>>   *posr = 0;   // relative/local interaction position vector<float(3)>
#endif
    // pulse shape vector<>
    vector<int>             *pdet = 0;   // detector id for pulse shape
    vector<float>           *ecore = 0;  // core energy
    vector<vector<int>>     *inter = 0;  // interaction id in G4 vector<>
  
    vector<vector<int>>     *pseg = 0;   // segment id in pulsedb
    vector<vector<int>>     *ngrid = 0;  // number of grid found around ppos
    vector<vector<int>>     *extrpl = 0;  // number of grid for extrapolation
    vector<vector<float>>   *core = 0;   // core pulse shape vector<float(121)>
    vector<vector<float>>   *spulse = 0; // segment pulse shape vector<float(4356)>

    int                      category; // 1: max 1 seg fired in a det, >1 det fired; 2: max >1 seg fired in a det
  };

  OBJ obj[NChain];

#ifdef ADDPS
  PSbasis *apsb;
#endif
  
  PS GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift);
  PS GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift, bool skipPS);
  atomic_int irun;
  atomic_int ievt;
  atomic<bool> kcout;
  double minchi2;
  int kInterrupt;
  mutex treemtx; // tree lock for threads read
  time_t start, stop;

};


#endif // #ifndef TREEREADERPULSE_H
