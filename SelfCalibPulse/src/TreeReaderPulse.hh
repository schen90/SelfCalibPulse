#ifndef TREEREADERPULSE_HH
#define TREEREADERPULSE_HH

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
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
  TreeReaderPulse();
  virtual ~TreeReaderPulse();

  virtual void Load(string configfile);
  virtual void MakeInit();
  virtual void MakeNoise();
  virtual void LoadNoise();

  virtual void Init(int iChain);
  virtual void InitEvtTree(int iChain);

  virtual void GeneratePSC(AGATA *agata);
  virtual void GeneratePSC(AGATA *agata, int nevts);
  virtual void GeneratePSC(AGATA *agata, int iconfig, int nevts);
  
  virtual void GeneratePSCLoop(int iChain, AGATA *agata,
			       int iconfig, int nentries);

  virtual void ProcessOneEvent(int iChain, AGATA *agata,
			       int iconfig, int run,
			       int ientry, int nentries);

  void SetMaxMemUsage(double value){ MaxMemUsage = value;}
  Double_t GetTotalSystemMemory();
  Double_t GetCurrentMemoryUsage();


private:
  Double_t MaxMemUsage = 100; // max memory usage %

  TChain* fChain[NChain]; //!pointer to the analyzed tree

#ifdef NOISE
  float noise[NOISE]; // base for noise
#endif

  int nConfig = 0;
  vector<float> fSourceE;
  vector<TVector3> fSourcePos;
  vector<string> path;
  vector<int> MinRun;
  vector<int> MaxRun;
  vector<int> Nevts;

  float SourceE;
  TVector3 SourcePos;

  struct OBJ{
    // Declaration of leaf types
    int                      ievent;     // event id
    vector<int>             *ndet = 0;   // detector id
    vector<int>             *g4seg = 0;  // segment id
    vector<float>           *energy = 0; // deposit energy
    
    vector<vector<float>>   *posa = 0;   // absolute/global interaction position vector<float(3)>
    vector<vector<float>>   *posr = 0;   // relative/local interaction position vector<float(3)>
    
    // pulse shape vector<>
    vector<int>             *pdet = 0;   // detector id for pulse shape
    vector<float>           *ecore = 0;  // core energy
    vector<vector<int>>     *inter = 0;  // interaction id in G4 vector<>
    
    vector<vector<int>>     *pseg = 0;   // segment id in pulsedb
    vector<vector<int>>     *ngrid = 0;  // number of grid found around ppos
    vector<vector<int>>     *extrpl = 0;  // number of grid for extrapolation
    vector<vector<float>>   *core = 0;   // core pulse shape vector<float(56)>
    vector<vector<float>>   *spulse = 0; // segment pulse shape vector<float(2016)>
    
    int                      category; // 1: max 1 seg fired in a det, >1 det fired; 2: max >1 seg fired in a det
  };
  OBJ obj[NChain];


  TFile *EvtFile[NChain];  
  TTree *EvtTree[NChain];
  
  struct EvtOBJ{
    int      iconfig;
    int      irun;
    int      ientry;
    int      nhits;

    int      interid;
    int      idet;
    int      iseg;
    float    depE;
    float    labpos[3];
    float    calpos[3];
    float    dist;
#ifdef NOISE
    int      noiseidx;
    int      noiseidxshift;
#endif
    float    sourcepos[3];
    float    sourceeng;
  };
  EvtOBJ evtobj[NChain];

  
#ifdef ADDPS
  PSbasis *apsb;
#endif

  PS GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift);
  PS GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift, bool skipPS);
  PS GetAPS(int iChain, AGATA *agata, int idet, int nidx, int nidxshift, bool skipPS, int &segidx);

  atomic_int irun;
  atomic_int ievt;
  atomic<bool> kcout;
  mutex treemtx; // tree lock for threads read
  time_t start, stop;

  atomic_int cPaths;
  int iter;
};


#endif // #ifndef TREEREADERPULSE_HH
