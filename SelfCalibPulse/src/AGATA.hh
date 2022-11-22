#ifndef AGATA_HH
#define AGATA_HH

#include <TVector3.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TTree.h>
#include <time.h>
#include <vector>
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMinuit.h"
#include "TMinuitMinimizer.h"

#include "Global.hh"
#include "AGATAgeo.hh"
#include "Hit.hh"
#include "EventHits.hh"
#include "HitCollection.hh"
#include "PSC.hh"
#include "Path.hh"
#include "Tracker.hh"

#ifdef MINUIT2
#include "Minuit2/Minuit2Minimizer.h"
#endif

using namespace std;
  
class AGATA {

public:
  AGATA(int detid);
  virtual ~AGATA();

  AGATAgeo* GetGeo(){ return agatageo;}

  // write/read files
  void WritePSCfiles(); // create Pulse Shape Collection files
  void WritePSCfiles(int detid); // create Pulse Shape Collection files, detid=-1 all dets
  void LoadPSCfiles(); // load all Pulse Shape Collection in memory
  void LoadPSCfiles(int detid); // load Pulse Shape Collection in memory, detid=-1 all dets
  void ClearPSCMem(); // clear PSC in memory
  void GetPSCstat(long long *PSCstat); // get Pulse Shape Collection state

  void WriteHCfiles();
  void WriteHCfiles(int detid);
  void LoadHCfiles();
  void LoadHCfiles(int detid);
  void ClearHCMem(); // clear HC in memory

  void WriteEvtHitsfiles(int detid);
  void Load(string configfile);
  void CombEvtHitsfiles();
  void LoadEvtHitsfiles2(int iconfig);
  void LoadEvtHitsconfigs();
  void ClearEvtHitsMem(); // clear EventHits in memory

  // Pulse Shape Collection
  void SetAddNewPSC(bool val){ 
    if(kAddNewPSC!=val){
      cout<<endl<<"SetAddNewPSC "<<val<<endl;
      kAddNewPSC = val;
    }
  }
  void SortPSC(); // sort PSC increasing ientry
  int  FindHC(int detid, int segid, int pscid);

  int   InitPSCandHC(int detid, int segid);

  void  FindInitZone(PS *aps, vector<int> *initzone);
  void  FindInitPSC(Hit *ahit, vector<int> *initpsc);

  int   AddPS(PS *aps, Hit *ahit);
  int   AddPStoPSC(PS *aps, Hit *ahit, int ipsc);

  float FindMaxDev(PS *aps, Hit *ahit);
  void  FindDevCut();

  void  FindDevSigma(PS *aps, Hit *ahit);
  void  CalcDevSigma();

  void  FindDivZone(PS *aps, PSC *apsc, vector<vector<int>> *divzone);
  int   AddPStoDiv(PS *aps, Hit *ahit);

  int   CheckPSinPSC(PS *aps, Hit *ahit);
  void  SetNSigma(int val){ nSigma = val; cout<<Form("nSigma = %.1f",nSigma)<<endl;}
  float GetNSigma(){ return nSigma;}

  void MakeCPulse();
  
  void  RemoveMotherPSC();
  void  RemoveSmallPSC(int minhits);

  void RemovePSfromPSC(PS *aps, Hit *ahit, int ipsc); // remove ahit from PSC(s)
  void RemovePSC(HitCollection *ahc);

  void Devseg(const float *apulse, const float *bpulse, float *dev); // calc deviation of two segment pulse shape
  Float_t Chi2seg(const float *spulse, const float *PSCspulse); // compare two segment pulse shape

  void CheckPSCstat(long long *PSCstat);

  // EventHits
  void SortEventHits(); // sort EventHits increasing iconfig, irun, ientry
  long long AddEventHits(EventHits* fHits);
  long long FindiEvtHit(int iconfig, int irun, int ientry, long long &istart);
  EventHits* FindEventHits(long long i){ return fEventHits->at(i);}
  long long GetEventHitsSize(){ return fEventHits->size();}

  // Tracking
  void TrackingLoop(); // tracking fHits
  void Tracking(int iter); // loop all events
  
  void RegisterPathswithHCs();

  // HCs pos optimize
  //Double_t WrappedEstimator(const double *par);
  void ExecFitLoop(int it);
  void ExecFit(int repeat);


  // PSA to assign initial pos
  void ReadPSAbasis();
  TVector3 GetPSpos(int detid, int segid, PS *aps);
  

  // Check Memory
  void SetMaxMemUsage(double value){ MaxMemUsage = value;}
  Double_t GetTotalSystemMemory();
  Double_t GetCurrentMemoryUsage();

  // run option
  void SetPSA(bool val){ kPSA = val;}

  void SetMaxNDiv(int val){ maxndiv = val;}
  
private:
  int Detid = -1;
  AGATAgeo *agatageo;

  Double_t MaxMemUsage = 50; // max memory usage %
  bool kPSA = false;
  
  // Pulse Shape Collection and Hit Collection
  TFile *pscfile[MaxNDets];
  TTree *psctree[MaxNDets][NSeg];

  int PSClimit[MaxNDets]; //PSC number limit
  vector<PSC*>* fPSC[MaxNDets][NSeg]; // Pulse Shape Collection storage
  vector<HitCollection*>* fHCs[MaxNDets][NSeg]; // Hit Collection storage
  vector<Int_t> freeHCs[MaxNDets][NSeg]; // idx of empty fHCs
  vector<HitCollection*>* fAllHCs; // Hit Collection storage
  atomic_int ihc;
  mutex PSCmtx[MaxNDets][NSeg]; // fPSC lock for threads
  mutex AllHCmtx;
  bool kAddNewPSC = true;
  float nSigma = 3.;

  vector<Int_t> HCMap[MaxNDets][NSeg];
  Int_t         HCstat[MaxNDets][NSeg][2]; // 0: fHCs size, 1: fHCs max idx
  
  // EventHits
  int nConfig = 0;
  vector<int> MinRun;
  vector<int> MaxRun;
  vector<EventHits*>* fEventHits; //fHits in Events for tracking
  atomic_int atomrun;
  atomic<long long> ievt;
  atomic<long long> NevtsTotal;
  mutex EvtHitsmtx;

#ifdef TRACKINGTREE
  TFile *Trfile;
  TTree *Trtree;
  int  Trnhits;
  bool TrSource;  // if start from source
  float TrSourceE;
  float TrSourcePos[3];
  bool TrCorrect;
  double TrFOM1;
  double TrFOM2;
  mutex Trtreemtx;
#endif
  
  // Path
  vector<Path*>* fPaths; // Path storage
  mutex Pathsmtx;

  // parameters for HC pos optimize
  Double_t fitlimit; // boundary for the fit parameters
  
  // PSA to assign initial pos
  vector<PSAbasis> fPSAbasis[NType][NSeg];
  
  // branch
  int   det;           // detector id
  int   seg;           // segment id
  int   index;         // index for PSCid
  int   nhits;         // number of hits in group

  float calpos[3];     // init selfcalib interaction position in labframe float(3)
  float cadpos[3];     // init selfcalib interaction position in detframe float(3)
  float calpos2[3];    // selfcalib interaction position in labframe float(3)
  float cadpos2[3];    // selfcalib interaction position in detframe float(3)

  float labpos[3];     // average interaction position in labframe float(3)
  float detpos[3];     // average interaction position in detframe float(3)  
  float dist;          // dist calpos - labpos
  float dist2;         // dist calpos2 - labpos

  float spulse[NSegCore][NSig];  // average pulse shape
  float devsigma[NSeg_comp];     // standard deviation of compared segment
  
  int npaths;

  // set tree branch
  void InitTreeWrite(TTree *tree);
  void InitTreeRead(TTree *tree);
  void ClosePSCFiles();
  
  Int_t NDets; // number of detectors from LookUpTable

  // counter
  atomic<long long> cPStotal;
  atomic<long long> cPSCtotal;
  atomic<long long> cPSCmem;
  atomic<long long> cPSCfile;
  atomic<long long> maxnhits;
  atomic<long long> cPaths;
  atomic<long long> cHits;
  atomic<long long> cHCs;
  atomic<long long> cPathsN[4];
  atomic<long long> maxndiv;

  time_t start, stop;

};

#endif // #ifndef AGATA_HH
