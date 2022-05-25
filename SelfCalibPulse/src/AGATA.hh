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
  void GetPSCstat(int *PSCstat); // get Pulse Shape Collection state

  void WriteHCfiles();
  void WriteHCfiles(int detid);
  void LoadHCfiles();
  void LoadHCfiles(int detid);
  void ClearHCMem(); // clear HC in memory

  void WriteEvtHitsfiles(int detid);
  void Load(string configfile);
  void CombEvtHitsfiles();
  void LoadEvtHitsfiles(int iconfig);
  void LoadEvtHitsfiles2(int iconfig);
  void LoadEvtHitsconfigs();
  void ClearEvtHitsMem(); // clear EventHits in memory

  // Pulse Shape Collection
  void SetAddNewPSC(bool val);
  bool IsNewPSC(); // if new PSC added

  void SortPSC(); // sort PSC increasing ientry
  Double_t AddPStoPSC(PS *aps, Hit *ahit, vector<int> &entrylist);
  void RemovePSfromPSC(PS *aps, Hit *ahit); // remove ahit from all PSCs
  void RemovePSC(HitCollection *ahc);
  void DividePSC(HitCollection *ahc, double factor);
  void CheckPSCs(int minpaths, int maxpaths);
  void ClearHitLevelMarker(int val);

  void SetPSClimit(int detid, int val){ PSClimit[detid] = val;}
  void SetMaxChi2(float value){ MaxChi2 = value;}
  void SetMaxChi2s(float val0, float val1, float val2){ MaxChi2s[0] = val0; MaxChi2s[1] = val1; MaxChi2s[2] = val2;}
  Float_t Dist(PS *aps, PSC *apsc); // calc distance of a pulse shape with a collection
  Float_t Chi2Fast3D(PS *aps, PSC *apsc, float *Chi2Limits, float *chi2s, bool kFast); // compare a pulse shape with a collection
  Float_t Chi2Fast(PS *aps, PSC *apsc, float Chi2Limit, bool kFast); // compare a pulse shape with a collection
  Float_t Chi2seg(const float *spulse, const float *PSCspulse); // compare two segment pulse shape

  // EventHits
  void SortEventHits(); // sort EventHits increasing iconfig, irun, ientry
  int AddEventHits(EventHits* fHits);
  int FindiEvtHit(int iconfig, int irun, int ientry, int &istart);
  EventHits* FindEventHits(int i){ return fEventHits->at(i);}
  int GetEventHitsSize(){ return fEventHits->size();}

  // Tracking
  void TrackingLoop(); // tracking fHits
  void Tracking(); // loop all events
  
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
  void SetWithPS(bool val){ kWithPS = val;}
  void SetGroupPos(bool val){ kGroupPos = val;}

  
private:
  int Detid = -1;
  AGATAgeo *agatageo;

  Double_t MaxMemUsage = 50; // max memory usage %
  bool kPSA = false;
  bool kWithPS = true;
  bool kGroupPos = false;
  
  float MaxChi2;
  float MaxChi2s[3];

  // Pulse Shape Collection and Hit Collection
  TFile *pscfile[MaxNDets];
  TTree *psctree[MaxNDets][NSeg];

  int PSClimit[MaxNDets]; //PSC number limit
  bool kAddNewPSC[MaxNDets][NSeg]; // flag if new PSC added
  vector<PSC> fPSC[MaxNDets][NSeg]; // Pulse Shape Collection storage
  vector<HitCollection*>* fHCs[MaxNDets][NSeg]; // Hit Collection storage
  vector<Int_t> freeHCs[MaxNDets][NSeg]; // idx of empty fHCs
  vector<Int_t> divHCs[MaxNDets][NSeg]; // idx of divided fHCs
  vector<HitCollection*>* fAllHCs; // Hit Collection storage
  atomic_int ihc;
  mutex PSCmtx[MaxNDets][NSeg]; // fPSC lock for threads
  mutex AllHCmtx;

  
  // EventHits
  int nConfig = 0;
  vector<int> MinRun;
  vector<int> MaxRun;
  vector<EventHits*>* fEventHits; //fHits in Events for tracking
  atomic_int atomrun;
  atomic_int ievt;
  atomic_int NevtsTotal;
  mutex EvtHitsmtx;

  
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

  int Marker;
  float chi2limit[3];  // max chi2 for group PS
  float calpos[3];     // init selfcalib interaction position in labframe float(3)
  float cadpos[3];     // init selfcalib interaction position in detframe float(3)
  float calpos2[3];    // selfcalib interaction position in labframe float(3)
  float cadpos2[3];    // selfcalib interaction position in detframe float(3)
#ifdef REALPOS
  float labpos[3];     // average interaction position in labframe float(3)
  float detpos[3];     // average interaction position in detframe float(3)  
  float dist;          // dist calpos - labpos
  float dist2;         // dist calpos2 - labpos
  float cpos[3];                 // first pulse shape position in detframe float(3)
#endif

  int   cpulsehits;              // gamma hits number of first pulse
  float cpulse[NSeg_comp][NSig]; // first pulse shape for comparison
  float segwgt[NSeg_comp];       // first pulse segment weight for comparison
  float spulse[NSegCore][NSig];  // average pulse shape

  int npaths;

  // set tree branch
  void InitTreeWrite(TTree *tree);
  void InitTreeRead(TTree *tree);
  void ClosePSCFiles();
  
  Int_t NDets; // number of detectors from LookUpTable

  // counter
  atomic_int cPStotal;
  atomic_int cPSCtotal;
  atomic_int cPSCmem;
  atomic_int cPSCfile;
  atomic_int maxnhits;
  atomic_int cPaths;
  atomic_int cHits;
  atomic_int cHCs;

  time_t start, stop;

};

#endif // #ifndef AGATA_HH
