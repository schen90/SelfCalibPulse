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

#include "Global.hh"
#include "AGATAgeo.hh"
#include "Hit.hh"
#include "EventHits.hh"
#include "PSC.hh"
#include "Path.hh"
#include "Tracker.hh"

using namespace std;

class AGATA {

public:
  AGATA();
  virtual ~AGATA();

  AGATAgeo* GetGeo(){ return agatageo;}

  // write/read files
  void InitPSC();                // initial PSC0 and PSC1
  void CopyPSC();                // Copy fPSC1 to fPSC0
  void ClearPSC1();              // Clear fPSC1
  void ReadPSCfiles(int detid);  // read file to PSC0
  void WritePSCfiles(int detid); // write PSC1 to file
  void GetPSCstat(int *PSCstat); // get Pulse Shape Collection state
  
  int FindPSC(int detid, int segid, float dpos[], int istart);

  int AddPStoPSC(PS *aps, int ipsc);

  // Fit Path
  int FitPath(Path *apath);
  
  // PSA to assign initial Hit pos
  TVector3 GetPSpos(int detid, int segid, PS *aps);
  Float_t Chi2seg(const float *spulse, const float *PSCspulse); // compare two segment pulse shape
  
  // Check Memory
  void SetMaxMemUsage(double value){ MaxMemUsage = value;}
  Double_t GetTotalSystemMemory();
  Double_t GetCurrentMemoryUsage();

private:
  AGATAgeo *agatageo;

  Double_t MaxMemUsage = 100; // max memory usage %
  
  // Pulse Shape Collection
  TFile *pscfile[MaxNDets];
  TTree *psctree[MaxNDets][NSeg];

  int   iPSC0;                         // PSC0 update times, 0: PSC0 empty
  vector<PSC*>* fPSC0[MaxNDets][NSeg]; // old Pulse Shape Collection for PSA
  vector<PSC*>* fPSC1[MaxNDets][NSeg]; // new Pulse Shape Collection for SelfCalib
  vector<TVector3> PSCpos[MaxNDets][NSeg]; // Pulse shape Collection position
  
  mutex PSCmtx[MaxNDets][NSeg]; // fPSC lock for threads


  // branch
  int   det;           // detector id
  int   seg;           // segment id
  int   nhits;         // number of hits in group
  
  float labpos[3];     // position in labframe float(3)
  float detpos[3];     // position in detframe float(3)

  float avelpos[3];     // average interaction position in labframe float(3)
  float avedpos[3];     // average interaction position in detframe float(3)

  float dist;          // dist avelpos - labpos
  
  float spulse[NSegCore][NSig];  // average pulse shape

  // set tree branch
  void InitTreeWrite(TTree *tree);
  void InitTreeRead(TTree *tree);

  Int_t NDets; // number of detectors from LookUpTable

  time_t start, stop;

  // counter
  atomic_int cPSCmem;
  atomic_int cPSCfile;
  atomic_int maxnhits;
  
};

#endif // #ifndef AGATA_HH
