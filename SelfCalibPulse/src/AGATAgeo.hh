#ifndef AGATAGEO_HH
#define AGATAGEO_HH

#include <TVector3.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TTree.h>
#include <time.h>
#include <vector>
#include "TMath.h"

#include "Global.hh"

using namespace std;
  
class AGATAgeo {

public:
  AGATAgeo();
  virtual ~AGATAgeo();

  Int_t GetNDets(){ return NDets;}
  void LoadGrid(Int_t itype, string gridfile = "pulsedb/LibTrap_A001.root");
  
  void LoadMatrix(string LookUpTable = "LookUp/CrystalPositionLookUpTable");
  TMatrixD Lab2DetPos(Int_t idet, TMatrixD LabPos);
  TMatrixD Det2LabPos(Int_t idet, TMatrixD DetPos);

  void MakeSegPos();
  TMatrixD GetSegPos(Int_t idet, Int_t seg){ return SegPos[idet][seg];}
  TMatrixD GetLocalSegPos(Int_t itype, Int_t seg){ return LocalSegPos[itype][seg];}

  void MakeSegMap(int neighbours);
  char *GetMask(int seg){ return hmask[seg];}

  bool CheckBounds(int detid, int segid, const double *lpos);

private:
  // grid pos and seg
  double GridRange[NType][3][2];
  int gridimap[NType][GridMaxSteps][GridMaxSteps][GridMaxSteps];
  vector<Int_t> GridSeg[NType];
  Int_t    NSegGrid[NType][NSeg];

  // transform matrix
  TMatrixD Rt[MaxNDets];
  TMatrixD Rt2[MaxNDets];
  TMatrixD Tr[MaxNDets];

  // segment position
  TMatrixD SegPos[MaxNDets][NSeg];
  TMatrixD LocalSegPos[NType][NSeg];

  char hmask[NSeg][NSegCore+1];
  
  Int_t NDets; // number of detectors from LookUpTable

  
};

#endif
