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

  void MakeSegmentMap();
  void GetNextSegs(Int_t iseg, Int_t *fseg);
  Float_t GetSegWeight(Int_t iseg, Int_t iiseg){ return SegWeight[iseg][iiseg];}
  
  bool CheckBounds(int detid, int segid, const double *lpos);
  void GetChi2sLimit(int detid, const double *dpos, float chi2slimit[]);

private:
  const int fstep = 2; // mm fine grid
  static const int GridMaxSteps = 50;

  // grid pos and seg
  double GridRange[NType][3][2];
  int gridimap[NType][GridMaxSteps][GridMaxSteps][GridMaxSteps];
  vector<Int_t> GridSeg[NType];
  Int_t    NSegGrid[NType][NSEGS];
  TMatrixD LocalSegPos[NType][NSEGS];

  float gridchi2smap[NType][GridMaxSteps][GridMaxSteps][GridMaxSteps][3];

  // transform matrix
  TMatrixD Rt[MaxNDets];
  TMatrixD Rt2[MaxNDets];
  TMatrixD Tr[MaxNDets];

  // segment position
  TMatrixD SegPos[MaxNDets][NSEGS];

  // neighbor segment
  Int_t NextSec[NSEGS][2]; // next sector
  Int_t NextSli[NSEGS][2]; // next slice
  Float_t SegWeight[NSEGS][NCHAN]; // weight in comparison
  
  Int_t NDets; // number of detectors from LookUpTable

  
};

#endif
