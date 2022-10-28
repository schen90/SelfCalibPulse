#ifndef HIT_HH
#define HIT_HH

#include <vector>
#include "TVector3.h"
#include <iostream>
#include "Global.hh"

using namespace std;

class Hit {

public:
  Hit(int hitdet, int hitseg, float E, TVector3 hitpos, TVector3 initpos); // keV and mm
  Hit(TVector3 sourcepos); // if hit is the source, det = -1, seg = -1
  Hit(); // if hit is the source, det = -1, seg = -1
  virtual ~Hit();

  void SetInterid(int val){ interid = val;}
  int GetInterid(){ return interid;}
  
  int GetDet(){ return det;}
  int GetSeg(){ return seg;}
  Float_t GetE(){ return depE;} // in keV

  TVector3 GetPosition(){ return TVector3(calpos[0],calpos[1],calpos[2]);}
  TVector3 GetRealPosition(){ return TVector3(labpos[0],labpos[1],labpos[2]);}

  void SetPSCID(int val){ PSCid = val;}
  void SetPosition(TVector3 pos){ calpos[0] = pos.X(); calpos[1] = pos.Y(); calpos[2] = pos.Z();}
  
#ifdef NOISE
  void SetNoiseIdx(int val){ noiseidx = val;}
  Int_t GetNoiseIdx(){ return noiseidx;}
  void SetNoiseIdxShift(int val){ noiseidxshift = val;}
  Int_t GetNoiseIdxShift(){ return noiseidxshift;}
#endif

private:
  int   det;
  int   seg;
  int   interid; // interaction id in a event

  Float_t depE; // deposit Energy in the hit, keV

  int     PSCid;      // selfcalib PSC id in seg
  Float_t calpos[3];  // selfcalib position in labframe
  Float_t labpos[3];  // simulated position in labframe

#ifdef NOISE
  int noiseidx;
  int noiseidxshift;
#endif
};

#endif /* HIT_HH */
