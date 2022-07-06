/*
 * Hit.hh
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#ifndef HIT_HH
#define HIT_HH

#include <vector>
#include "TVector3.h"
#include <iostream>
#include "Global.hh"
#include "HitCollection.hh"

using namespace std;

class HitCollection;

class Hit {

public:
  Hit(int hitdet, int hitseg, float E, TVector3 hitpos, TVector3 initpos); // keV and mm
  Hit(TVector3 sourcepos); // if hit is the source, det = -1, seg = -1
  Hit(); // if hit is the source, det = -1, seg = -1
  virtual ~Hit();

  void SetLevel(int val){ level = val;}
  int GetLevel(){ return level;}
  
  int GetDet(){ return det;}
  int GetSeg(){ return seg;}
  
  void ClearHitCollection(){ hitCollections->clear();}
  void AddHitCollection(HitCollection* hc){ hitCollections->push_back(hc);}
  bool InHitCollection(HitCollection* hc){
    for(int i=0; i<hitCollections->size(); i++)
      if(hc==hitCollections->at(i))
	return true;
    return false;
  }
  void RemoveHitCollection(HitCollection* hc){
    for(int i=0; i<hitCollections->size(); i++)
      if(hc==hitCollections->at(i))
	hitCollections->erase(hitCollections->begin()+i);
  }
  vector<HitCollection*>* GetHitCollections(){ return hitCollections;}
  Int_t hasHitCollection(){ return hitCollections->size();}
  Int_t hasgoodHCs(int thres);

  Float_t GetE(){ return depE;} // in keV
  
  TVector3 GetPosition(){ return TVector3(calpos[0],calpos[1],calpos[2]);}
  TVector3 GetRealPosition(){ return TVector3(labpos[0],labpos[1],labpos[2]);}

  void CalcAveHCsPosition(); // calc. calpos

  void SetDetectorID(Int_t val){ det=val;}
  void SetSegmentID(Int_t val){ seg=val;}

  Int_t GetDetectorID(){ return det;}
  Int_t GetSegmentID(){ return seg;}

#ifdef NOISE
  void SetNoiseIdx(int val){ noiseidx = val;}
  Int_t GetNoiseIdx(){ return noiseidx;}
  void SetNoiseIdxShift(int val){ noiseidxshift = val;}
  Int_t GetNoiseIdxShift(){ return noiseidxshift;}
#endif
  
private:
  int level; // 1:init  0:cannot group with others  2:check for divide

  int   det;
  int   seg;
  vector<HitCollection*>* hitCollections; //one hit can be assigned to several hitcollections, pos set as the ave of all hitcollections

  Float_t depE; // deposit Energy in the hit, keV
  Float_t calpos[3];  // selfcalib position in labframe, set as the average of all PSC
  Float_t labpos[3];  // simulated position in labframe

#ifdef NOISE
  int noiseidx;
  int noiseidxshift;
#endif
};

#endif /* HIT_HH */
