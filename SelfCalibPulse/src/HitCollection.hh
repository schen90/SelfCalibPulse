/*
 * HitCollection.hh
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#ifndef HITCOLLECTION_HH
#define HITCOLLECTION_HH

#include <vector>
#include "TVector3.h"
#include "Hit.hh"
#include "Path.hh"
#include <thread>

class Hit;
class Path;

using namespace std;

class HitCollection {

public:
  HitCollection(Int_t detid, Int_t segid, Int_t pscid, Float_t *lpos, Float_t *cpos);
  virtual ~HitCollection();

  int GetDet(){ return det;}
  int GetSeg(){ return seg;}
  int GetPid(){ return PSCid;}
  int GetGid(){ return gid;}

  void SetGid(int val){ gid = val;}

  void SetPosition(Float_t* posval){ for(int i=0; i<3; i++) calpos[i]=posval[i];}
  void SetInitPosition(Float_t* posval){ for(int i=0; i<3; i++) initpos[i]=posval[i];}
  void SetRealPosition(Float_t* posval){ for(int i=0; i<3; i++) labpos[i]=posval[i];}

  void SetFitPosition(Float_t* posval){
#ifdef NTHREADS2
    FitThreadID = this_thread::get_id();
#endif
    for(int i=0; i<3; i++) fitpos[i]=posval[i];
  }

  TVector3 GetPosition(){ return TVector3(calpos[0], calpos[1], calpos[2]);}
  TVector3 GetInitPosition(){ return TVector3(initpos[0],initpos[1],initpos[2]);}
  TVector3 GetRealPosition(){ return TVector3(labpos[0],labpos[1],labpos[2]);}

  TVector3 GetFitPosition(){
#ifdef NTHREADS2
    if(this_thread::get_id()!=FitThreadID)
      return TVector3(calpos[0], calpos[1], calpos[2]);
#endif
    return TVector3(fitpos[0], fitpos[1], fitpos[2]);
  }
  
  void AddHit(Hit* hit){  hits->push_back(hit);}
  void RemoveHit(Hit* hit){
    for(int i=0; i<hits->size(); i++){
      if(hit==hits->at(i)) hits->erase(hits->begin()+i);
    }
  }
  vector<Hit*>* GetHits(){ return hits;}
  unsigned int GetSize(){ return hits->size();}

  void AddPath(Path* p){ paths->push_back(p);}
  vector<Path*>* GetPaths(){ return paths;}
  void AddPathHit(Hit* hit){ phits->push_back(hit);}
  vector<Hit*>* GetPathHits(){ return phits;}

  void LockPaths();
  void UnlockPaths();

  void RegisterWithHits();

  void Clear(){
    for(int ix=0; ix<3; ix++){
      labpos[ix] = 0;
    }
    hits->clear();
    paths->clear();
    phits->clear();
    return;
  }
  
#ifdef NTHREADS2
  //mutex mtx; // lock for threads
  thread::id FitThreadID;
#endif

  float MaxChi2s[3]; // chi2 range for the group
  int Marker;
  
private:
  int   det;
  int   seg;
  int   PSCid; // id in fHCs
  int   gid; // global id in fAllHCs

  float calpos[3];  // calib position in lab frame
  float initpos[3]; // initial position in lab frame
  float labpos[3];  // real position in lab frame

  float fitpos[3];  // position used in fitting
  
  vector<Hit*>* hits; // a group of hits with similar PS
  vector<Path*>* paths; //
  vector<Hit*>* phits; // the hit belong to the paths
};

#endif /* HITCOLLECTION_HH */
