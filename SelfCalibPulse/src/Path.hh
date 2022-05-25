/*
 * Path.hh
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#ifndef PATH_HH
#define PATH_HH

#include "Hit.hh"
#include "HitCollection.hh"
#include "Global.hh"
#include <mutex>

class Hit;
class HitCollection;

using namespace std;

class Path { // three consecutive hits from a track, the smallest set to compare Compton angle

public:
  Path(Hit* i0, Hit* i1, Hit* i2, Float_t incE, Float_t depE, Float_t incER, Float_t depER);
  virtual ~Path();

  void init(Float_t threshold); // if angleDiff % > threshold, ignore
  void initHC();
  bool isIgnored() { return ignore;}

  Hit* GetHit0() { return hit0;}
  Hit* GetHit1() { return hit1;}
  Hit* GetHit2() { return hit2;}

  //vector<HitCollection*>* GetHCs0() { return hit0->GetHitCollections();}
  //vector<HitCollection*>* GetHCs1() { return hit1->GetHitCollections();}
  //vector<HitCollection*>* GetHCs2() { return hit2->GetHitCollections();}

  Float_t CalcChi2();
  Float_t CalcChi2(Hit *fithit);

  void LockPath(){ mtx.lock();}
  void UnlockPath(){ mtx.unlock();}
  
  void RegisterWithHCs();

private:
  mutex mtx;
  static Float_t weight_norm;

  Hit* hit0;
  Hit* hit1;
  Hit* hit2;

#ifndef SHORT
  Float_t incidentEnergy_exact; // not really needed
  Float_t depositedEnergy_exact; // not really needed
  Float_t comptonAngle_exact; // not really needed

  Float_t incidentEnergy; //keV
  Float_t depositedEnergy; //keV
#endif
  Float_t comptonAngle;
  
  Float_t weight;
  bool ignore;

  Float_t angleDiff;
};

#endif /* PATH_HH */
