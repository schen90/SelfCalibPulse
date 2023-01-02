/*
 * Hit.cc
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#include "Hit.hh"

Hit::Hit(int hitdet, int hitseg, float E, TVector3 initpos) {
  level = 1;
  det = hitdet;
  seg = hitseg;
  depE = E;
  calpos[0] = initpos.X(); calpos[1] = initpos.Y(); calpos[2] = initpos.Z(); 

  hitCollections = new vector<HitCollection*>();
}

Hit::Hit(TVector3 sourcepos) {  // if hit is the source
  level = 1;
  det = -1;
  seg = -1;
  depE = 0;
  calpos[0] = sourcepos.X(); calpos[1] = sourcepos.Y(); calpos[2] = sourcepos.Z();

  hitCollections = new vector<HitCollection*>();  
}

Hit::Hit() {  // if hit is the source
  Hit(TVector3(0,0,0));
}

Hit::~Hit() {
  // do not delete the hitCollections vector
}


Int_t Hit::hasgoodHCs(int thres){
  int ngoodhcs = 0;
  for(HitCollection* ahc : *hitCollections)
    if(ahc->GetSize()>thres) ngoodhcs++;
  return ngoodhcs;
}


void Hit::CalcAveHCsPosition(int thres) {
  if(det<0) return; // skip source hit

  if (unlikely(hitCollections->size() == 0)) {

    return;
    
  } else if (likely(hitCollections->size() == 1)) {// why likely????
    TVector3 tmppos = hitCollections->at(0)->GetPosition();
    calpos[0] = tmppos.X(); calpos[1] = tmppos.Y(); calpos[2] = tmppos.Z();
    return;

  } else {
    TVector3 average(0, 0, 0);
    double weight = 0;
    for(HitCollection* ahc : *hitCollections){
      if(ahc->GetSize()>thres){
#define WEIGHT
#ifndef WEIGHT
	average += ahc->GetPosition();
	weight++;
#else
	average += ahc->GetPosition() * ahc->GetSize();
	weight += ahc->GetSize();
#endif
      }
    }
    if(weight>0) average *= 1.0 / weight;
    else         average = hitCollections->at(0)->GetPosition();
    calpos[0] = average.X(); calpos[1] = average.Y(); calpos[2] = average.Z();
    return;
  }
}
