/*
 * HitCollection.cc
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#include "HitCollection.hh"

HitCollection::HitCollection(Int_t detid, Int_t segid, Int_t pscid,
			     Float_t* lpos, Float_t* cpos) {
  det = detid;
  seg = segid;
  PSCid = pscid;

  for(int i=0; i<3; i++){
    labpos[i] = lpos[i];
    calpos[i] = cpos[i];
    initpos[i] = cpos[i];
  }
  //calpos = cpos;
  
  hits = new vector<Hit*>();
  paths = new vector<Path*>();
  phits = new vector<Hit*>();

  Marker = 1;
}

HitCollection::~HitCollection() {
}

void HitCollection::RegisterWithHits() {
  for (UInt_t i = 0; i < hits->size(); i++) {
    hits->at(i)->AddHitCollection(this);
  }
}

void HitCollection::LockPaths(){
  for(Path* path : *paths) path->LockPath();
}

void HitCollection::UnlockPaths(){
  for(Path* path : *paths) path->UnlockPath();
}
