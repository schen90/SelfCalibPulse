/*
 * Hit.cc
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#include "Hit.hh"

Hit::Hit(int hitdet, int hitseg, float E, TVector3 pos, TVector3 initpos) {
  level = 1;
  det = hitdet;
  seg = hitseg;
  depE = E;
  labpos[0] = pos.X();     labpos[1] = pos.Y();     labpos[2] = pos.Z();
  calpos[0] = initpos.X(); calpos[1] = initpos.Y(); calpos[2] = initpos.Z(); 

  hitCollections = new vector<HitCollection*>();
}

Hit::Hit(TVector3 sourcepos) {  // if hit is the source
  level = 1;
  det = -1;
  seg = -1;
  depE = 0;
  labpos[0] = sourcepos.X(); labpos[1] = sourcepos.Y(); labpos[2] = sourcepos.Z();
  calpos[0] = labpos[0];     calpos[1] = labpos[1];     calpos[2] = labpos[2];

  hitCollections = new vector<HitCollection*>();  
}

Hit::Hit() {  // if hit is the source
  Hit(TVector3(0,0,0));
}

Hit::~Hit() {
  // do not delete the hitCollections vector
}

void Hit::CalcAveHCsPosition() {
  if(det<0) return; // skip source hit

  if (unlikely(hitCollections->size() == 0)) {
    // keep original position at segment center
    //cout << "Error: not part of any HitCollection. Coords: "
    //     << labpos.x() << " " << labpos.y() << " " << labpos.z()
    //	   << " " << det << " " << seg << " " << endl;
    //calpos = labpos;
    return;
    
  } else if (likely(hitCollections->size() == 1)) {// why likely????
    TVector3 tmppos = hitCollections->at(0)->GetPosition();
    calpos[0] = tmppos.X(); calpos[1] = tmppos.Y(); calpos[2] = tmppos.Z();
    return;

  } else {
    TVector3 average(0, 0, 0);
    for (unsigned int i = 0; i < hitCollections->size(); i++) {
      average += hitCollections->at(i)->GetPosition();
    }
    average *= 1.0 / hitCollections->size();
    calpos[0] = average.X(); calpos[1] = average.Y(); calpos[2] = average.Z();
    return;
  }
}
