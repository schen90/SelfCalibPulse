#include "Hit.hh"

Hit::Hit(int hitdet, int hitseg, float E, TVector3 pos, TVector3 initpos) {
  det = hitdet;
  seg = hitseg;
  depE = E;
  labpos[0] = pos.X();     labpos[1] = pos.Y();     labpos[2] = pos.Z();
  calpos[0] = initpos.X(); calpos[1] = initpos.Y(); calpos[2] = initpos.Z();
}

Hit::Hit(TVector3 sourcepos) {  // if hit is the source
  det = -1;
  seg = -1;
  depE = 0;
  labpos[0] = sourcepos.X(); labpos[1] = sourcepos.Y(); labpos[2] = sourcepos.Z();
  calpos[0] = labpos[0];     calpos[1] = labpos[1];     calpos[2] = labpos[2];
}

Hit::Hit() {  // if hit is the source
  Hit(TVector3(0,0,0));
}

Hit::~Hit() {
}


