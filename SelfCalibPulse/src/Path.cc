#include "TMath.h"
#include <iostream>
#include "Path.hh"
#include "Global.hh"

Float_t CalcComptonAngle(Float_t Ein, Float_t Edep){
  float ME_keV = mec2*1000;
  Float_t cosa = 1. + ME_keV / Ein - ME_keV / (Ein - Edep);
  if(unlikely( cosa > 1 ))
    cosa = 1;
  else if(unlikely( cosa < -1 ))
    cosa = -1;

  return (180. / TMath::Pi() * acos(cosa));
}


Float_t GetAngle(TVector3 v, TVector3 w){
  return (180. / TMath::Pi() * (v.Angle(w)));
}


Path::Path(Hit* i0, Hit* i1, Hit* i2,
           Float_t incE, Float_t depE,
           Float_t incER, Float_t depER) {

  hit0 = i0;
  hit1 = i1;
  hit2 = i2;

  incidentEnergy_exact = incE;
  depositedEnergy_exact = depE;
  comptonAngle_exact = CalcComptonAngle(incE, depE);

  incidentEnergy = incER;
  depositedEnergy = depER;

  comptonAngle = CalcComptonAngle(incER, depER);

  ignore = false;
}

Path::~Path() {
}

void Path::init(Float_t threshold) {
  TVector3 pos0 = hit0->GetPosition();
  TVector3 pos1 = hit1->GetPosition();
  TVector3 pos2 = hit2->GetPosition();

  TVector3 link1 = pos1 - pos0;
  TVector3 link2 = pos2 - pos1;

  Float_t foundAngle = GetAngle(link1, link2);
  angleDiff = foundAngle - comptonAngle;

  Float_t diff = angleDiff / comptonAngle; // diff in %
  
  if (fabs(diff) > threshold || link1.Mag2() < 1 || link2.Mag2() < 1) { // not real Compton
    ignore = true;
  } else {
    ignore = false;
  }
}


Float_t Path::CalcChi2() {
  if(ignore)  return 0;

  TVector3 pos0 = hit0->GetPosition();
  TVector3 pos1 = hit1->GetPosition();
  TVector3 pos2 = hit2->GetPosition();

  Float_t foundAngle = GetAngle(pos1 - pos0, pos2 - pos1);
  angleDiff = foundAngle - comptonAngle;

  Float_t diff = angleDiff; // diff in degree

  return fabs(diff);
}
