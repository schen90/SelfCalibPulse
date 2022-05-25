/*
 * Path.cc
 *
 *  Created on: Jun 7, 2014
 *      Author: schigum
 */

#include "TMath.h"
#include <iostream>
#include "Path.hh"
#include "Global.hh"

Float_t Path::weight_norm = 0;

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

#ifndef SHORT
  incidentEnergy_exact = incE;
  depositedEnergy_exact = depE;
  comptonAngle_exact = CalcComptonAngle(incE, depE);

  incidentEnergy = incER;
  depositedEnergy = depER;
#endif
  comptonAngle = CalcComptonAngle(incER, depER);

  weight = 0;
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

void Path::initHC() {
  vector<HitCollection*>* hc0 = hit0->GetHitCollections();
  vector<HitCollection*>* hc1 = hit1->GetHitCollections();
  vector<HitCollection*>* hc2 = hit2->GetHitCollections();

  TVector3 pos0 = hit0->GetPosition();
  TVector3 pos1 = hit1->GetPosition();
  TVector3 pos2 = hit2->GetPosition();

  Double_t dist1 = (pos0 - pos1).Mag();
  Double_t dist2 = (pos1 - pos2).Mag();

  Float_t hc0w = 0;
  for (HitCollection* hc : *hc0) hc0w += hc->GetSize();
  hc0w = log(hc0w / hc0->size()) + 1; // why this way of weight????

  Float_t hc1w = 0;
  for (HitCollection* hc : *hc1) hc1w += hc->GetSize();
  hc1w = log(hc1w / hc1->size()) + 1;

  Float_t hc2w = 0;
  for (HitCollection* hc : *hc2) hc2w += hc->GetSize();
  hc2w = log(hc2w / hc2->size()) + 1;

  if (weight_norm == 0) {
    weight_norm = sqrt(dist1 * dist2 * hc0w * hc1w * hc2w);
  }
  weight = sqrt(dist1 * dist2 * hc0w * hc1w * hc2w) / weight_norm; // start weight = 1
}

Float_t Path::CalcChi2() {
  Hit *ahit;
  return CalcChi2(ahit);
}

Float_t Path::CalcChi2(Hit *fithit) {
  if(ignore)  return 0;

  TVector3 pos0 = hit0->GetPosition();
  TVector3 pos1 = hit1->GetPosition();
  TVector3 pos2 = hit2->GetPosition();

  // if used in fit
  if(!(!fithit)){
    vector<HitCollection*>* hcs = fithit->GetHitCollections();

    TVector3 average(0, 0, 0);
    if (unlikely(hcs->size() == 0)) {
      cout << "Error: not part of any HitCollection. Coords: " << endl;
    
    } else if (likely(hcs->size() == 1)) {// why likely????
      average = hcs->at(0)->GetFitPosition();

    } else {
      for (unsigned int i = 0; i < hcs->size(); i++) {
	average += hcs->at(i)->GetFitPosition();
      }
      average *= 1.0 / hcs->size();
    }

    if(fithit==hit0) pos0 = TVector3(average.X(), average.Y(), average.Z());
    if(fithit==hit1) pos1 = TVector3(average.X(), average.Y(), average.Z());
    if(fithit==hit2) pos2 = TVector3(average.X(), average.Y(), average.Z());
  }

  
  Float_t foundAngle = GetAngle(pos1 - pos0, pos2 - pos1);
  angleDiff = foundAngle - comptonAngle;

  Float_t diff = angleDiff; // diff in degree

  return fabs(diff);
}

void Path::RegisterWithHCs(){
  vector<HitCollection*>* hcs;

  if(!(hit0->GetDetectorID()<0)){ // not source hit
    hcs = hit0->GetHitCollections();
    for(HitCollection* hc : *hcs){ hc->AddPath(this); hc->AddPathHit(hit0); }
  }

  hcs = hit1->GetHitCollections();
  for(HitCollection* hc : *hcs){ hc->AddPath(this); hc->AddPathHit(hit1);}

  hcs = hit2->GetHitCollections();
  for(HitCollection* hc : *hcs){ hc->AddPath(this); hc->AddPathHit(hit2);}

}
