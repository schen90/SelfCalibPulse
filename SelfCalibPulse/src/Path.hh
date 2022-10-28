#ifndef PATH_HH
#define PATH_HH

#include "Hit.hh"
#include "Global.hh"

class Hit;

using namespace std;

class Path { // three consecutive hits from a track, the smallest set to compare Compton angle

public:
  Path(Hit* i0, Hit* i1, Hit* i2, Float_t incE, Float_t depE, Float_t incER, Float_t depER);
  virtual ~Path();

  void init(Float_t threshold); // if angleDiff % > threshold, ignore
  bool isIgnored() { return ignore;}

  Float_t GetComptonAngle() { return comptonAngle;}

  Hit* GetHit0() { return hit0;}
  Hit* GetHit1() { return hit1;}
  Hit* GetHit2() { return hit2;}
  
  Float_t CalcChi2();

private:
  Hit* hit0;
  Hit* hit1;
  Hit* hit2;

  Float_t incidentEnergy_exact; // not really needed
  Float_t depositedEnergy_exact; // not really needed
  Float_t comptonAngle_exact; // not really needed
  
  Float_t incidentEnergy; //keV
  Float_t depositedEnergy; //keV

  Float_t comptonAngle;

  bool ignore;

  Float_t angleDiff;
};

#endif /* PATH_HH */
