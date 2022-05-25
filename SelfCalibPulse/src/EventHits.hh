#ifndef EVENTHITS_HH
#define EVENTHITS_HH

#include <vector>
#include "TVector3.h"
#include <iostream>
#include "Global.hh"
#include "Hit.hh"

using namespace std;

class EventHits;

class EventHits {

public:
  EventHits(float sourceE, TVector3 sourcePos){
    fhits = new vector<Hit*>();
    SourcePos[0] = sourcePos.X();
    SourcePos[1] = sourcePos.Y();
    SourcePos[2] = sourcePos.Z();
    SourceE = sourceE;
  }
  virtual ~EventHits(){
    for(Hit* ah : *fhits) delete ah;
    delete fhits;
  }

  void Add(Hit *ah){ fhits->push_back(ah);}
  Hit* Get(int i){ return fhits->at(i);}
  int Size(){ return fhits->size();}
  vector<Hit*>* GetfHits(){ return fhits;}
  
  TVector3 GetSourcePos(){ return TVector3(SourcePos[0],SourcePos[1],SourcePos[2]);}
  float GetSourceE(){ return SourceE;}
  
  void SetIdx(int val1, int val2, int val3){ icfg=val1; irun=val2; ievt=val3;}
  void GetIdx(int &val1, int &val2, int &val3){ val1=icfg; val2=irun; val3=ievt;}
  
private:
  // idx
  int icfg;
  int irun;
  int ievt;
  
  vector<Hit*>* fhits; // hits in one event

  float SourcePos[3]; // source position in lab frame
  float SourceE; // source energy keV, -1 for unknown E

};

#endif /* EVENTHITS_HH */
