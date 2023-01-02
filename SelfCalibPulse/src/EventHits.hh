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
  EventHits(vector<float> sourceE, vector<TVector3> sourcePos){
    fhits = new vector<Hit*>();
    SourcePos = sourcePos;
    SourceE = sourceE;
    bestis = -1;
  }
  virtual ~EventHits(){
    for(Hit* ah : *fhits) delete ah;
    delete fhits;
  }

  void Add(Hit *ah){ fhits->push_back(ah);}
  Hit* Get(int i){ return fhits->at(i);}
  int  Size(){ return fhits->size();}
  vector<Hit*>* GetfHits(){ return fhits;}

  int      GetNSource(){ return SourceE.size();}
  TVector3 GetSourcePos(int i){ return SourcePos[i];}
  float    GetSourceE(int i){ return SourceE[i];}

  void SetBestis(int val){ bestis = val;}
  int  GetBestis(){ return bestis;}
  
  void SetIdx(int val1, int val2, int val3){ icfg=val1; irun=val2; ievt=val3;}
  void GetIdx(int &val1, int &val2, int &val3){ val1=icfg; val2=irun; val3=ievt;}

#ifdef DIFFTOTE
  float Etot;
#endif
  
private:
  // idx
  int icfg;
  int irun;
  int ievt;
  
  vector<Hit*>* fhits; // hits in one event

  vector<TVector3> SourcePos; // source position in lab frame
  vector<float> SourceE; // source energy keV, -1 for unknown E
  int bestis;

};

#endif /* EVENTHITS_HH */
