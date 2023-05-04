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
    bestis.push_back(-1);
  }
  virtual ~EventHits(){
    for(Hit* ah : *fhits) delete ah;
    delete fhits;
  }

  void Add(Hit *ah){ fhits->push_back(ah); sign.push_back(-1);}
  Hit* Get(int i){ return fhits->at(i);}
  int  Size(){ return fhits->size();}
  vector<Hit*>* GetfHits(){ return fhits;}
  vector<int> GetSign(){ return sign;}

  int      GetNSource(){ return SourceE.size();}
  TVector3 GetSourcePos(int i){ return SourcePos[i];}
  float    GetSourceE(int i){ return SourceE[i];}

  void SignClust(int iclust, int ihit){ sign[ihit] = iclust;}
  int  GetClust(int ihit){ return sign[ihit];}

  int  GetNClust(){ return bestis.size();}
  void SetBestis(int iclust, int val){ 
    while( iclust>bestis.size()-1 ) bestis.push_back(-1);
    bestis[iclust] = val;
  }
  int  GetBestis(int iclust){
    if( iclust>bestis.size()-1 ){
      cerr<<"iclust = "<<iclust<<" over ClustNumebr = "<<bestis.size()<<" !!!"<<endl;
      return -999;
    }
    return bestis[iclust];
  }
  
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
  vector<int> sign; // sign to clust

  vector<TVector3> SourcePos; // source position in lab frame
  vector<float> SourceE; // source energy keV, -1 for unknown E
  vector<int> bestis; // best match source id for each cluster

};

#endif /* EVENTHITS_HH */
