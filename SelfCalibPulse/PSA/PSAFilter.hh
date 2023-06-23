#ifndef PSAFILTER_HH
#define PSAFILTER_HH

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <string.h>

#include "Global.hh"

using namespace std;

class PSAFilter{
public:
  PSAFilter(string basisfile);
  ~PSAFilter();

  void MakeSegmentMap(int neighbours);
  void ReadPSAbasis(string basisfile);
  int  NearestPoint(float px, float py, float pz, int iseg, float& dist);

  void ProcessOneEvent(pointExp &PE);
  void PrepareEvent(pointExp &PE, bool doReset);
  void InitialSolution(pointExp &PE, bool keepPrevious);
  int  ProcessTheEvent(pointExp &PE);
  int  SearchFullGrid(pointExp &PE, int netChSeg, float netChEner);
  void AddToSolution(pointExp &PE, pointPsa &bestPoint, float scaleFact);
  float GetChi2(pointExp &PE);

  int  MakeSearchWave(pointExp &PE);
  int  MakeLocalMask(pointExp &PE);


  void GetPtPos(int netChSeg, int id, float *pos){
    if(id<0){
      for(int ix=0; ix<3; ix++) pos[ix] = averPt[netChSeg].pos[ix];
    }else{
      for(int ix=0; ix<3; ix++) pos[ix] = segPts[netChSeg][id].pos[ix];
    }
  }

private:

  const int fstep = 2; // mm fine grid

  vector<pointPsa> segPts[NSEGS];
  pointPsa averPt[NSEGS];
  int      averPtlist[NSEGS];
  char hmask[NSEGS][NCHAN+1];

  float Chi2InnerLoop(const float *pReal, const float *pBase, float bScale);

};

#endif
