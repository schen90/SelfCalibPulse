#ifndef TRACKER_HH
#define TRACKER_HH

#include <vector>
#include <TVector3.h>

#include "Global.hh"
#include "Hit.hh"

using namespace std;

class Hit;

class Tracker {

public:
  Tracker(vector<Hit*>* hits);
  Tracker(vector<Hit*>* hits, double E);
  Tracker(vector<Hit*>* hits, double E, TVector3 sourcepos);
  virtual ~Tracker();

  vector<int> GetTrack(){ return track;}
  
  void SetSourcePos(float *posval){ sPos.SetXYZ(posval[0],posval[1],posval[2]);}
  void SetAlfaRed(double val){ alfared = val;}
  void OFTtracking();
  double Ftot(vector<TVector3> &pos, vector<int> &intid, vector<double> &energy,
	      double etotale, vector<int> &order);

  void Simpletracking();
  double CalcMinChi2(vector<TVector3> &pos, vector<int> &intid, vector<double> &energy,
		     double etotale, vector<int> &order);

  
private:
  vector<Hit*>* fHits; // Hits from one event
  int nhits;

  double EGamma; // keV gamma energy for calibration, -1 means unknown energy
  
  vector<int> track; // the most likely order of Hits

  TVector3 sPos; // cm source position
  
  double e[MaxNDets]; // MeV
  TVector3 Pos[MaxNDets]; // cm relative to source position
  double r[MaxNDets][MaxNDets], r_ge[MaxNDets][MaxNDets]; // cm
  double thetatest[MaxNDets][MaxNDets];

  // tracking parameters
  double radius; // center to front surface in cm
  double radius2; // center to back surface in cm
  double eres; // energy resolution MeV
  int kmax;
  double alfamin;
  double alfared; // threshold for angle acceptance
  double deltaalfa;
  int maxinteraction = 7; // max number of interaction in a cluster
  int nopair;
  double minprobtrack;
  int interactionpair;
  double probtestpair;
  double sigma_thet;

  
  void swap(vector<double> &v,int m, int l){
    double temp;

    temp = v[m];
    v[m] = v[l];
    v[l] = temp;
  }

  void swap(vector<TVector3> &v,int m, int l){
    TVector3 temp;

    temp = v[m];
    v[m] = v[l];
    v[l] = temp;
  }

  void swap(vector<int> &v, int m, int l){
    int temp;
    temp = v[m];
    v[m] = v[l];
    v[l] = temp;
  }

  void swap(int v[],int m, int l){
    int temp;

    temp = v[m];
    v[m] = v[l];
    v[l] = temp;
  }

};

#endif
