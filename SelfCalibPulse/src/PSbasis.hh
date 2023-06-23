#ifndef PSBASIS_HH
#define PABASIS_HH

#include "TVector3.h"
#include "TMatrixD.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Global.hh"

class PSbasis { // input PSbasis

public:
  PSbasis();
  PSbasis(int detid);
  virtual ~PSbasis();

  void ReadPSbasis();
  Int_t GetPS(int itype, TMatrixD pos, double energy, int &seg, TMatrixD &spulse);

  
private:
  const int fstep = 2; // mm fine grid
  static const int GridMaxSteps = 50;

  int Detid = -1;
  
  bool kextrapol = true;

  double range[NType][3][2];
  int imap[NType][GridMaxSteps][GridMaxSteps][GridMaxSteps];

  double threshold = 20; // threshold 20 keV for fired segment
  double threshold2 = 300; // threshold 300 keV for Compton scattering

  // db
  vector<Int_t>    dbseg[NType];
  vector<TMatrixD> dbpos[NType];
  vector<TMatrixD> dbspulse[NType]; // core at the end


  void swap(vector<double> &v,int m, int l){
    double temp;

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

};

#endif /* PSBASIS_HH  */
