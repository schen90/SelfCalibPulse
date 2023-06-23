#ifndef PSC_HH
#define PSC_HH

#include <vector>

// structure for Pulse Shape Collection
class PSC {

public:
  PSC(Int_t detid, Int_t segid){ det=detid; seg=segid;}
  ~PSC(){;}

  int   det;
  int   seg;
  int   index;         // index for PSCid
  int   nhits;         // number of hits in group

  int   divzone[NCOMP]; // divided zone
  int   divdir[3];     // previous divide direction
  vector<int> dividx;  // idx of daughter PSCs

  float labpos[3];     // average interaction position in labframe float(3)
  float detpos[3];     // average interaction position in detframe float(3)

  float devabscut[NCOMP][2];
  //float maxdevabs[NCOMP];     // max deviation (abs) value
  vector<float> devabs[NCOMP];     // all deviation (abs) values
  vector<float> dev[NCOMP];    // all deviation (abs) values
  float devsigma[NCOMP];       // standard deviation of compared segment

  float spulse[NCHAN][BSIZE];  // average pulse shape

  int   segcmp[NCOMP];         // segments for comparison
  float cpulse[NCOMP][BSIZE];  // average pulse shape for comparison
};

#endif /* PSC_HH */
