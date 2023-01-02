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

  int   divzone[NSeg_comp]; // divided zone
  vector<int> dividx;  // idx of daughter PSCs

  float devabscut[NSeg_comp][2];
  //float maxdevabs[NSeg_comp];     // max deviation (abs) value
  vector<float> devabs[NSeg_comp];     // all deviation (abs) values
  vector<float> dev[NSeg_comp];     // all deviation (abs) values
  float devsigma[NSeg_comp];      // standard deviation of compared segment

  float spulse[NSegCore][NSig];   // average pulse shape

  int   segcmp[NSeg_comp];        // segments for comparison
  float cpulse[NSeg_comp][NSig];  // average pulse shape for comparison
};

#endif /* PSC_HH */
