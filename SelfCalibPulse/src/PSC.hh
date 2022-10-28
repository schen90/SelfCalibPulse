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

  float labpos[3];     // average interaction position in labframe float(3)
  float detpos[3];     // average interaction position in detframe float(3)

  int   segcmp[NSeg_comp];        // segments for comparison
  float devabscut[NSeg_comp];
  //float maxdevabs[NSeg_comp];     // max deviation (abs) value
  vector<float> devabs[NSeg_comp];     // all deviation (abs) values
  float spulse[NSegCore][NSig];   // average pulse shape
};

#endif /* PSC_HH */
