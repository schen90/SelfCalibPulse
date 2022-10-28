#ifndef PSC_HH
#define PSC_HH

#include <vector>

// structure for Pulse Shape Collection
class PSC {

public:
  PSC(Int_t detid, Int_t segid){
    det=detid;
    seg=segid;
    nhits=0;

    for(int ix=0; ix<3; ix++)
      avelpos[ix] = avedpos[ix] = 0;

    for(int iseg=0; iseg<NSegCore; iseg++)
      for(int isig=0; isig<NSig; isig++)
	spulse[iseg][isig] = 0;
  }
  
  ~PSC(){;}

  int   det;
  int   seg;
  int   nhits;         // number of hits in group
  
  float labpos[3];     // position in labframe float(3)
  float detpos[3];     // position in detframe float(3)

  float avelpos[3];    // average interaction position in labframe float(3)
  float avedpos[3];    // average interaction position in detframe float(3)  

  float spulse[NSegCore][NSig];   // average pulse shape
};

#endif /* PSC_HH */
