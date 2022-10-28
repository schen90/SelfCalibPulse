#ifndef GLOBAL_HH
#define GLOBAL_HH

#include <stdlib.h>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <TMath.h>
#include <TVector3.h>
#include "TRandom.h"
#include "TSystem.h"

using namespace std;

#define GB 1073741824. // size of 1GB 1024*1024*1024
#define MaxMemoryUsage 100.
#define NTHREADS 6
#define NTHREADS2 30
#define ONECLUST // make one cluster in tracking
#define TRACKINGTREE // ouput tracking results
#define WITHPS //comment out to save memory
#define PSA // PSA to assign initial pos
#define MINUIT2
//#define MULTISEG // include multi-segment events
#define CHECKTRACK
//#define SINGLEHIT
#define REALPOS
#define ADDPS // input G4Tree noPS, addPS from db
#define NOISE 1000000
#define PSCEMIN 900. // keV, PSC greate with PS Energy > PSCEMIN
#define MINHITS 4    // min nhits for a good HC
#define MAXHITS 200  // half of max nhit for a HC
#define SHORT

// agata
#define NType 3
#define GridDist 2. // 2mm grid
#define GridMaxSteps 50
#define MaxNDets 180
#define NSli 6
#define NSec 6
#define NSeg 36
#define NSegCore 37
#define NSeg_comp 6
#define NSig 56
#define NSig_comp 56
#define LOOP_SSE4_seg 14
#define LOOP_SSE8_seg 7

//#define SSE_M256
#define CHI2        1
#define CHI2_SQ     0
#define CHI2_CHI2   1
#define CHI2_CHI2_2 2
#define CHI2_ABS    3
#define CHI2_FABS   4
#define CHI2_SQRT   5
#define CHI2_2SQRT  6

#define PI 3.1415927
#define rho_ge 5.32 /* g/cm3 */
#define Z_ge 32
#define A_ge 74
#define N_av 6.022e23
#define r_0 2.8179e-15
#define mec2 0.510998910 /* electron mass MeV */
#define Alpha 0.0072993 /* fine structure constant 1/137 */

#define SQ(x) ((x)*(x))
#define CB(x) ((x)*(x)*(x))
#define likely(x)    __builtin_expect(!!(x), 1)
#define unlikely(x)  __builtin_expect(!!(x), 0)


// structure for Pules Shape
struct PS{
  int   det;
  int   seg;
  int   nhits;
  int   interid;  // interaction id in a event
  vector<float> hiteng; // hit energy
  float energy; // core energy
#ifdef REALPOS
  float labpos[3]; // lab position
  float detpos[3]; // det position
#endif
  float apulse[NSeg_comp][NSig]; // adapted pulse shape for comparison
  float segwgt[NSeg_comp];       // weight for comparison
  float opulse[NSegCore][NSig]; // original pulse shape
};

typedef struct PS PS;

// structure for Pulse Shape Collection
struct PSC{
  int   det;           // detector id start from 0
  int   seg;           // segment id start from 0
  int   index;         // index for PSCid
  int   nhits;         // number of hits in group

  //float calpos[3];     // selfcalib interaction position in labframe float(3)
  //float cadpos[3];     // selfcalib interaction position in detframe float(3)
#ifdef REALPOS
  float labpos[3];     // average interaction position in labframe float(3)
  float detpos[3];     // average interaction position in detframe float(3)  
  float cpos[3];       // first pulse shape position in detframe float(3)
#endif

#ifdef WITHPS
  int   cpulsehits;              // gamma hits number of first pulse
  float cpulse[NSeg_comp][NSig];  // first pulse shape for comparison
  float segwgt[NSeg_comp];        // weight for comparison
  float spulse[NSegCore][NSig];  // average pulse shape
#endif
};

typedef struct PSC PSC;

// structure for PSA basis
struct PSAbasis{
  float pos[3];                  // pos in det frame
  float spulse[NSeg_comp][NSig]; // pulse shape for comparison
};

typedef struct PSAbasis PSAbasis;

#endif // ifndef GLOBAL_HH
