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
#define NTHREADS2 40
#define ONECLUST // make one cluster in tracking
#define TRACKINGTREE // ouput tracking results
//#define PSA // PSA to assign initial pos
#define RADIUS0 2 // mm, initial PSC size from PSA
#define MINUIT2
#define DIFFTOTE 10 // total energy match source, maxdiff 10keV
//#define MULTISEG // include multi-segment events
//#define CHECKTRACK // check track using OFT tacking
//#define SINGLEHIT
#define PSCEMIN 0. // keV, PSC greate with PS Energy > PSCEMIN
#define MINHITS 10    // min nhits for a good HC
#define MAXHITS 500  // max nhit for a HC
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
#define NSig 60
#define NSig_comp 56
#define LOOP_SSE4_seg 14
#define LOOP_SSE8_seg 7

//#define SSE_M256
#define CHI2        3
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
  float energy; // core energy

  float opulse[NSegCore][NSig]; // original pulse shape
};

typedef struct PS PS;

// structure for PSA basis
struct PSAbasis{
  float pos[3];                  // pos in det frame
  float spulse[NSeg_comp][NSig]; // pulse shape for comparison
};

typedef struct PSAbasis PSAbasis;

#endif // ifndef GLOBAL_HH
