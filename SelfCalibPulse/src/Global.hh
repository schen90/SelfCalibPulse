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
#include <iostream>
#include <x86intrin.h>

using namespace std;

#define GB 1073741824. // size of 1GB 1024*1024*1024
#define MaxMemoryUsage 100.
#define NTHREADS 50
#define ONECLUST // make one cluster in tracking
//#define NOISE 1000000
//#define TRACKINGTREE // ouput tracking results
//#define MULTISEG // include multi-segment events
//#define ADDPS // input G4Tree noPS, addPS from db

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

#define FIXED_METRIC_NONE   0
#define FIXED_METRIC_ABSVAL 1
#define FIXED_METRIC_SQUARE 2
#define FIXED_METRIC_1SQRT  3
#define FIXED_METRIC_2SQRT  4

#define FIXED_METRIC FIXED_METRIC_ABSVAL

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
  
  float labpos[3]; // lab position
  float detpos[3]; // det position
  
  float opulse[NSegCore][NSig]; // original pulse shape
};

typedef struct PS PS;

#endif // ifndef GLOBAL_HH
