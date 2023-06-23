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
#define NTHREADS 10
#define NTHREADS2 80
//#define ONECLUST // make one cluster in tracking
#define TRACKINGTREE // ouput tracking results
//#define PSA // PSA to assign initial pos
#define RADIUS0 2 // mm, initial PSC size from PSA
#define MINUIT2
#define DIFFTOTE 10 // total energy match source, maxdiff 10keV; negative value means only accept such events
//#define MULTISEG // include multi-segment events
//#define CHECKTRACK // check track using OFT tacking
//#define SINGLEHIT
#define ADDPS // input G4Tree noPS, addPS from db
#define NOISE 1000000
#define PSCEMIN 300. // keV, PSC greate with PS CoreEnergy > PSCEMIN
#define MINHITS 10    // min nhits for a good HC
#define MAXHITS 1000  // max nhit for a HC
#define SHORT

// agata
const int MaxNDets = 180;
const int NType = 3;
const int NSLIC = 6;
const int NSECT = 6;
const int NSEGS = NSLIC*NSECT;
const int NCHAN = NSEGS + 1;
const int INDCC = NCHAN - 1;  // 36 index of CC
const int TSTEP = 10; // time step of the experimental data and of the internal basis signals (ns)

const int BSIZE = 56; // number of points of the internal basis
const int BZERO =  0; // "T=0" sample of the calculated traces

const int DSIZE = 60; // length of experimental data
const int DZERO = 11; // "T=0" sample of the experimental traces

const int NCOMP = 6;

const int LOOP_SAMP = 40;
const int LOOP_4SSE = LOOP_SAMP/4;
const int LOOP_8SSE = LOOP_SAMP/8;

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
  int   nhits;
  int   interid;  // interaction id in a event
  vector<float> hiteng; // hit energy
  float energy; // core energy

  float labpos[3]; // lab position
  float detpos[3]; // det position

  float opulse[NCHAN][BSIZE];   // original pulse shape
  float cpulse[NCOMP][BSIZE];  // adapted pulse shape for comparison
};

typedef struct PS PS;

// structure for PSA basis
struct PSAbasis{
  float pos[3];                  // pos in det frame
  float spulse[NCOMP][BSIZE]; // pulse shape for comparison
};

typedef struct PSAbasis PSAbasis;


// structure for configfile
struct Config{
  int              NSource;
  vector<float>    fSourceE;
  vector<TVector3> fSourcePos;
  string           path;
  int              MinRun;
  int              MaxRun;
  long long        Nevts;
  vector<int>      skipdetid;
};


#endif // ifndef GLOBAL_HH
