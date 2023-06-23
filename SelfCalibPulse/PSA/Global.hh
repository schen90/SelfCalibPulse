#ifndef GLOBAL_HH
#define GLOBAL_HH

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;

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

const int LOOP_SAMP = 40;          // number of samples used to calculate FoM for each segment
const int LOOP_4SSE = LOOP_SAMP/4; // processed in groups of 4

#define CHI2        3
#define CHI2_SQ     0
#define CHI2_CHI2   1
#define CHI2_CHI2_2 2
#define CHI2_ABS    3
#define CHI2_FABS   4
#define CHI2_SQRT   5
#define CHI2_2SQRT  6

#define GB 1073741824. // size of 1GB 1024*1024*1024

/*
GeDataBaseLNL2022={
 '00A': ['00', '$PSABASE/LibTrap_A006.dat'], '00B': ['01', '$PSABASE/LibTrap_B005.dat'], '00C': ['02', '$PSABASE/LibTrap_C001.dat'],
 '01A': ['03', '$PSABASE/LibTrap_A011.dat'], '01B': ['04', '$PSABASE/LibTrap_B006.dat'], '01C': ['05', '$PSABASE/LibTrap_C012.dat'],
 '02A': ['06', '$PSABASE/LibTrap_A016.dat'], '02B': ['07', '$PSABASE/LibTrap_B017.dat'], '02C': ['08', '$PSABASE/LibTrap_C013.dat'],
 '04A': ['12', '$PSABASE/LibTrap_A004.dat'], '04B': ['13', '$PSABASE/LibTrap_B004.dat'], '04C': ['14', '$PSABASE/LibTrap_C010.dat'],
 '05A': ['15', '$PSABASE/LibTrap_A001.dat'], '05B': ['16', '$PSABASE/LibTrap_B001.dat'], '05C': ['17', '$PSABASE/LibTrap_C006.dat'],
 '06A': ['18', '$PSABASE/LibTrap_A008.dat'], '06B': ['19', '$PSABASE/LibTrap_B009.dat'], '06C': ['20', '$PSABASE/LibTrap_C014.dat'],
 '07A': ['21', '$PSABASE/LibTrap_A014.dat'], '07B': ['22', '$PSABASE/LibTrap_B010.dat'], '07C': ['23', '$PSABASE/LibTrap_C016.dat'],
 '08A': ['24', '$PSABASE/LibTrap_A002.dat'], '08B': ['25', '$PSABASE/LibTrap_B007.dat'], '08C': ['26', '$PSABASE/LibTrap_C007.dat'],
 '09A': ['27', '$PSABASE/LibTrap_A017.dat'], '09B': ['28', '$PSABASE/LibTrap_B018.dat'], '09C': ['29', '$PSABASE/LibTrap_C018.dat'],
 '10A': ['30', '$PSABASE/LibTrap_A013.dat'], '10B': ['31', '$PSABASE/LibTrap_B015.dat'], '10C': ['32', '$PSABASE/LibTrap_C011.dat'],
 '11A': ['33', '$PSABASE/LibTrap_A010.dat'], '11B': ['34', '$PSABASE/LibTrap_B011.dat'], '11C': ['35', '$PSABASE/LibTrap_C009.dat'],
 '13A': ['39', '$PSABASE/LibTrap_A018.dat'], '13B': ['40', '$PSABASE/LibTrap_B012.dat'], '13C': ['41', '$PSABASE/LibTrap_C019.dat'],
 '14A': ['42', '$PSABASE/LibTrap_A009.dat'], '14B': ['43', '$PSABASE/LibTrap_B020.dat'], '14C': ['44', '$PSABASE/LibTrap_C005.dat'],
}
*/

const int NDets = 39; //39

string ADLpath = "ADL/";
int    DetId[NDets]   = {  0,  1,  2,
			   3,  4,  5,
			   6,  7,  8,
			  12, 13, 14,
			  15, 16, 17,
			  18, 19, 20,
			  21, 22, 23,
			  24, 25, 26,
			  27, 28, 29,
			  30, 31, 32,
			  33, 34, 35,
			  39, 40, 41,
			  42, 43, 44};

string ADLfile[NDets] = {"LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root",
			 "LibTrap_A001.root", "LibTrap_B001.root", "LibTrap_C001.root"};


// structure for PSA basis
struct pointPsa{
  float pos[3];              // pos in det frame
  int   netChSeg;
  float Amp[NCHAN][BSIZE];   // pulse shape for comparison
};


// structure for experimental traces
struct pointExp{
  int     numNetCharges;
  int     netChargeSegnum[NSEGS]; // normally given in the order of energy release
  float   netChargeEnergy[NSEGS];
  float   tAmp[NCHAN][BSIZE]; // the total original experimental data
  float   wAmp[NCHAN][BSIZE]; // the experimental data to be used and modified by the grid search ==> residuals
  float   sAmp[NCHAN][BSIZE]; // experimental trace manipulated for being used by the "chi2" loops
  float   rAmp[NCHAN][BSIZE]; // the "fitted" trace

  int     resPt[NSEGS];  // the best solution point for the numNetCharges segments
  float   resFac[NSEGS]; // relative amplitude of the point

  bool    isValid;
  bool    isInitialized;
  float   baseScale;          // temporarily used to manage the chi2 mapping
  int     netChSeg;
  char    localMask[NCHAN];
  int     bestPt;
  float   chi2min;                    // not very meaningful outside the grid-search

  // Produces sAmplitude from wAmplitude.
  void MakeSearchWave(char *mask){
    memset(sAmp, 0, sizeof(sAmp));
    for(int iSegm = 0; iSegm < NCHAN; iSegm++){
      if(mask[iSegm] != '0'){
        float *sA = sAmp[iSegm]; // the wave to be composed
        float *wA = wAmp[iSegm]; // taking the data from this
        for(int kk=0; kk<BSIZE; kk++){
          sA[kk] = wA[kk];
        }
      }
    }
  }

  // Modify wAmplitude (float) adding the base point passed in pPsa
  void AddBaseTrace(pointPsa &pPsa, const float fact) {
    int ptNetCharge = pPsa.netChSeg;
    for(int iSegm = 0; iSegm < NCHAN; iSegm++){ // segments and core
      float *realTrace = wAmp[iSegm]; // from the experimental data
      float *baseTrace;               // signal basis samples
      baseTrace = pPsa.Amp[iSegm];   // take them as they are

      for(int kk=0; kk<BSIZE; kk++) {
        (*realTrace++) += (*baseTrace++)*fact;
      }
    }

  }

};



Double_t GetTotalSystemMemory(){
  Long64_t pages = sysconf(_SC_PHYS_PAGES);
  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) pages)*page_size;
}

Double_t GetCurrentMemoryUsage(){
  Int_t tSize, resident, share;
  ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) resident) * page_size;
}



#endif
