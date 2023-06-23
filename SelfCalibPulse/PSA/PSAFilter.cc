#include <stdlib.h>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <TMath.h>
#include <TVector3.h>
#include "TRandom.h"
#include "TSystem.h"
#include <chrono>
#include <x86intrin.h>

#include "PSAFilter.hh"

using namespace std;

PSAFilter::PSAFilter(string basisfile){
  MakeSegmentMap(2);

  // read signal basis
  ReadPSAbasis(basisfile);
}


void PSAFilter::MakeSegmentMap(int neighbours){
  cout<<"Make Segment Map...";
  memset(hmask, '0', sizeof(hmask));
  // '1' self, '2' neighbour, '9' CC
  for(int iseg=0; iseg<NSEGS; iseg++){
    int isec = iseg/NSLIC;
    int isli = iseg%NSLIC;

    for(int jseg=0; jseg<NSEGS; jseg++){
      int jsec = jseg/NSLIC;
      int jsli = jseg%NSLIC;

      int distV = abs(jsli - isli);
      int distH = abs(jsec - isec);
      distH = min(distH, abs(jsec - isec + NSECT));
      distH = min(distH, abs(jsec - isec - NSECT));
      int mdist = distV + distH;
      if(mdist<=neighbours && distH<neighbours && distV<neighbours){
	hmask[iseg][jseg] = (iseg==jseg) ? '1' : '2';
      }
    }

    hmask[iseg][NCHAN-1] = '9';
    hmask[iseg][NCHAN] = 0; // to close each line as a string
  }
  cout<<endl;
  return;
}


void PSAFilter::ReadPSAbasis(string basisfile){

  double MemTotalGB = GetTotalSystemMemory()/GB;

  TFile *fdb = new TFile(basisfile.c_str());
  if(!fdb->IsOpen()){
    cerr<<"cannot find ADLfile "<<basisfile.c_str()<<endl;
  }

  TTree *dbtree = (TTree *)fdb->Get("tree");

  Int_t   dbsegi;
  Float_t dbposi[3];
  Float_t dbspulsei[BSIZE*NCHAN];

  dbtree->SetBranchAddress(    "seg",   &dbsegi);
  dbtree->SetBranchAddress(    "pos",    dbposi);
  dbtree->SetBranchAddress( "spulse", dbspulsei);
  int npoint = dbtree->GetEntriesFast();

  for(int iseg=0; iseg<NSEGS; iseg++)
    averPt[iseg].pos[0] = averPt[iseg].pos[1] = averPt[iseg].pos[2] = 0;

  for(int ipoint=0; ipoint<npoint; ipoint++){
    dbtree->GetEntry(ipoint);

    pointPsa Pt;
    for(int ix=0; ix<3; ix++){
      Pt.pos[ix] = dbposi[ix];
      averPt[dbsegi].pos[ix] += dbposi[ix];
    }
    Pt.netChSeg = dbsegi;

    for(int iseg=0; iseg<NCHAN; iseg++){
      for(int isig=0; isig<BSIZE; isig++){
	Pt.Amp[iseg][isig] = dbspulsei[iseg*BSIZE+isig];
      }
    }

    segPts[dbsegi].push_back(Pt);
  }

  for(int iseg=0; iseg<NSEGS; iseg++){
    for(int ix=0; ix<3; ix++){
      averPt[iseg].pos[ix] = averPt[iseg].pos[ix] / segPts[iseg].size();
    }
    float cdist = fstep*0.49f; // slightly less than half fine step
    averPtlist[iseg] = NearestPoint(averPt[iseg].pos[0],averPt[iseg].pos[1],averPt[iseg].pos[2],iseg,cdist);
  }

  double MemUsageGB = GetCurrentMemoryUsage()/GB;
  cout<<"load "<<npoint<<" points from "<<basisfile
      <<Form(" Mem: %.1f/%.1f",MemUsageGB,MemTotalGB)<<endl;

  fdb->Close();

  return;
}


inline int PSAFilter::NearestPoint(float px, float py, float pz, int iseg, float& dist){
  int npts = segPts[iseg].size();

  int cnear = -1;
  float cprec = dist*dist;

  float cdist = (float)1.e20;
  for(int jj = 0; jj < npts; jj++){
    float dx = segPts[iseg][jj].pos[0] - px;
    float dy = segPts[iseg][jj].pos[1] - py;
    float dz = segPts[iseg][jj].pos[2] - pz;
    float d2 = dx*dx + dy*dy + dz*dz;
    if(d2 < cdist){
      cdist = d2;
      cnear = jj;
      if(cdist < cprec)
        break;
    }
  }
  dist = (float)sqrt(cdist);
  return cnear;
}



void PSAFilter::ProcessOneEvent(pointExp &PE){

  for(int rpt=0; rpt<4; rpt++){
    //////////////////////////////////////
    // 1  prepare event for this iteration
    //////////////////////////////////////

    PrepareEvent(PE, (rpt==0));

    if(!PE.isValid){
      return;
    }
    memset(PE.rAmp, 0, sizeof(PE.rAmp)); // Clear the "fitted" trace
    InitialSolution(PE, (rpt!=0) );  // first time: center of the fired segment; otherwise: position of previous iteration

    /////////////////////////////////////
    // 2  call the 1-hit grid-search algo
    /////////////////////////////////////

    int rv = ProcessTheEvent(PE);

  }

}


void PSAFilter::PrepareEvent(pointExp &PE, bool doReset){
  if(doReset){
    // find hit segments
    int    numsegs  = PE.numNetCharges;
    int   *segOrder = PE.netChargeSegnum;
    float *eneOrder = PE.netChargeEnergy;

    if(numsegs<1 || numsegs>12){
      PE.isValid = false;
      return;
    }

    PE.netChSeg = -1;
    PE.isValid = true;
    if(numsegs>1){
      // order the net-charge segments according to the released energy (not done for GridSearch==SegCenter)
      for(int n1=numsegs-1; n1>=0; n1--){
        for(int n2=0; n2<n1; n2++){
          if(eneOrder[n2]<eneOrder[n2+1]){
            swap(segOrder[n2], segOrder[n2+1]);
            swap(eneOrder[n2], eneOrder[n2+1]);
          }
        }
      }
    }

  }else{
    PE.netChSeg = -1;
    PE.isValid = true;
  }

  // save tAmplitude into wAmplitude, where is can be dynamically modified by the search algorithms
  memcpy(PE.wAmp, PE.tAmp, sizeof(PE.wAmp));
}


// place everything at the center or keep what already found in the previous minimization
void PSAFilter::InitialSolution(pointExp &PE, bool keepPrevious){
  for(int snum = 0; snum < PE.numNetCharges; snum++){
    int   netChSeg  = PE.netChargeSegnum[snum];
    float netChEner = PE.netChargeEnergy[snum];
    float scaleFact = netChEner;

    int bestPt = (keepPrevious && PE.isInitialized) ? PE.resPt[snum] : averPtlist[netChSeg];

    PE.resPt[snum]  = bestPt;
    PE.resFac[snum] = 1.f;

    PE.netChSeg = netChSeg;   // now working on this

    // subtract the result of this point from PF.wAmplitude
    pointPsa &bestPoint = segPts[netChSeg][bestPt];
    PE.AddBaseTrace(bestPoint, -scaleFact);
  }

  PE.isInitialized = true;
}


int PSAFilter::ProcessTheEvent(pointExp &PE){

  int chiDone = 0;

  if(!PE.isValid) return chiDone;

  for(int snum = 0; snum < PE.numNetCharges; snum++){
    int   netChSeg  = PE.netChargeSegnum[snum];
    float netChEner = PE.netChargeEnergy[snum];

    float scaleFact = netChEner;
    PE.baseScale   = scaleFact;
    PE.netChSeg    = netChSeg; // now working on this

    // add-back the initial/previous solution to PE.wAmplitude
    int bestPt = PE.resPt[snum];
    pointPsa &bestPoint = segPts[netChSeg][bestPt];
    PE.AddBaseTrace(bestPoint, scaleFact);

    // prepare sAmplitude from the subset of active segments
    int nActive = MakeSearchWave(PE);

    PE.chi2min = float(1.e30);
    PE.bestPt = -1;

    int rv = SearchFullGrid(PE, netChSeg, netChEner);
    chiDone += abs(rv);


    pointPsa &bestPoint1 = segPts[netChSeg][PE.bestPt];

    // The best point is normalized to the experimental amplitude and accumulated in rAmplitude
    AddToSolution(PE, bestPoint1, scaleFact);

    // Subtract the result of this point from PF.wAmplitude, which is the the spectrum of residuals.
    PE.AddBaseTrace(bestPoint1, -scaleFact);
    PE.resPt[snum] = PE.bestPt;
    PE.resFac[snum] = 1.f;

  }  // loop over the net-charge segments

  return chiDone;
}


int PSAFilter::SearchFullGrid(pointExp &PE, int netChSeg, float netChEner){
  char *lMask = PE.localMask;

  int   bestPt  = 0;
  float chi2min = PE.chi2min;
  int   chiDone = 0;

  int iPts = 0;
  int nPts = segPts[netChSeg].size();

  for(; iPts < nPts; iPts++){
    chiDone++;
    float baseScale = PE.baseScale;
    float chi2 = 0;
    for(int iseg=0; iseg<NCHAN; iseg++){
      if(lMask[iseg] != '0'){
        float realTrace[BSIZE];
        memcpy(realTrace, PE.sAmp[iseg], sizeof(realTrace));
        float *baseTrace = segPts[netChSeg][iPts].Amp[iseg];
        chi2 += Chi2InnerLoop(realTrace, baseTrace, baseScale);
        if(chi2 > chi2min) break;
      }
    } // end loop over the segments
    if(chi2 < chi2min){
      bestPt  = iPts;
      chi2min = chi2;
    }
  }// end loop over the base points iPts

  PE.bestPt  = bestPt;
  PE.chi2min = chi2min;

  return chiDone;
}


// Add the properly scaled and shifted bestPoint to PF.rAmp
void PSAFilter::AddToSolution(pointExp &PE, pointPsa &bestPoint, float scaleFact){
  float * realTrace;
  float * baseTrace;
  for(int iseg=0; iseg<NCHAN; iseg++){
    realTrace = PE.rAmp[iseg];
    baseTrace = bestPoint.Amp[iseg];
    for(int nn=0; nn<BSIZE; nn++){
      (*realTrace++) += (*baseTrace++)*scaleFact;
    }
  }
}


float PSAFilter::GetChi2(pointExp &PE){
  float chi2 = 0;
  for(int iseg=0; iseg<NCHAN; iseg++){
    float realTrace[BSIZE];
    memcpy(realTrace, PE.tAmp[iseg], sizeof(realTrace));
    float baseTrace[BSIZE];
    memcpy(baseTrace, PE.rAmp[iseg], sizeof(baseTrace));
    chi2 += Chi2InnerLoop(realTrace, baseTrace, 1);
  }
  return chi2;
}


int PSAFilter::MakeSearchWave(pointExp &PE){

  // select the active segments
  int nActive = MakeLocalMask(PE);

  // calculate sAmplitude
  PE.MakeSearchWave(PE.localMask);

  return nActive;
}


int PSAFilter::MakeLocalMask(pointExp &PE){
  // set mask of segments for the search loop
  char *localMask = PE.localMask;
  strcpy(localMask,hmask[PE.netChSeg]);
  int sMult = PE.numNetCharges;
  // removing common segments with the other net charges
  for(int ii=0; ii<sMult; ii++){
    int ncseg = PE.netChargeSegnum[ii];
    if(ncseg!=PE.netChSeg){
      // net-charge segments are always removed
      localMask[ncseg] = '0';
    }
  }

  int ns = 0;
  for(int ii=0; ii<NCHAN; ii++){
    if(localMask[ii]!='0')
      ns++;
  }
  return ns;
}


inline float PSAFilter::Chi2InnerLoop(const float *apulse, const float *bpulse, float bScale){

  const __m128 masks = _mm_set1_ps(-0.0f);
  const __m128 zeros = _mm_setzero_ps();
  const __m128 misig = _mm_set1_ps(0.01); // min sigma

  __m128 loopFac    = _mm_set_ps(bScale, bScale, bScale, bScale);

  __m128* realtrace = (__m128*)apulse;
  __m128* basetrace = (__m128*)bpulse;

  __m128 wave = _mm_setzero_ps();
  __m128 diff = _mm_setzero_ps();
  __m128 sigm = _mm_setzero_ps();
  __m128 chis = _mm_setzero_ps();

  for(int nn=0; nn<LOOP_4SSE; nn++){
    wave = _mm_mul_ps(loopFac, basetrace[nn]);
    diff = _mm_sub_ps(realtrace[nn], wave);
#if   CHI2 == CHI2_SQ
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff)); // square
#elif CHI2 == CHI2_CHI2
    //sigm = _mm_andnot_ps(masks, wave);
    sigm = _mm_max_ps(_mm_andnot_ps(masks, wave), misig);    // chi2
    chis = _mm_add_ps(chis, _mm_div_ps(_mm_mul_ps(diff,diff), sigm)); // chi2
#elif CHI2 == CHI2_CHI2_2
    sigm = _mm_max_ps(_mm_andnot_ps(masks, wave), sigm); // chi2_2
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff));               // chi2_2
#elif CHI2 == CHI2_ABS
    chis = _mm_add_ps(chis, _mm_andnot_ps(masks, diff)); // ABS value
#elif CHI2 == CHI2_FABS
    sigm = _mm_max_ps(_mm_andnot_ps(masks, wave), misig);         // fABS
    chis = _mm_add_ps(chis, _mm_div_ps(_mm_andnot_ps(masks, diff), sigm)); // fABS
#elif CHI2 == CHI2_SQRT
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_andnot_ps(masks, diff))); // sqrt
#elif CHI2 == CHI2_2SQRT
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_sqrt_ps(_mm_andnot_ps(masks, diff)))); // 2sqrt
#else
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff)); // square
#endif
  }

#if   CHI2 == CHI2_CHI2_2
  sigm = _mm_max_ps(sigm, misig); // chi2_2
  chis = _mm_div_ps(chis, sigm);  // chi2_2
#endif

  chis = _mm_hadd_ps(_mm_hadd_ps(chis, zeros), zeros);

  float chi2 = _mm_cvtss_f32(chis);

  return chi2;
}

