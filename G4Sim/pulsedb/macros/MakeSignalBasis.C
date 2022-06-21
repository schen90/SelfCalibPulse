// convert pulse shape basis: add cross talk and preAmp responses

#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

using namespace std;

const int NSLIC = 6;
const int NSECT = 6;
const int NSEGS = NSLIC*NSECT;
const int NCHAN = NSEGS + 1;
const int INDCC = NCHAN - 1;  // 36 index of CC
const int TSTEP = 10; // time step of the experimental data and of the internal basis signals (ns)
const int BSIZE = 56;
const int basis_tzero = 10;
const int ZeroOffset = 8;

// Xtalk
float fDiffXTFactor = 10;
float xTalkProp[NSEGS][NSEGS];
float xTalkDiff[NSEGS][NSEGS];
float xTalkAver[NSEGS];
bool  segBroken[NSEGS+2]; // segBroken[NSEGS] always false; segBroken[NSEG+1] true if any broken segment is present
bool  xTalkCorr;

// preAmp
bool  fTauSliceGiven = true;
bool  fTauSegmentGiven = false;
float fTauSlice[7] = {35,35,35,35,35,35,35}; // preamp response of the 6 slices and the CC
float fTauSegment[NSEGS+1];
float tauSGCC[NSEGS+1]; //// the 36 segments and the CC (and some spare)
float fTauDecay = 0; // fall time in ns of the preamplifiers (default is not to correct basis for this)
float fSmooth = 0.1f; // final gaussian smoothing of the signal basis (0==no, default is 0.1)

struct pointPsa{
  int   ientry;
  float x, y, z;
  int   netChSeg;
  float Amp[NCHAN][BSIZE];

  void addXT(float *fwave1, float *factors1, float *fwave2, float *factors2){
    for(int sg = 0; sg < NCHAN; sg++){ // including CC
      float vv1 = factors1[sg];
      float vv2 = factors2[sg];
      if(sg==netChSeg){
	if(vv1){
	  for(int nn = 0; nn < BSIZE; nn++){
	    Amp[sg][nn] = Amp[sg][nn]*vv1;
	  }
	}
      }else{
	for(int nn = 0; nn < BSIZE; nn++){
	  Amp[sg][nn] += fwave1[nn]*vv1 + fwave2[nn]*vv2;
	}
      }
    }
  }


  // delta->exp, one value for each segment and the CC
  void convDeltaToExp(const float *Tau, float tDecayLong, float smooth) {
    for(int sg = 0; sg < NCHAN; sg++){ // including CC
      float *psg = Amp[sg];
      float tau1 = Tau[sg];
      if(tau1 > 0){
	float aa = (float)exp(-TSTEP/tau1);
	float gg = 1.f-aa;
	float vv = 0;
	for(int nn = 0; nn < BSIZE; nn++){
	  vv = gg*Amp[sg][nn] + aa*vv;
	  Amp[sg][nn] = vv;
	}
      }
      if(tDecayLong > 0){
	float tau2 = tDecayLong;
	float aa2 = (float)exp(-TSTEP/tau2);
	float vv1 = 0;
	float vv2 = 0;
	for(int nn = 0; nn < BSIZE; nn++){
	  vv2 = (Amp[sg][nn] - vv1) + aa2*vv2;
	  vv1 = Amp[sg][nn];
	  Amp[sg][nn] = vv2;
	}
      }
      if(smooth>0 && smooth<0.5f){
	// gaussian smoothing 1-8-1 --> att. ~60% at 50 MHz
	// assuming continuation of the edge values
	float s1 = smooth;
	float s2 = (1-2*smooth);
	float v2 = Amp[sg][0];
	float v1 = 0;
	float v0 = 0;
	psg = Amp[sg];
	for(int ii = 1; ii < BSIZE; ii++){
	  v0 = v1; v1 = v2; v2 = psg[ii];
	  psg[ii-1] = s1*v0 + s2*v1 + s1*v2;
	}
	v0 = v1; v1 = v2;
	psg[BSIZE-1] = s1*v0 + s2*v1 + s1*v2;
	{  // twice ??
	  float v2 = Amp[sg][0];
	  float v1 = 0;
	  float v0 = 0;
	  psg = Amp[sg];
	  for(int ii = 1; ii < BSIZE; ii++){
	    v0 = v1; v1 = v2; v2 = psg[ii];
	    psg[ii-1] = s1*v0 + s2*v1 + s1*v2;
	  }
	  v0 = v1; v1 = v2;
	  psg[BSIZE-1] = s1*v0 + s2*v1 + s1*v2;
	}
      }
    }
  }
  
};

vector<pointPsa> *segPts[NSEGS];

void ReadBasis(string basisfile);
void WriteBasis(string basisfile);

// Xtalk
bool ReadXtalkCoeffs(string cname);
bool CalcXtalkCoeffs(float *xTprop, float xfactor);
bool XtalkCorrection(int badSegment);

// preAmp responses
void ExponentialResponse(const float *TauSGCC, float TauDecay, float smooth);
void ExperimentalResponse(){;}
int PreampResponse(bool exponential);
  

void MakeSignalBasis(string basisfile){
  // read basis
  ReadBasis(Form("ori/%s",basisfile.c_str()));

  // write ori basis
  WriteBasis(Form("before/%s",basisfile.c_str()));
  
  //----------------------------------------------
  // Read and apply x-talk corrections
  //----------------------------------------------
  xTalkCorr = false;
  // read (direct) cross talk to modify the signal basis
  string fname = "param/xdir_1325-1340.cal";
  bool readOK = ReadXtalkCoeffs(fname);
  if(!readOK) return;

  // pass them to SignalBasis
  bool calcOK = CalcXtalkCoeffs(NULL, fDiffXTFactor);
  if(!calcOK) return;

  // apply the cross-talk correction to the basis signals
  XtalkCorrection(-1);

  // write basis with Xtalk
  WriteBasis(Form("Xtalk/%s",basisfile.c_str()));


  //----------------------------------------------
  // add preAmp response
  //----------------------------------------------
  // Response function (true==exponential, false=hardwired_experimental)
  PreampResponse(true);

  
  // Write basis
  WriteBasis(basisfile);

  return;
}


void ReadBasis(string basisfile){
  for(int sg = 0; sg < NSEGS; sg++){
    segPts[sg] = new vector<pointPsa>;
  }

  // read basis
  TFile *fbasis = new TFile(basisfile.c_str());
  if(!fbasis->IsOpen()){
    cerr<<"cannot find basisfile "<<basisfile<<endl;
    return;
  }
  cout<<"read pulse signal basis from "<<basisfile<<endl;

  const int nsig = 121; // npoints in 5ns step
  /*
  Int_t seg;
  Double_t pos[3];
  Double_t core[121];
  Double_t spulse[4356];
  TTree *basistree = (TTree *)fbasis->Get("tree");
  basistree->SetBranchAddress("seg",&seg);
  basistree->SetBranchAddress("pos",pos);
  basistree->SetBranchAddress("core",core);
  basistree->SetBranchAddress("spulse",spulse);
  */
  Int_t seg;
  Float_t pos[3];
  Float_t spulse[4477];
  TTree *basistree = (TTree *)fbasis->Get("tree");
  basistree->SetBranchAddress("seg",&seg);
  basistree->SetBranchAddress("pos",pos);
  basistree->SetBranchAddress("spulse",spulse);

  int npoint = basistree->GetEntriesFast();

  for(int ipoint=0; ipoint<npoint; ipoint++){
    if(ipoint%1000==0) cout<<"\r ipoint = "<<ipoint<<flush;
    basistree->GetEntry(ipoint);

    seg = seg-1; // input file start from 1

    pointPsa Pt;
    Pt.ientry = ipoint;
    Pt.netChSeg = seg;
    Pt.x = pos[0];  Pt.y = pos[1];  Pt.z = pos[2];

    for(int iseg=0; iseg<NSEGS; iseg++){
      for(int isig=0; isig<ZeroOffset; isig++) Pt.Amp[iseg][isig] = 0;
      for(int isig=0; isig<BSIZE-ZeroOffset; isig++){
        //Pt.Amp[iseg][isig] = spulse[iseg*nsig+isig];
        Pt.Amp[iseg][ZeroOffset+isig] = (spulse[iseg*nsig+basis_tzero+isig*2] + spulse[iseg*nsig+basis_tzero+isig*2+1])/2;
      }
    }

    for(int isig=0; isig<ZeroOffset; isig++) Pt.Amp[NCHAN-1][isig] = 0;
    for(int isig=0; isig<BSIZE-ZeroOffset; isig++){
      //Pt.Amp[NCHAN-1][isig] = core[isig];
      //Pt.Amp[NCHAN-1][isig] = (core[basis_tzero+isig*2] + core[basis_tzero+isig*2+1])/2;
      Pt.Amp[NCHAN-1][ZeroOffset+isig] = (spulse[NSEGS*nsig+basis_tzero+isig*2] + spulse[NSEGS*nsig+basis_tzero+isig*2+1])/2;
    }

    segPts[seg]->push_back(Pt);
  }

  cout<<"\r load "<<npoint<<" points from basis"<<endl;
  fbasis->Close();
  return;
}


void WriteBasis(string basisfile){
  // write basis
  cout<<"write pulse signal basis to "<<basisfile<<endl;
  TFile *fbasis = new TFile(basisfile.c_str(),"RECREATE");
  Int_t seg;
  Float_t pos[3];
  Float_t spulse[BSIZE*NCHAN];
  TTree *basistree = new TTree("tree","pulse shape db with XT and preAmp");
  basistree->Branch("seg",&seg);
  basistree->Branch("pos",pos,"pos[3]/F");
  basistree->Branch("spulse",spulse,Form("spulse[%d]/F",BSIZE*NCHAN));

  int ientry=0;
  int index[NSEGS]; for(int sg=0; sg<NSEGS; sg++) index[sg]=0;
  bool knext = true;
  while(knext){
    if(ientry%1000==0) cout<<"\r ientry = "<<ientry<<flush;
    knext = false;
    for(int sg=0; sg<NSEGS; sg++){
      if(index[sg]>=segPts[sg]->size()) continue;
      if(segPts[sg]->at(index[sg]).ientry!=ientry) continue;

      seg = sg;
      pos[0] = segPts[sg]->at(index[sg]).x;
      pos[1] = segPts[sg]->at(index[sg]).y;
      pos[2] = segPts[sg]->at(index[sg]).z;

      for(int iseg=0; iseg<NCHAN; iseg++)
	for(int isig=0; isig<BSIZE; isig++)
	  spulse[iseg*BSIZE+isig] = segPts[sg]->at(index[sg]).Amp[iseg][isig];

      /*
      for(int isig=0; isig<BSIZE; isig++)
	core[isig] = segPts[sg]->at(index[sg]).Amp[NCHAN-1][isig];
      */

      basistree->Fill();
      index[sg]++;
      ientry++;
      knext=true;
    }
  }
  
  cout<<"\r write "<<ientry<<" points to basis"<<endl;
  fbasis->cd();
  basistree->Write();
  fbasis->Close();
  return;
}

bool ReadXtalkCoeffs(string cname){
  int   nn, seggate, segseen;
  float xx;
  FILE *tfp;

  const float nonvalid = -1.e20f;

  for(segseen = 0; segseen < NSEGS; segseen++){
    for(seggate = 0; seggate < NSEGS; seggate++){
      xTalkProp[segseen][seggate] = nonvalid;
    }
  }

  if( (tfp = fopen(cname.c_str(), "r")) == NULL ){
    cout << "Error opening " << cname << endl;
    return false;
  }
  cout << "Reading cross-talk coefficients from " << cname << endl;
  nn = 0;
  while(fscanf(tfp, "%d %d %f", &segseen, &seggate, &xx) == 3){
    if(seggate < 0 || seggate >= NSEGS || segseen < 0 || segseen >= NSEGS){
      printf("\nWrong channel number: %d %d %f\n", seggate, segseen, xx);
      return false;
    }
    xTalkProp[segseen][seggate] = xx;

    nn++;
  }
  fclose(tfp);

  if(nn != NSEGS*NSEGS){
    cout << "Expecting " << NSEGS*NSEGS << " values;  found only " << nn << endl;
    return false;
  }

  // testing that all were really present
  int notgiven = 0;
  for(segseen = 0; segseen < NSEGS; segseen++){
    for(seggate = 0; seggate < NSEGS; seggate++){
      if(xTalkProp[segseen][seggate] <= nonvalid){
	printf("%3d %3d not given\n", segseen, seggate);
	notgiven++;
      }
    }
  }
  if(notgiven){
    cout << notgiven << " not given or repeated values" << endl;
    return false;
  }
  cout << "Found all " << NSEGS*NSEGS << " expected values" << endl;

  return true;
}


bool CalcXtalkCoeffs(float *xTprop, float xfactor){
  if(xTprop)
    memcpy(xTalkProp, xTprop, sizeof(xTalkProp));

  // should alse be possible to "break" a segment from the command line
  // This part WRONG because THE DIAGONAL ELEMENTS OF BROKEN SEGMENT ARE NO MORE ZERO !!!
  segBroken[NSEGS+1] = false;
  for(int nn = 0; nn < NSEGS; nn++){
    segBroken[nn] = xTalkProp[nn][nn] ? false : true;
    segBroken[NSEGS+1] |= segBroken[nn];
  }
  segBroken[NSEGS] = false; // core

  // second index refers to hit segment      (aggressor)
  // first  index refers to affected segment (victim)

  // average effect on a segment produced by "distant" aggressors
  // could recover a few more values by not fully rejecting the left and right sector
  memset(xTalkAver, 0, sizeof(xTalkAver));
  for(int n1 = 0; n1 < NSEGS; n1++){   // victim
    if(segBroken[n1]) continue;        // exclude broken segments
    float vv = 0;
    int   nn = 0;
    int   sector1 = n1/NSLIC;
    //int   slice1  = n1%NSLIC;
    for(int n2 = 0; n2 < NSEGS; n2++){ // aggressors
      if(segBroken[n2]) continue;
      int   sector2 = n2/NSLIC;
      //int   slice2  = n2%NSLIC;
      int secdist = (sector2>sector1) ? (sector2-sector1) : (sector2-sector1);
      if( secdist==0 || secdist==1 || secdist==(NSECT-1) )
	continue;
      vv += xTalkProp[n1][n2];
      nn++;
    }
    if(nn)
      xTalkAver[n1] = vv/nn;
  }

  // differential cross talk
  memset(xTalkDiff, 0, sizeof(xTalkDiff));

#define XTDTYPE 3
#if XTDTYPE == 0
  // average of non neighbours
  for(int n2 = 0; n2 < NSEGS; n2++){
    for(int n1 = 0; n1 < NSEGS; n1++){
      xTalkDiff[n1][n2] = xTalkAver[n1];
    }
  }
#elif XTDTYPE == 1
  // difference between average and actual
  for(int n2 = 0; n2 < NSEGS; n2++){
    for(int n1 = 0; n1 < NSEGS; n1++){
      xTalkDiff[n1][n2] = xTalkAver[n2]-xTalkProp[n1][n2];
    }
  }
#elif XTDTYPE == 2
  // same as proportional
  memcpy(xTalkDiff, xTalkProp, sizeof(xTalkDiff));
#elif XTDTYPE == 3
  // same as propotional, keeping only direct neighbours
  memcpy(xTalkDiff, xTalkProp, sizeof(xTalkDiff));
  for(int n1 = 0; n1 < NSEGS; n1++){   // victim
    int sector1 = n1/NSLIC;
    int slice1  = n1%NSLIC;
    for(int n2 = 0; n2 < NSEGS; n2++){ // aggressor
      int   sector2 = n2/NSLIC;
      int   slice2  = n2%NSLIC;
      float fact = 1.f;
      if(sector2==sector1){  // same column
	if(abs(slice2-slice1) != 1)
	  fact = 0;          // keep only +-=1
      }
      else if(abs(sector1-sector2)==1 ||
	      abs(sector1-sector2)==5){ // side column
	if(slice2 != slice1)
	  fact = 0;                     // keep only same slice
      }
      else{
	fact = 0;
      }
      xTalkDiff[n1][n2] *= fact;
    }
  }
#else
#endif

  // scaling as requested
  for(int n1 = 0; n1 < NSEGS; n1++){
    for(int n2 = 0; n2 < NSEGS; n2++){
      xTalkDiff[n1][n2] *= xfactor;
    }
  }

  // diagonals and broken (??) segments seg to 0
  for(int n1 = 0; n1 < NSEGS; n1++){
    for(int n2 = 0; n2 < NSEGS; n2++){
      if(n1==n2 || segBroken[n1] || segBroken[n2])
	xTalkDiff[n1][n2] = 0;
    }
  }

  // change sign of victims in the other triple-preamp of the aggressor sector
  // this does not really belong here: the matrix has to be provided externally!!
  for(int n2 = 0; n2 < NSEGS; n2 += 6){
    int sect2 = n2/NSLIC;
    int slic2 = n2%NSLIC;
    for(int n1 = 0; n1 < NSEGS; n1++){
      int sect1 = n1/NSLIC;
      int slic1 = n1%NSLIC;
      if(sect1 == sect2){
	if(slic1/3 != slic2/3)
	  xTalkDiff[n1][n2] *= -1.f;
	else{
	  int dist = (slic1==slic2) ? 1 : abs(slic2-slic1);
	  xTalkDiff[n1][n2] *= -10.f/dist;
	}
      }
    }
  }

  xTalkCorr = true;

  return true;
}


bool XtalkCorrection(int badSegment){
  if(!xTalkCorr)
    return false;

  float xvalue1[BSIZE] = {0}; // proportional xTalk signal
  float xvalue2[BSIZE] = {0}; // differential xTalk signal
  float xvalueD[BSIZE] = {0}; // reference signal for xvalue2
  float factor1[NCHAN] = {0}; // proportional xTalk weight
  float factor2[NCHAN] = {0}; // differential xTalk weight

  int count = 0;
  for(int ii=0; ii<NSEGS; ii++){ // only for the segments

    vector<pointPsa> *pPt = segPts[ii];
    int           npoints = pPt->size();
    if(npoints <= 0) continue;

    int          netChSeg = pPt->at(0).netChSeg;

    for(int nn = 0; nn < NSEGS; nn++){
      factor1[nn] = xTalkProp[nn][netChSeg]; // as we are reading the xdir matrix
      factor2[nn] = xTalkDiff[nn][netChSeg];
    }
    factor1[INDCC] = 0; // not for CC ??
    factor2[INDCC] = 0; // not for CC ??

    for(int jj=0; jj<npoints; jj++){
      // proportional xTalk signal, average of NC+CC
      for(int nn = 0; nn < BSIZE; nn++){
	xvalue1[nn] = (pPt->at(jj).Amp[netChSeg][nn] + pPt->at(jj).Amp[INDCC][nn])/2;
      }
      // reference signal for differential xTalk
      for(int nn = 0; nn < BSIZE; nn++){
	xvalueD[nn] = pPt->at(jj).Amp[netChSeg][nn];  // netCharge
      }
      // and its derivative as differential xTalk signal
      for(int nn = 0; nn < BSIZE-1; nn++){
	xvalue2[nn] = xvalueD[nn+1]-xvalueD[nn];
      }
      xvalue2[BSIZE-1] = 0;
      // make the corrections
      pPt->at(jj).addXT(xvalue1, factor1, xvalue2, factor2);
      count++;
    } // for(int jj=0; jj<npoints; jj++, pPt++)
    
  } // for(int ii=0; ii<NSEGS; ii++)

  return true;
}


// Filter all basis signals with an exponential response of tau.
// TauSGCC contains the decay constants for the 36 segments and the core
void ExponentialResponse(const float *TauSGCC, float TauDecay, float smooth){
  int count = 0;
  for(int ii=0; ii<NSEGS; ii++){
    int npoints = segPts[ii]->size();
    for(int jj=0; jj<npoints; jj++){
      segPts[ii]->at(jj).convDeltaToExp(TauSGCC, TauDecay, smooth);
      count++;
    }
  }
}



int PreampResponse(bool exponential){
  if(exponential){
    // convolution with an exponential response using tauDefault and/or values from PsaFilter.conf
    if(fTauSliceGiven){
      for(int sector = 0; sector < NSECT; sector++){
	for(int slice = 0; slice < NSLIC; slice++){
	  tauSGCC[sector*NSLIC+slice] = fTauSlice[slice];
	}
      }
      tauSGCC[NSEGS] = fTauSlice[NSLIC];
    }
    if(fTauSegmentGiven){
      for(int ns = 0; ns <= NSEGS; ns++){
	if(fTauSegment[ns] > 0)
	  tauSGCC[ns] = fTauSegment[ns];
      }
    }
    ExponentialResponse(tauSGCC, fTauDecay, fSmooth);  // ns
    cout << "  Decay constants (ns) of the exponential filter applied to the signal basis:";
    streamsize oldprec = cout.precision(1);
    cout << fixed;
    for(int slice = 0; slice < NSLIC; slice++){
      cout << "\n slice " << setw(1) << slice << " ";
      for(int sector = 0; sector < NSECT; sector++){
	cout << setw(8) << tauSGCC[sector*NSLIC+slice];
      }
    }
    cout << "\n core   " << setw(8) << tauSGCC[NSEGS] << endl;
    if(fSmooth > 0 && fSmooth < .5f)
      cout << "  Signals basis smoothed with an " << fSmooth << "-" << 1-2*fSmooth << "-" << fSmooth << " kernel" << endl;
    cout.precision(oldprec);
    cout.unsetf(ios_base::floatfield);
  }
  else{
    // alternatively, convolution with the experimental response derived from the pulser (kernel stored in fBasis)
    ExperimentalResponse();
    cout << "Signal basis modified using a hard-wired experimental response" << endl;
  }
  
  return 0;
}


#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc>1){
    MakeSignalBasis(string(argv[1]));
  }else{
    //MakeSignalBasis("pulseA.root");
    MakeSignalBasis("LibTrap_A001.root");
  }
  return 0;
}
#endif
