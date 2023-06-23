#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <TMath.h>
#include <TVector3.h>
#include "TRandom.h"
#include "TSystem.h"
#include <chrono>
#include <time.h>

#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "Global.hh"
#include "PSAFilter.cc"

using namespace std;

void AnaData(string tracefile, string psafile){

  // read tracefile
  struct OBJ{
    Int_t     EntryID;
    Float_t   SegTraces[DSIZE*NSEGS];
    Float_t   CoreTraces[DSIZE*2];
    Float_t   SegE[NSEGS];
    Float_t   CoreE[2];
    Float_t   CoreT[2];
    Int_t     CrystalId;
    ULong64_t CrystalTS;
  };
  OBJ obj;

  TFile *fin = new TFile(tracefile.c_str());
  if(!fin->IsOpen()){
    cerr<<"cannot find tracefile "<<tracefile<<endl;
    return;
  }
  TTree *intree = (TTree *)fin->Get("TreeMaster");
  intree->SetBranchAddress("EntryID",   &obj.EntryID);
  intree->SetBranchAddress("SegTraces",  obj.SegTraces);
  intree->SetBranchAddress("CoreTraces", obj.CoreTraces);
  intree->SetBranchAddress("SegE",       obj.SegE);
  intree->SetBranchAddress("CoreE",      obj.CoreE);
  intree->SetBranchAddress("CoreT",      obj.CoreT);
  intree->SetBranchAddress("CrystalId", &obj.CrystalId);
  intree->SetBranchAddress("CrystalTS", &obj.CrystalTS);
  int nentries = intree->GetEntriesFast();
  cout<<"read "<<nentries<<" events from "<<tracefile<<endl;


  // output psafile
  TFile *fout = new TFile(psafile.c_str(),"RECREATE");
  TTree *anatree = new TTree("tree","PSA tree");
  Int_t   numNetCharges;
  Int_t   seg;
  Float_t NetCh;
  Float_t pos[3];
  Float_t segpos[3];
  Float_t chi2;
  anatree->Branch("EntryID",   &obj.EntryID,   "EntryID/I");
  anatree->Branch("SegE",       obj.SegE,      "SegE[36]/F");
  anatree->Branch("CoreE",      obj.CoreE,     "CoreE[2]/F");
  anatree->Branch("CrystalId", &obj.CrystalId, "CrystalId/I");

  anatree->Branch("numNetCharges", &numNetCharges, "numNetCharges/I");
  anatree->Branch("seg",   &seg,   "seg/I");
  anatree->Branch("NetCh", &NetCh, "NetCh/F");
  anatree->Branch("pos",    pos,   "pos[3]/F");
  anatree->Branch("segpos", segpos,"segpos[3]/F");

  anatree->Branch("chi2",  &chi2,  "chi2/F");


  // PSA
  PSAFilter *fpsa[NDets];
  for(int idet=0; idet<NDets; idet++){
    cout<<"PSA for CrystalId "<<DetId[idet]<<" : ";
    fpsa[idet] = new PSAFilter(ADLpath+ADLfile[idet]);
  }

  // clock
  time_t start, stop;
  time(&start);

  // loop all entries
  int ientry;
  for(ientry=0; ientry<nentries; ientry++){
    if(ientry%1000==0)
      cout<<"\r finish "<<ientry<<" / "<<nentries<<" entries..."<<flush;

    intree->GetEntry(ientry);

    // check the basis loaded
    int idet = -1;
    for(int id=0; id<NDets; id++){
      if(DetId[id]==obj.CrystalId){
	idet = id;
	break;
      }
    }
    if(idet<0) continue;

    /////////////////////////////////////////
    // read data
    /////////////////////////////////////////
    pointExp pE;
    // read SegTrace
    for(int iseg=0; iseg<NSEGS; iseg++){
      float last = 0;
      int ann = 0;
      int bnn = DZERO;
      for(; ann<BSIZE; ann++, bnn++){
	if(bnn < DSIZE){
	  last = obj.SegTraces[DSIZE*iseg+bnn];
	}
	pE.tAmp[iseg][ann] = last;
      }
    }
    // read CoreTrace
    {
      float last = 0;
      int ann = 0;
      int bnn = DZERO;
      for(; ann<BSIZE; ann++, bnn++){
	if(bnn < DSIZE){
	  last = obj.CoreTraces[bnn];
	}
	pE.tAmp[INDCC][ann] = last;
      }
    }

    // find NetCharge segments
    int numsegs = 0;
    for(int iseg=0; iseg<NSEGS; iseg++){
      if(obj.SegE[iseg]>0){
	pE.netChargeSegnum[numsegs] = iseg;
	pE.netChargeEnergy[numsegs] = obj.SegE[iseg];
	numsegs++;
      }
    }
    pE.numNetCharges = numsegs;
    pE.netChSeg = -1;
    pE.bestPt = -1;
    pE.chi2min = float(1.e30);
    pE.isValid = false;
    pE.isInitialized = false;

    ////////////////////////////////////
    // grid search
    ////////////////////////////////////
    fpsa[idet]->ProcessOneEvent(pE);


    ////////////////////////////////////
    // output results
    ////////////////////////////////////
    numNetCharges = pE.numNetCharges;
    chi2 = fpsa[idet]->GetChi2(pE);
    for(int snum=0; snum<numNetCharges; snum++){
      seg   = pE.netChargeSegnum[snum];
      NetCh = pE.netChargeEnergy[snum];

      int bestPt = pE.resPt[snum];
      fpsa[idet]->GetPtPos(seg, bestPt, pos);
      fpsa[idet]->GetPtPos(seg, -1, segpos); //segment center

      anatree->Fill();
    }


  }// end of loop entries
  cout<<"\r finish "<<ientry<<" / "<<nentries<<" entries.  "<<endl;

  time(&stop);
  printf("============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  cout<<"write to "<<psafile<<" ..."<<endl;
  fout->cd();
  anatree->Write();
  fout->Close();
  fin->Close();

  return;
}


#ifndef __CINT__
int main(int argc, char *argv[]){
  if(argc!=3){
    cout << "Usage: \n"
	 << argv[0]
	 << " trace-file psa-file\n";
    return 0;
  }

  AnaData(string(argv[1]), string(argv[2]));

  return 0;
}
#endif
