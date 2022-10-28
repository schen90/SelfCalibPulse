#ifndef AGATA_CC
#define AGATA_CC

#include <numeric>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <x86intrin.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TMath.h>
#include "TVector3.h"
#include "TMatrixD.h"
#include <vector>

#include "AGATA.hh"
#include "AGATAgeo.hh"

using namespace std;

AGATA::AGATA(){

  agatageo = new AGATAgeo();
  NDets = agatageo->GetNDets();

  for(int idet=0; idet<NDets; idet++){
    for(int iseg=0; iseg<NSeg; iseg++){
      fPSC0[idet][iseg] = new vector<PSC*>();
      fPSC1[idet][iseg] = new vector<PSC*>();
    }
  }
  
  cPSCmem  = 0;
  cPSCfile = 0;
  maxnhits = 0;
  
  //InitPSC();

}

AGATA::~AGATA(){
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// Pulse Shape Collections
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
void AGATA::InitPSC(){
  cout<<"\e[1;31m Init PSC points: \e[0m"<<endl;

  // Get positions for signal basis
  string dbfile[3] = {"G4Sim/pulsedb/LibTrap_A001.root",
                      "G4Sim/pulsedb/LibTrap_B001.root",
                      "G4Sim/pulsedb/LibTrap_C001.root"};

  vector<Int_t>    dbseg[NType];
  vector<TMatrixD> dbpos[NType];

  for(int itype=0; itype<NType; itype++){

    TFile *fdb = new TFile(dbfile[itype].c_str());
    if(!fdb->IsOpen()){
      cerr<<"cannot find dbfile "<<dbfile[itype]<<endl;
      return;
    }

    TTree *dbtree = (TTree *)fdb->Get("tree");
    Int_t dbsegi;
    Float_t dbposi[3];
    dbtree->SetBranchAddress("seg",&dbsegi);
    dbtree->SetBranchAddress("pos",dbposi);
    int npoint = dbtree->GetEntriesFast();

    TMatrixD tmppos(3,1);
    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);

      for(int i=0; i<3; i++) tmppos(i,0)=dbposi[i];

      dbseg[itype].push_back(dbsegi); // start from 0
      dbpos[itype].push_back(tmppos);
    }

    cout<<"load "<<npoint<<" points from "<<dbfile[itype]<<endl;
    fdb->Close();
  }

  // Greate PSC for every detector
  int idet, itype, npoint;
  float MemUsageGB, MemTotalGB;
  cPSCmem = 0;
  for(idet=0; idet<NDets; idet++){
    itype = idet%3;
    npoint = dbseg[itype].size();

    for(int ipoint=0; ipoint<npoint; ipoint++){
      int iseg = dbseg[itype][ipoint];
      TMatrixD DetPos = dbpos[itype][ipoint];
      TMatrixD LabPos = agatageo->Det2LabPos(idet,DetPos);

      PSC* apsc0 = new PSC(idet, iseg);
      PSC* apsc1 = new PSC(idet, iseg);
      for(int ix=0; ix<3; ix++){
	apsc0->labpos[ix] = apsc1->labpos[ix] = LabPos(ix,0);
	apsc0->detpos[ix] = apsc1->detpos[ix] = DetPos(ix,0);
      }
      fPSC0[idet][iseg]->push_back(apsc0);
      fPSC1[idet][iseg]->push_back(apsc1);
      PSCpos[idet][iseg].push_back( TVector3( LabPos(0,0), LabPos(1,0), LabPos(2,0) ) );
      
      cPSCmem++;
      if(cPSCmem%10000==0){
	MemUsageGB = GetCurrentMemoryUsage()/GB;
	MemTotalGB = GetTotalSystemMemory()/GB;

	cout<<"\r Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	    <<" PSC-"<<cPSCmem<<" x2.."
	    <<" create PSC for det "<<idet<<": "<<npoint<<" points..."<<flush;
      }

    }
  }
  MemUsageGB = GetCurrentMemoryUsage()/GB;
  MemTotalGB = GetTotalSystemMemory()/GB;

  cout<<"\r Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
      <<" PSC-"<<cPSCmem<<" x2.."
      <<" create PSC for det "<<idet<<": "<<npoint<<" points..."<<endl;

  iPSC0 = 0;
  
  return;
}


// Copy fPSC1 to fPSC0
void AGATA::CopyPSC(){
  cout<<"\e[1;33m Update Signal Basis for PSA \e[0m"<<endl;
  for(int idet=0; idet<NDets; idet++){
    for(int iseg=0; iseg<NSeg; iseg++){
      int nc = fPSC0[idet][iseg]->size();
      for(int ic=0; ic<nc; ic++){
	fPSC0[idet][iseg]->at(ic)->nhits = fPSC1[idet][iseg]->at(ic)->nhits;

	copy_n(fPSC1[idet][iseg]->at(ic)->avelpos, 3, fPSC0[idet][iseg]->at(ic)->avelpos);
	copy_n(fPSC1[idet][iseg]->at(ic)->avedpos, 3, fPSC0[idet][iseg]->at(ic)->avedpos);

	for(int iiseg=0; iiseg<NSegCore; iiseg++){
          copy_n(fPSC1[idet][iseg]->at(ic)->spulse[iiseg], NSig, fPSC0[idet][iseg]->at(ic)->spulse[iiseg]);
        }
	
      }
      
    }
  }
  iPSC0++;
  return;
}


// Clear fPSC
void AGATA::ClearPSC1(){
  cout<<"\e[1;33m Clear fPSC1 \e[0m"<<endl;
  float zero3[3] = {0,0,0};
  float zerosig[NSig]; for(int isig=0; isig<NSig; isig++) zerosig[isig]=0;
  
  for(int idet=0; idet<NDets; idet++){
    for(int iseg=0; iseg<NSeg; iseg++){
      for(PSC *apsc : *fPSC1[idet][iseg]){
	apsc->nhits = 0;
	
	copy_n(zero3, 3, apsc->avelpos);
	copy_n(zero3, 3, apsc->avedpos);

	for(int iiseg=0; iiseg<NSegCore; iiseg++){
          copy_n(zerosig, NSig, apsc->spulse[iiseg]);
        }

      }
    }
  }

  maxnhits = 0;
  return;
}


// Read Pulse Shape Collection to fPSC0 for PSA
void AGATA::ReadPSCfiles(int detid){
  vector <int> idlist; // detid list to write
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  string pscfilesname = "PSCfiles/Det";

  cout<<"\e[1;33m Read Pulse Shape Collections from "<<pscfilesname;
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;

  int idet, iseg, nentries;
  float MemUsageGB, MemTotalGB;
  cPSCmem = 0;
  for(int i=0; i<idlist.size(); i++){
    idet = idlist[i];
    pscfile[idet] = new TFile(Form("%s%04d.root",pscfilesname.c_str(),idet));

    for(iseg=0; iseg<NSeg; iseg++){
      psctree[idet][iseg] = (TTree*)pscfile[idet]->Get(Form("tree%d",iseg));

      InitTreeRead(psctree[idet][iseg]);

      nentries = psctree[idet][iseg]->GetEntriesFast();

      for(int i=0; i<nentries; i++){

	psctree[idet][iseg]->GetEntry(i);

	int idx = FindPSC(det, seg, detpos, i);
	if(idx<0) continue;

	PSC *apsc = fPSC0[det][seg]->at(idx);
	apsc->nhits = nhits;

	copy_n(avelpos, 3, apsc->avelpos);
	copy_n(avedpos, 3, apsc->avedpos);

	for(int iiseg=0; iiseg<NSegCore; iiseg++){
          copy_n(spulse[iiseg], NSig, apsc->spulse[iiseg]);
        }

	cPSCmem++;
	if(cPSCmem%10000==0){
          MemUsageGB = GetCurrentMemoryUsage()/GB;
          MemTotalGB = GetTotalSystemMemory()/GB;

	  cout<<"\r Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
	      <<" PSC-"<<cPSCmem<<".."
	      <<" read PSC for det "<<idet<<": "<<nentries<<" entries..."<<flush;
        }

      }
    }
    
  }
  MemUsageGB = GetCurrentMemoryUsage()/GB;
  MemTotalGB = GetTotalSystemMemory()/GB;

  cout<<"\r Mem "<<Form("%.1f/%.1f",MemUsageGB,MemTotalGB)<<"GB.."
      <<" PSC-"<<cPSCmem<<".."
      <<" read PSC for det "<<idet<<": "<<nentries<<" entries..."<<endl;

  iPSC0 = 1;
  return;
}


// write Pulse Shape Collection fPSC1 to files
void AGATA::WritePSCfiles(int detid){

  vector <int> idlist; // detid list to write
  for(int idet=0; idet<NDets; idet++){
    if(detid>-1 && idet!=detid) continue;
    idlist.push_back(idet);
  }
  string pscfilesname0 = "./PSCfiles/tmp/Det";
  string pscfilesname = "./PSCfiles/Det";

  cout<<"\e[1;33m Create PSCfiles for "<<pscfilesname;
  cout<<idlist[0];  if(idlist.size()>1) cout<<" ~ "<<idlist.back();
  cout<<" ... \e[0m"<<endl;

  cPSCfile = 0;
  for(int idet : idlist){

    pscfile[idet] = new TFile(Form("%s%04d.root",pscfilesname0.c_str(),idet),"RECREATE");

    for(int iseg=0; iseg<NSeg; iseg++){
      cout<<"\r writing tree for Det "<<idet<<" Seg "<<iseg<<" ..."<<flush;
      psctree[idet][iseg] = new TTree(Form("tree%d",iseg),Form("pulse shape tree for det%04d seg%d",idet,iseg));

      InitTreeWrite(psctree[idet][iseg]);

      for(PSC *apsc : *fPSC1[idet][iseg]){

	det = apsc->det;
	seg = apsc->seg;
	nhits = apsc->nhits;

	copy_n(apsc->labpos, 3, labpos);
	copy_n(apsc->detpos, 3, detpos);

	copy_n(apsc->avelpos, 3, avelpos);
	copy_n(apsc->avedpos, 3, avedpos);

	dist = 0;
	for(int ix=0; ix<3; ix++) dist += SQ(avelpos[ix]-labpos[ix]);
        dist = sqrt(dist);

	for(int iiseg=0; iiseg<NSegCore; iiseg++){
          copy_n(apsc->spulse[iiseg], NSig, spulse[iiseg]);
        }
	
	psctree[idet][iseg]->Fill();
	cPSCfile++;
      }
      pscfile[idet]->cd();
      psctree[idet][iseg]->Write();
    }// end of loop seg

    pscfile[idet]->Close();
    gROOT->ProcessLine(Form(".!mv -f %s%04d.root %s%04d.root", pscfilesname0.c_str(),idet, pscfilesname.c_str(),idet));
  }

  cout<<"\e[1;33m Write "<<cPSCfile<<" PSC to files \e[0m"<<endl;
  
  return;
}


void AGATA::GetPSCstat(int *PSCstat){
  PSCstat[0] = cPSCmem;
  PSCstat[1] = cPSCfile;
  PSCstat[2] = maxnhits;
  return;
}


int AGATA::FindPSC(int detid, int segid, float dpos[], int istart){
  int idx = -1;
  int ndx = fPSC0[detid][segid]->size();
  
  for(int i=0; i<ndx; i++){
    int itmp = (istart+i)%ndx;
    PSC* apsc = fPSC0[detid][segid]->at(itmp);

    if(fabs(dpos[0] - apsc->detpos[0])>0.2) continue;
    if(fabs(dpos[1] - apsc->detpos[1])>0.2) continue;
    if(fabs(dpos[2] - apsc->detpos[2])>0.2) continue;

    idx = itmp;
    break;
  }

  return idx;
}


int AGATA::AddPStoPSC(PS *aps, int ipsc){
  int detid = aps->det;
  int segid = aps->seg;

#ifdef NTHREADS
  lock_guard<mutex> lock(PSCmtx[detid][segid]); // lock PSC for one segment
#endif
  
  if(ipsc > fPSC1[detid][segid]->size()-1){
    cout<<Form("cannot find ipsc = %d ; fPSC1[%d][%d].size() = %d",ipsc,detid,segid,(int)fPSC1[detid][segid]->size())<<endl;
    return -1;
  }

  double nhitstmp = (double)fPSC1[detid][segid]->at(ipsc)->nhits;

  // calc average position
  TMatrixD pos(3,1); for(int ix=0; ix<3; ix++) pos(ix,0) = aps->labpos[ix];
  TMatrixD LabPos(3,1);
  for(int ix=0; ix<3; ix++) LabPos(ix,0) = fPSC1[detid][segid]->at(ipsc)->avelpos[ix];
  LabPos = 1./(nhitstmp+1)*(nhitstmp*LabPos + pos);
  TMatrixD DetPos = agatageo->Lab2DetPos(detid,LabPos);

  for(int ix=0; ix<3; ix++){
    fPSC1[detid][segid]->at(ipsc)->avelpos[ix] = LabPos(ix,0);
    fPSC1[detid][segid]->at(ipsc)->avedpos[ix] = DetPos(ix,0);
  }

  // calc average pulse
  for(int iseg=0; iseg<NSegCore; iseg++){
    for(int isig=0; isig<NSig; isig++){
      fPSC1[detid][segid]->at(ipsc)->spulse[iseg][isig] =
        1./(nhitstmp+1)*(fPSC1[detid][segid]->at(ipsc)->spulse[iseg][isig]*nhitstmp +
                         aps->opulse[iseg][isig]);
    }
  }

  fPSC1[detid][segid]->at(ipsc)->nhits++;

  int tmpnhits = fPSC1[detid][segid]->at(ipsc)->nhits;
  if(tmpnhits>maxnhits) maxnhits=tmpnhits;

  return tmpnhits;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// Optimize Hit position
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
int AGATA::FitPath(Path *apath){
  int ipos = -1;

  float ComptonAngle = apath->GetComptonAngle();

  TVector3 pos0 = apath->GetHit0()->GetPosition();
  TVector3 pos1 = apath->GetHit1()->GetPosition();
  TVector3 pos2 = apath->GetHit2()->GetPosition();

  float foundAngle = 180. / TMath::Pi() * (pos1 - pos0).Angle(pos2 - pos1);
  float angleDiff = fabs(foundAngle - ComptonAngle) + 1.;

  int detid = apath->GetHit1()->GetDet();
  int segid = apath->GetHit1()->GetSeg();

  int npoint = PSCpos[detid][segid].size();

  for(int ipoint=0; ipoint<npoint; ipoint++){
    pos1 = PSCpos[detid][segid][ipoint];
    foundAngle = 180. / TMath::Pi() * (pos1 - pos0).Angle(pos2 - pos1);
    float tmpDiff = fabs(foundAngle - ComptonAngle);

    if(tmpDiff<angleDiff){
      ipos = ipoint;
      angleDiff = tmpDiff;
    }
  }

  return ipos;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// PSA to assign initial Hit pos
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
TVector3 AGATA::GetPSpos(int detid, int segid, PS *aps){
  int ipos = -1;

  if(iPSC0>0){
    char localMask[NSegCore];
    strcpy(localMask,agatageo->GetMask(segid));

    //compare pulse shape
    int npoint = fPSC0[detid][segid]->size();
    float minchi2 = 1e10;
    float asegpulse[NSig_comp], bsegpulse[NSig_comp];

    for(int ipoint=0; ipoint<npoint; ipoint++){
      float chi2=0;
      for(int iseg=0; iseg<NSegCore; iseg++){
	if(localMask[iseg] != '0'){
	  copy_n(aps->opulse[iseg], NSig_comp, asegpulse);
	  copy_n(fPSC0[detid][segid]->at(ipoint)->spulse[iseg], NSig_comp, bsegpulse);
	  float tmpchi2 = Chi2seg(asegpulse, bsegpulse);
	  chi2 += tmpchi2; // sum

	  if(chi2>minchi2) break; // interrupt if exceed minchi2 
	}
      }

      if(chi2<minchi2){
	minchi2 = chi2;
	ipos = ipoint;
      }

    }
  }

  TVector3 pos;
  if(ipos>-1){
    // initial pos from PSA
    TMatrixD LabPos(3,1);
    for(int ix=0; ix<3; ix++) LabPos(ix,0) = fPSC0[detid][segid]->at(ipos)->labpos[ix];
    pos.SetXYZ(LabPos(0,0), LabPos(1,0), LabPos(2,0));
  }else{
    // initial pos at segment center
    TMatrixD SegPos = agatageo->GetSegPos(detid,segid);
    pos.SetXYZ(SegPos(0,0), SegPos(1,0), SegPos(2,0));
  }

  return pos;
}


Float_t AGATA::Chi2seg(const float *apulse, const float *bpulse){

  const __m128 masks = _mm_set1_ps(-0.0f);
  const __m128 zeros = _mm_setzero_ps();
  
  __m128* realtrace = (__m128*)apulse;
  __m128* basetrace = (__m128*)bpulse;

  __m128 diff = _mm_setzero_ps();
  __m128 chis = _mm_setzero_ps();

  for(int nn=0; nn<LOOP_SSE4_seg; nn++){
    diff = _mm_sub_ps(realtrace[nn], basetrace[nn]);
#if   FIXED_METRIC == FIXED_METRIC_ABSVAL   // pow(|d|,1  )
    chis = _mm_add_ps(chis, _mm_andnot_ps(masks, diff));
#elif FIXED_METRIC == FIXED_METRIC_SQUARE   // pow( d ,2  )
    chis = _mm_add_ps(chis, _mm_mul_ps(diff,diff));
#elif FIXED_METRIC == FIXED_METRIC_1SQRT    // pow(|d|,1/2)
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_andnot_ps(masks, diff)));
#elif FIXED_METRIC == FIXED_METRIC_2SQRT    // pow(|d|,1/4)
    chis = _mm_add_ps(chis, _mm_sqrt_ps(_mm_sqrt_ps(_mm_andnot_ps(masks, diff))));
#else
# error Inconsistency in the definition of the distance metric when using the SSE versions
#endif
  }

  chis = _mm_hadd_ps(_mm_hadd_ps(chis, zeros), zeros);

  float chi2 = _mm_cvtss_f32(chis);

  return chi2;
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// tree and file
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
void AGATA::InitTreeWrite(TTree *tree){
  tree->Branch("det",&det);
  tree->Branch("seg",&seg);
  tree->Branch("nhits",&nhits);

  tree->Branch("labpos",labpos,"labpos[3]/F");
  tree->Branch("detpos",detpos,"detpos[3]/F");

  tree->Branch("avelpos",avelpos,"avelpos[3]/F");
  tree->Branch("avedpos",avedpos,"avedpos[3]/F");

  tree->Branch("dist",&dist);

  tree->Branch("spulse",spulse,Form("spulse[%d][%d]/F",NSegCore,NSig));
}

void AGATA::InitTreeRead(TTree *tree){
  tree->SetBranchAddress("det",&det);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("nhits",&nhits);

  tree->SetBranchAddress("labpos",labpos);
  tree->SetBranchAddress("detpos",detpos);

  tree->SetBranchAddress("avelpos",avelpos);
  tree->SetBranchAddress("avedpos",avedpos);

  tree->SetBranchAddress("spulse",spulse);
}


//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
// check memory
//oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo---oooo0000oooo
Double_t AGATA::GetTotalSystemMemory(){
  Long64_t pages = sysconf(_SC_PHYS_PAGES);
  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) pages)*page_size;
}

Double_t AGATA::GetCurrentMemoryUsage(){
  Int_t tSize, resident, share;
  ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();

  Long64_t page_size = sysconf(_SC_PAGE_SIZE);
  return ((Double_t) resident) * page_size;
}


#endif // #ifndef AGATA_CC
