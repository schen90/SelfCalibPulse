//./macros/MakeDataG4_noPS_88Y G4inputfile outputfile

#include "TRint.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TInterpreter.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

using namespace std;

const int Ntype = 3;
const int nsig = 56;
const int nseg = 36;
const int nchan = nseg+1;

double griddist = 2; // 2mm grid
double range[Ntype][3][2];
const Int_t MaxSteps = 50;
int imap[Ntype][MaxSteps][MaxSteps][MaxSteps];

const int NCrystals = 180;
TMatrixD Rt[NCrystals];
TMatrixD Tr[NCrystals];

double threshold = 15; // threshold 15 keV for fired segment
double threshold2 = 300; // threshold 200 keV for core energy (trigger)

void swap(vector<double> &v, int m, int l){
  double temp;
  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
}

void swap(vector<int> &v, int m, int l){
  int temp;
  temp = v[m];
  v[m] = v[l];
  v[l] = temp;
}

void MakeDataG4_noPS_88Y(string G4inputfile, string outputfile){

  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");
  
  // read db
  vector<Int_t>    dbseg[Ntype];
  vector<TMatrixD> dbpos[Ntype];

  string dbfile[3] = {"pulsedb/LibTrap_A001.root","pulsedb/LibTrap_B001.root","pulsedb/LibTrap_C001.root"};

  for(int itype=0; itype<Ntype; itype++){
    for(int ix=0; ix<3; ix++){
      range[itype][ix][0] = 1000;
      range[itype][ix][1] = -1000;
    }

    for(int ix=0; ix<MaxSteps; ix++)
      for(int iy=0; iy<MaxSteps; iy++)
	for(int iz=0; iz<MaxSteps; iz++)
	  imap[itype][ix][iy][iz] = -1;
    
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

    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);
      for(int ix=0; ix<3; ix++){
        if(dbposi[ix]<range[itype][ix][0]) range[itype][ix][0] = dbposi[ix];
        if(dbposi[ix]>range[itype][ix][1]) range[itype][ix][1] = dbposi[ix];
      }
    }

    TMatrixD tmppos(3,1);
    int idx[3];
    for(int ipoint=0; ipoint<npoint; ipoint++){
      dbtree->GetEntry(ipoint);

      for(int i=0; i<3; i++){
	tmppos(i,0)=dbposi[i];
	idx[i] = (int)((dbposi[i]-range[itype][i][0]) / griddist + 0.5);
	if(idx[i]<0 || idx[i]>=MaxSteps){
	  cerr<<"grid point outside Map range!!!"<<endl;
	  return;
	}
      }

      dbseg[itype].push_back(dbsegi);
      dbpos[itype].push_back(tmppos);
      imap[itype][idx[0]][idx[1]][idx[2]] = ipoint;
    }
  
    cout<<"load "<<npoint<<" points from "<<dbfile[itype]<<endl;
    fdb->Close();
  }
  
  // find input LookUpTable
  ifstream fin;
  int dummy_i;
  fin.open("LookUp/CrystalPositionLookUpTable");
  if(!fin){cerr<<"Cannot find LookUp/CrystalPositionLookUpTable"<<endl; return;}
  cout<<"\e[1;32m find CrystalPositionLookUpTable... \e[0m"<<endl;

  int ir;
  for(int i=0; i<NCrystals; i++){
    ir = -1;
    fin >> ir >> dummy_i;
    if(ir<0 || ir>=NCrystals){ ir=i-1; break;}

    Tr[ir].ResizeTo(3,1);  Tr[ir].Zero();
    Rt[ir].ResizeTo(3,3);  Rt[ir].Zero();
    for(int it=0; it<3; it++) fin >> Tr[ir](it,0);
    for(int it=0; it<3; it++){
      fin >> dummy_i;
      for(int it2=0; it2<3; it2++) fin >> Rt[ir](it,it2);
    }
    Rt[ir].Invert();    // change to rot from world frame -> detector frame
  }
  fin.close();
  cout<<"read position for "<<ir+1<<" detectors "<<endl;

  // read data
  FILE *fp;
  char buffer[250];
  fp = fopen(G4inputfile.c_str(),"r");
  if(!fp){  printf("could not open file %s\n",G4inputfile.c_str()); return;}
  printf("\e[1;33m opened file %s \e[0m\n",G4inputfile.c_str());

  while(!feof(fp)){
    memset(buffer,0L,sizeof(buffer));
    fgets(buffer,150,fp);

    if(strncmp(buffer,"$",1)==0) break;
    
    if(feof(fp)){  cerr<<"cannot find data!!!"<<endl; return;}
  }

  
  // output
  TFile *fout = new TFile(outputfile.c_str(),"RECREATE");
  TTree *tree = new TTree("tree","simulate data");

  // G4 vector<>
  int                      ievent; // event id
  vector<int>              ndet;   // detector id
  vector<int>              g4seg;  // segment id
  vector<float>            energy; // deposit energy
  vector<vector<float>>    posa;   // absolute/global interaction position vector<float(3)>
  vector<vector<float>>    posr;   // relative/local interaction position vector<float(3)>

  // pulse shape vector<>
  vector<int>              pdet;   // detector id for pulse shape
  vector<float>            ecore;  // core energy
  vector<vector<int>>      inter;  // interaction id in G4 vector<>

  vector<vector<int>>      iseg;   // fired segment id
  vector<vector<float>>    eseg;   // fired segment energy

  vector<vector<int>>      pseg;   // segment id in pulsedb
  vector<vector<int>>      ngrid;  // number of grid found around ppos
  
  int                      fold;
  bool                     kOneSeg; // at least one det has only one seg fired

  int                      category; // 1: max 1 seg fired in a det, >1 det fired; 2: max >1 seg fired in a det
  
  tree->Branch("ievent",&ievent);
  tree->Branch("ndet",&ndet);
  tree->Branch("g4seg",&g4seg);
  tree->Branch("energy",&energy);
  tree->Branch("posa",&posa);
  tree->Branch("posr",&posr);

  tree->Branch("pdet",&pdet);
  tree->Branch("ecore",&ecore);
  tree->Branch("inter",&inter);

  tree->Branch("iseg",&iseg);
  tree->Branch("eseg",&eseg);

  tree->Branch("pseg",&pseg);
  tree->Branch("ngrid",&ngrid);

  tree->Branch("fold",&fold);
  tree->Branch("kOneSeg",&kOneSeg);

  tree->Branch("category",&category);


  bool                     kskip;

  // clock
  time_t start, stop;
  time(&start);

  int counter=0;
  int readflag = 0;
  int iff, nprev;
  double tmpe, tmpx, tmpy, tmpz;
  ievent = -1; 
  while(!feof(fp)){ // loop events in G4inputfile
    //if(ievent>100) break;
    if(ievent%1000==0)
      cout<<"\r finish "<<ievent<<" events..."<<flush;

    if(fgets(buffer,sizeof(buffer),fp)==NULL){
      printf("end of file\n");
      readflag = 1;
      ievent++;
    }

    sscanf(buffer,"%d",&iff);
    if(iff<-100) continue; //skip -101, -102

    if(iff>-99)
      sscanf(buffer,"%d %lf %lf %lf %lf %d",&iff,&tmpe,&tmpx,&tmpy,&tmpz,&nprev);// read data

    if(iff==-1){
      counter++; // gamma mult counter
    }

    /***************** read interaction points and energies  ****************/
    if(iff>-1){
      ndet.push_back(iff);
      g4seg.push_back(nprev);
      energy.push_back(tmpe);
      TMatrixD mposa(3,1);
      mposa(0,0)=tmpx;  mposa(1,0)=tmpy;  mposa(2,0)=tmpz;
      TMatrixD mposr(3,1);
      mposr = Rt[iff]*(mposa-Tr[iff]);
      vector<float> tmpposa;
      vector<float> tmpposr;
      for(int i=0; i<3; i++){
	tmpposa.push_back(mposa(i,0));
	tmpposr.push_back(mposr(i,0));
      }
      posa.push_back(tmpposa); // global position
      posr.push_back(tmpposr); // local position
    }

    /****************** finish read all gammas in one event ****************/
    if(iff==-100){ readflag=1; counter=1; ievent++;}

    /************************************************************************/
    /*************************** start sim pulse ****************************/
    /************************************************************************/
    if(readflag==1 && energy.size()>0){
      kskip = false;

      // find interaction points in each detector
      for(int i=0; i<ndet.size(); i++){
	int itype = ndet[i]%3;
	//if(itype!=0) continue; // only get pulse db for type 0 det so far

	int id;
	for(id=0; id<pdet.size(); id++) if(ndet[i]==pdet[id]) break;
	if(id==pdet.size()){
	  pdet.push_back(ndet[i]); // fired det id
	  ecore.push_back(0);
	  inter.push_back(vector<int>());
	}
	ecore[id] += energy[i];
	inter[id].push_back(i);
      }
      

      // remove det with ecore<threshold2
      for(int i=0; i<pdet.size(); i++){
	if(ecore[i]<threshold2){
          pdet.erase(pdet.begin()+i);
          ecore.erase(ecore.begin()+i);
          inter.erase(inter.begin()+i);
          i--;
        }
      }

      fold = pdet.size();

      //**************************************************//
      // calculate pulse shape for each detector
      //**************************************************//
      string tmpcout;
      for(int idet=0; idet<pdet.size(); idet++){ // loop detector
	int itype = pdet[idet]%3;

	pseg.push_back(vector<int>()); // fired db seg id in a det
	ngrid.push_back(vector<int>()); // ngrid of every interaction
	
	for(int it=0; it<inter[idet].size(); it++){// loop interactions in a det
	  TMatrixD tmppos(3,1);
	  int idx[3];
	  for(int iaxis=0; iaxis<3; iaxis++){
	    tmppos(iaxis,0)=posr[inter[idet][it]][iaxis];
	    idx[iaxis] = (int)((tmppos(iaxis,0)-range[itype][iaxis][0]) / griddist + 0.5);
	  }

	  // find grid within 2mm cube
	  vector <int> ip; // id list of close grid point for interpolation
	  vector <double> ipdist;
	  double mindist = 5;
	  int minip = -1;
	  int tmpseg = -1;
	  int npoint = dbpos[itype].size();

	  int idxrange = 1;
	  //for(int ipoint=0; ipoint<npoint; ipoint++){
	  for(int ix=idx[0]-idxrange; ix<=idx[0]+idxrange; ix++)
	    for(int iy=idx[1]-idxrange; iy<=idx[1]+idxrange; iy++)
	      for(int iz=idx[2]-idxrange; iz<=idx[2]+idxrange; iz++){
		if(ix<0 || ix>=MaxSteps) continue;
		if(iy<0 || iy>=MaxSteps) continue;
		if(iz<0 || iz>=MaxSteps) continue;

		int ipoint = imap[itype][ix][iy][iz];
		if(ipoint<0) continue;

		if(fabs(tmppos(0,0)-dbpos[itype][ipoint](0,0))>2) continue;
		if(fabs(tmppos(1,0)-dbpos[itype][ipoint](1,0))>2) continue;
		if(fabs(tmppos(2,0)-dbpos[itype][ipoint](2,0))>2) continue;
		double disttmp = sqrt((tmppos-dbpos[itype][ipoint]).Sqr().Sum());
		ip.push_back(ipoint);
		ipdist.push_back(disttmp);
		if(disttmp<mindist){
		  mindist = disttmp;
		  minip = ipoint;
		  tmpseg = dbseg[itype][ipoint]; // seg of closest grid point
		}

	      }

	  if(tmpseg<0){ // remove evt if interaction point outside grid map
	    kskip = true;
	    inter[idet].erase(inter[idet].begin()+it);
	    it--;
	    continue;
	  }
	  pseg[idet].push_back(tmpseg); // segment list for interactions
	  
	  // remove grid from different segment
	  for(int i=0; i<ip.size();){
	    if(dbseg[itype][ip[i]]!=tmpseg){
	      ip.erase(ip.begin()+i);
	    }else{
	      i++;
	    }
	  }
	  
	  ngrid[idet].push_back(ip.size());
	  
	}// end of loop interactions
      }// end of loop pdet


      //**************************************************//
      // find number of fired seg
      //**************************************************//

      int maxfiredseg = 0;
      int nfireddet = 0;
      int nComptondet = 0;
      kOneSeg = false;

      for(int idet=0; idet<pdet.size(); idet++){ // loop detectors

	vector<int>   segid; // fired segment
	vector<float> sege;  // energy deposit in fired segment

	for(int it=0; it<inter[idet].size(); it++){
	  segid.push_back(pseg[idet][it]);
	  sege.push_back(energy[inter[idet][it]]);
	}

	// sum energy in segment
	for(int i=0; i<segid.size(); i++){
	  for(int j=i+1; j<segid.size(); j++){
	    if(segid[i]==segid[j]){
	      sege[i]+=sege[j];
	      sege.erase(sege.begin()+j);
	      segid.erase(segid.begin()+j);
	      j--;	      
	    }
	  }
	}

	// remove seg with sege<threshold
	for(int i=0; i<sege.size(); i++){
	  if(sege[i]<threshold){
	    sege.erase(sege.begin()+i);
	    segid.erase(segid.begin()+i);
	    i--;
	  }
	}

	iseg.push_back(segid);
	eseg.push_back(sege);

	int nfiredseg = sege.size();
	int nComptonseg = 0;
	for(int i=0; i<sege.size(); i++){
	  if(sege[i]>threshold2){
	    nComptonseg++;
	  }
	}// end of loop seg

        if(nfiredseg==1 && nComptonseg==1) kOneSeg = true;
	
	if(nfiredseg>maxfiredseg) maxfiredseg = nfiredseg;
	if(nfiredseg>0) nfireddet++;
	if(nComptonseg>0) nComptondet++;
	
      }// end of loop det

      if(!kskip){
	category = 0;
	//if(maxfiredseg==1 && nfireddet>1) category = 1;
	if(maxfiredseg==1 && nComptondet>1) category = 1;
        if(maxfiredseg>1 && nComptondet>1 && kOneSeg) category = 2;
      
	//if(category>0) tree->Fill(); // fill only evt with >1 seg fired
	//if(category==1) tree->Fill(); // fill only evt with max 1 seg fired in a det, >1 det fired
	//if(category==1 || category==2) tree->Fill(); //
	if(fold>1) tree->Fill();

	//if(category==1) cout<<tmpcout;
      }

      /***************************** init *****************************/
      ndet.clear();
      g4seg.clear();
      energy.clear();
      posa.clear();
      posr.clear();

      pdet.clear();
      ecore.clear();
      inter.clear();

      iseg.clear();
      eseg.clear();
      
      pseg.clear();
      ngrid.clear();

    } // end of one event
    readflag = 0;
    
  }
  cout<<"\r finish "<<ievent<<" events..."<<endl;
  fclose(fp);
  
  time(&stop);
  printf("============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  cout<<"Write to "<<outputfile<<endl;
  
  fout->cd();
  tree->Write();
  fout->Close();

  return;
}

#ifndef __CINT__
int main(int argc, char *argv[]){

  if(argc>2){
    MakeDataG4_noPS_88Y( string(argv[1]), string(argv[2]));
  }else if(argc>1){
    MakeDataG4_noPS_88Y( string(argv[1]), "rootfiles/noPS/G4SimData0000.root" );
  }else{
    MakeDataG4_noPS_88Y( "trunk/GammaEvents.0000", "rootfiles/noPS/G4SimData0000.root" );
  }

  return 0;
}
#endif
