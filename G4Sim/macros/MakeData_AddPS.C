//./macros/MakeData_AddPS inputfile outputfile

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

bool knoise = false;
bool kextrapol = true;
bool kgridip = true;

const int Ntype = 1;
const int nsig = 56;
const int nseg = 36;
const int nchan = nseg+1;

double griddist = 2; // 2mm grid
double range[Ntype][3][2];
const Int_t MaxSteps = 50;
int imap[Ntype][MaxSteps][MaxSteps][MaxSteps];

double threshold = 20; // threshold 20 keV for fired segment
double threshold2 = 300; // threshold 300 keV for Compton scattering

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

void MakeData_AddPS(string inputfile, string outputfile){

  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");
  
  // read db
  vector<Int_t>    dbseg[Ntype];
  vector<TMatrixD> dbpos[Ntype];
  vector<TMatrixD> dbcore[Ntype];
  vector<TMatrixD> dbspulse[Ntype];

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
    Float_t dbspulsei[nsig*nchan];

    dbtree->SetBranchAddress("seg",&dbsegi);
    dbtree->SetBranchAddress("pos",dbposi);
    dbtree->SetBranchAddress("spulse",dbspulsei);
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
    TMatrixD tmpcore(nsig,1);
    TMatrixD tmpspulse(nsig*nseg,1);
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

      for(int i=0; i<nsig; i++) tmpcore(i,0)=dbspulsei[nseg*nsig+i];
      for(int iseg=0; iseg<nseg; iseg++){
	for(int i=0; i<nsig; i++)
	  tmpspulse(iseg*nsig+i,0)=dbspulsei[iseg*nsig+i];
      }

      dbseg[itype].push_back(dbsegi);
      dbpos[itype].push_back(tmppos);
      imap[itype][idx[0]][idx[1]][idx[2]] = ipoint;

      dbcore[itype].push_back(tmpcore);
      dbspulse[itype].push_back(tmpspulse);
    }
  
    cout<<"load "<<npoint<<" points from "<<dbfile[itype]<<endl;
    fdb->Close();
  }
  
  //--ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo--
  // read data
  TFile *fin = new TFile(inputfile.c_str());
  if(!fin->IsOpen()){  printf("could not open file %s\n",inputfile.c_str()); return;}
  printf("\e[1;33m opened file %s \e[0m\n",inputfile.c_str());

  TTree *oldtree = (TTree *)fin->Get("tree");
  // read in branch
  int                      oievent;     // event id
  vector<int>             *ondet = 0;   // detector id
  vector<int>             *og4seg = 0;  // segment id
  vector<float>           *oenergy = 0; // deposit energy
  vector<vector<float>>   *oposa = 0;   // absolute/global interaction position vector<float(3)>
  vector<vector<float>>   *oposr = 0;   // relative/local interaction position vector<float(3)>
  
  // pulse shape vector<>
  vector<int>             *opdet = 0;   // detector id for pulse shape
  vector<float>           *oecore = 0;  // core energy
  vector<vector<int>>     *ointer = 0;  // interaction id in G4 vector<>
  
  oldtree->SetBranchAddress("ievent",&oievent);
  oldtree->SetBranchAddress("ndet",&ondet);
  oldtree->SetBranchAddress("g4seg",&og4seg);
  oldtree->SetBranchAddress("energy",&oenergy);
  oldtree->SetBranchAddress("posa",&oposa);
  oldtree->SetBranchAddress("posr",&oposr);

  oldtree->SetBranchAddress("pdet",&opdet);
  oldtree->SetBranchAddress("ecore",&oecore);
  oldtree->SetBranchAddress("inter",&ointer);

  //--ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo----ooo0000ooo--

  
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

  vector<vector<int>>      pseg;   // segment id in pulsedb
  vector<vector<int>>      ngrid;  // number of grid found around ppos
  vector<vector<int>>      extrpl; // if use extrapolation
  vector<vector<float>>    core;   // core pulse shape vector<float(56)> 
  vector<vector<float>>    spulse; // segment pulse shape vector<float(2016)>
  vector<vector<int>>      gridip; // grid points used for spulse
  vector<vector<float>>    gridwgt; // weight grid points used for spulse
  
  int                      category; // 1: max 1 seg fired in a det, >1 det fired; 2: max >1 seg fired in a det
  bool                     kskip;
  
  tree->Branch("ievent",&ievent);
  tree->Branch("ndet",&ndet);
  tree->Branch("g4seg",&g4seg);
  tree->Branch("energy",&energy);
  tree->Branch("posa",&posa);
  tree->Branch("posr",&posr);

  tree->Branch("pdet",&pdet);
  tree->Branch("ecore",&ecore);
  tree->Branch("inter",&inter);

  tree->Branch("pseg",&pseg);
  tree->Branch("ngrid",&ngrid);
  if(kextrapol){
    tree->Branch("extrpl",&extrpl);
  }
  tree->Branch("core",&core);
  tree->Branch("spulse",&spulse);
  if(kgridip){
    tree->Branch("gridip",&gridip);
    tree->Branch("gridwgt",&gridwgt);
  }
  
  tree->Branch("category",&category);

  // clock
  time_t start, stop;
  time(&start);

  int nentries = oldtree->GetEntriesFast();
  int ievt;
  for(ievt=0; ievt<nentries; ievt++){
    if(ievt%1000==0)  cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<flush;
    oldtree->GetEntry(ievt);

    /***************** read interaction points and energies  ****************/
    ievent = oievent;
    for(int i=0; i<ondet->size(); i++){
      ndet.push_back(ondet->at(i));
      g4seg.push_back(og4seg->at(i));
      energy.push_back(oenergy->at(i));
      posa.push_back(oposa->at(i)); // global position
      posr.push_back(oposr->at(i)); // local position
    }

    /************************************************************************/
    /*************************** start sim pulse ****************************/
    /************************************************************************/
    if(energy.size()>0){
      kskip = false;

      // find interaction points in each detector
      for(int i=0; i<opdet->size(); i++){
	int itype = opdet->at(i)%3;
	if(itype!=0) continue; // only get pulse db for type 0 det so far
	
	pdet.push_back(opdet->at(i));
	ecore.push_back(oecore->at(i));
	inter.push_back(ointer->at(i));
      }
      
      //**************************************************//
      // calculate pulse shape for each detector
      //**************************************************//
      string tmpcout;
      for(int idet=0; idet<pdet.size(); idet++){ // loop detector
	int itype = pdet[idet]%3;

	pseg.push_back(vector<int>()); // fired db seg id in a det
	ngrid.push_back(vector<int>()); // ngrid of every interaction
	if(kextrapol) extrpl.push_back(vector<int>()); // extrapolate npoint of every interaction

	if(kgridip){
	  gridip.push_back(vector<int>()); // grid ip list for every detector
	  gridwgt.push_back(vector<float>()); // weight of grid list for every detector
	}
	core.push_back(vector<float>(nsig)); // core pulse shape for every detector
	spulse.push_back(vector<float>(nsig*nseg)); // seg pulse shape for every detector
	
	for(int it=0; it<inter[idet].size(); it++){// loop interactions in a det
	  TMatrixD tmppos(3,1);
	  int idx[3];
	  for(int iaxis=0; iaxis<3; iaxis++){
	    tmppos(iaxis,0)=posr[inter[idet][it]][iaxis];
            idx[iaxis] = (int)((tmppos(iaxis,0)-range[itype][iaxis][0]) / griddist + 0.5);
	  }

	  // find grid within 2mm cube
	  vector <int> ip; // id list of close grid point for interpolation
	  vector <int> ip2; // id list of close grid point 4mm for extrapolation
	  vector <double> ipdist;
	  vector <double> ip2dist;
	  vector <int> exflag;
	  double mindist = 5;
	  int minip = -1;
	  int tmpseg = -1;
	  int npoint = dbpos[itype].size();

	  int idxrange = 2;
	  //for(int ipoint=0; ipoint<npoint; ipoint++){
	  for(int ix=idx[0]-idxrange; ix<=idx[0]+idxrange; ix++)
            for(int iy=idx[1]-idxrange; iy<=idx[1]+idxrange; iy++)
              for(int iz=idx[2]-idxrange; iz<=idx[2]+idxrange; iz++){
                if(ix<0 || ix>=MaxSteps) continue;
                if(iy<0 || iy>=MaxSteps) continue;
                if(iz<0 || iz>=MaxSteps) continue;

                int ipoint = imap[itype][ix][iy][iz];
                if(ipoint<0) continue;

		if(kextrapol){
		  if(fabs(tmppos(0,0)-dbpos[itype][ipoint](0,0))>4) continue;
		  if(fabs(tmppos(1,0)-dbpos[itype][ipoint](1,0))>4) continue;
		  if(fabs(tmppos(2,0)-dbpos[itype][ipoint](2,0))>4) continue;
		  double disttmp = sqrt((tmppos-dbpos[itype][ipoint]).Sqr().Sum());
		  ip2.push_back(ipoint);
		  ip2dist.push_back(disttmp);
		}
	    
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

		if(kextrapol){ ip2.pop_back(); ip2dist.pop_back();}
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
	      ipdist.erase(ipdist.begin()+i);
	    }else{
	      exflag.push_back(0);
	      i++;
	    }
	  }
	  // sort according to dist
	  for(int i=0; i<ip.size(); i++){
	    for(int j=i+1; j<ip.size(); j++){
	      if(ipdist[i] > ipdist[j]){
		swap(ipdist, i, j);
		swap(ip, i, j);
	      }
	    }
	  }
	  
	  // extrapolation ---------------------------------------------------
	  if(kextrapol){
	    extrpl[idet].push_back(0);
	    if(ip.size()<8){
	      // remove grid from different segment
	      for(int i=0; i<ip2.size();){
		if(dbseg[itype][ip2[i]]!=tmpseg){
		  ip2.erase(ip2.begin()+i);
		  ip2dist.erase(ip2dist.begin()+i);
		}else{
		  i++;
		}
	      }
	      // sort according to dist
	      for(int i=0; i<ip2.size(); i++){
		for(int j=i+1; j<ip2.size(); j++){
		  if(ip2dist[i] > ip2dist[j]){
		    swap(ip2dist, i, j);
		    swap(ip2, i, j);
		  }
		}
	      }	      
	      // add extrapolation points
	      for(int i=0; i<ip2.size();i++){
		ip.push_back(ip2[i]);
		exflag.push_back(1);
	      }
	    }
	  }
	  // -----------------------------------------------------------------

	  if(mindist==0){

	    ngrid[idet].push_back(1);
	    for(int isig=0; isig<nsig; isig++) core[idet][isig]+=dbcore[itype][minip](isig,0)*energy[inter[idet][it]];
	    for(int iseg=0; iseg<nseg; iseg++)
	      for(int i=0; i<nsig; i++)
		spulse[idet][iseg*nsig+i]+=dbspulse[itype][minip](iseg*nsig+i,0)*energy[inter[idet][it]];
	    if(kgridip){
	      gridip[idet].push_back(minip);
	      gridwgt[idet].push_back(1);
	    }

	  }else{// more than 1 grid

	    // prepare spulse list
	    vector <TMatrixD> corelist;
	    vector <TMatrixD> spulselist;
	    vector <TMatrixD> poslist;
	    vector <vector<int>> igridlist;
	    vector <vector<float>> gridwgtlist;
	    for(int i=0; i<ip.size(); i++){
	      corelist.push_back(dbcore[itype][ip[i]]);
	      spulselist.push_back(dbspulse[itype][ip[i]]);
	      poslist.push_back(dbpos[itype][ip[i]]);
	      igridlist.push_back(vector<int>(1));
	      igridlist[i][0] = i;
	      gridwgtlist.push_back(vector<float>(1));
	      gridwgtlist[i][0] = 1;
	    }

	    // linear combine pulse
	    for(int iaxis=0; iaxis<3; iaxis++){ //combine iaxis direction
	      int iaxis1 = (iaxis+1)%3;
	      int iaxis2 = (iaxis+2)%3;

	      for(int i=0; i<poslist.size(); ){
		for(int j=i+1; j<poslist.size(); ){
		  bool kcombine = true;
		  if(fabs(poslist[i](iaxis1,0)-poslist[j](iaxis1,0))>0.1) kcombine = false;
		  if(fabs(poslist[i](iaxis2,0)-poslist[j](iaxis2,0))>0.1) kcombine = false;

		  if(kcombine){
		    double factori = (tmppos(iaxis,0)-poslist[j](iaxis,0)) / (poslist[i](iaxis,0)-poslist[j](iaxis,0));
		    double factorj = (tmppos(iaxis,0)-poslist[i](iaxis,0)) / (poslist[j](iaxis,0)-poslist[i](iaxis,0));
		    corelist[i]   = factori*corelist[i] + factorj*corelist[j];
		    spulselist[i] = factori*spulselist[i] + factorj*spulselist[j];
		    poslist[i](iaxis,0) = tmppos(iaxis,0);

		    if(factori==0){
		      igridlist[i].clear();
		      gridwgtlist[i].clear();
		    }else{
		      for(int ig=0; ig<gridwgtlist[i].size(); ig++)
			gridwgtlist[i][ig] = gridwgtlist[i][ig]*factori;
		    }
		    
		    if(factorj!=0)
		      for(int ig=0; ig<igridlist[j].size(); ig++){
			igridlist[i].push_back(igridlist[j][ig]);
			gridwgtlist[i].push_back(gridwgtlist[j][ig]*factorj);
		      }

		    corelist.erase(corelist.begin()+j);
		    spulselist.erase(spulselist.begin()+j);
		    poslist.erase(poslist.begin()+j);
		    igridlist.erase(igridlist.begin()+j);
		    gridwgtlist.erase(gridwgtlist.begin()+j);

		  }else{
		    j++;
		  }
		}// end of loop j
		
		if(!kextrapol) poslist[i](iaxis,0) = tmppos(iaxis,0); // reduced interpolation
		i++;
	      }// end of loop i
	    }// end of loop iaxis

	    if(poslist.size()>1){

	      if(!kextrapol){
		cerr<<"something wrong..."<<endl;

	      }else{ // if use extrapolation

		while(poslist.size()>0){
		  if(fabs(tmppos(0,0)-poslist[0](0,0))<0.01 &&
		     fabs(tmppos(1,0)-poslist[0](1,0))<0.01 &&
		     fabs(tmppos(2,0)-poslist[0](2,0))<0.01)
		    break;

		  corelist.erase(corelist.begin());
		  spulselist.erase(spulselist.begin());
		  poslist.erase(poslist.begin());
		  igridlist.erase(igridlist.begin());
		  gridwgtlist.erase(gridwgtlist.begin());
		}

		if(poslist.size()==0){ // cannot get simulated position in XYZ combine
		  vector<vector<int>> comb;

		  for(int i=0; i<ip.size(); i++){
		    corelist.push_back(dbcore[itype][ip[i]]);
		    spulselist.push_back(dbspulse[itype][ip[i]]);
		    poslist.push_back(dbpos[itype][ip[i]]);

		    igridlist.push_back(vector<int>(1));
		    igridlist[i][0] = i;

		    gridwgtlist.push_back(vector<float>(1));
		    gridwgtlist[i][0] = 1;

		    comb.push_back(vector<int>(3));
		    comb[i][0] = comb[i][1] = comb[i][2] = -1;
		  }

		  /*
		  // output extrapolation
		  tmpcout += (string)Form("\n poslist.size = %d \n",poslist.size());
		  tmpcout += (string)Form("simpos : %.2f  %.2f  %.2f \n",tmppos(0,0),tmppos(1,0),tmppos(2,0));
		  tmpcout += "grids:\n";
		  for(int i=0; i<poslist.size(); i++)
		    tmpcout += (string)Form("%d       %.2f  %.2f  %.2f \n",exflag[i],poslist[i](0,0),poslist[i](1,0),poslist[i](2,0));
		  tmpcout += "\n";
		  */
		  
		  for(int i=0; i<poslist.size(); i++){ // loop points

		    int n = poslist.size();

		    for(int iaxis=0; iaxis<3; iaxis++){
		      if(comb[i][iaxis]>=0) continue; // i used in iaxis

		      int iaxis1 = (iaxis+1)%3;
		      int iaxis2 = (iaxis+2)%3;

		      for(int j=i+1; j<n; j++){
			if(comb[j][iaxis]>=0) continue; // j used in iaxis

			if(fabs(poslist[i](iaxis1,0)-poslist[j](iaxis1,0))>0.1) continue;
			if(fabs(poslist[i](iaxis2,0)-poslist[j](iaxis2,0))>0.1) continue;
			comb[j][iaxis] = 0;

			if(comb[i][iaxis]<0){ // combine i,j, push back
			  comb[i][iaxis] = 0;

			  TMatrixD tmpposlist(3,1);
			  tmpposlist(iaxis,0) = tmppos(iaxis,0);
			  tmpposlist(iaxis1,0) = poslist[i](iaxis1,0);
			  tmpposlist(iaxis2,0) = poslist[i](iaxis2,0);

			  bool kadd = true;
			  for(int jj=0; jj<n; jj++){
			    if(fabs(tmpposlist(iaxis,0)-poslist[jj](iaxis,0))>0.01) continue;
			    if(fabs(tmpposlist(iaxis1,0)-poslist[jj](iaxis1,0))>0.01) continue;
			    if(fabs(tmpposlist(iaxis2,0)-poslist[jj](iaxis2,0))>0.01) continue;
			    kadd = false; // point exist
			  }
			  if(!kadd) continue;
			  
			  // add new point
			  double factori = (tmppos(iaxis,0)-poslist[j](iaxis,0)) / (poslist[i](iaxis,0)-poslist[j](iaxis,0));
			  double factorj = (tmppos(iaxis,0)-poslist[i](iaxis,0)) / (poslist[j](iaxis,0)-poslist[i](iaxis,0));
			  TMatrixD tmpcorelist = factori*corelist[i] + factorj*corelist[j];
			  TMatrixD tmpspulselist = factori*spulselist[i] + factorj*spulselist[j];

			  vector<int> tmpigridlist;
			  vector<float> tmpgridwgtlist;
			  if(factori!=0)
			    for(int ig=0; ig<igridlist[i].size(); ig++){
			      tmpigridlist.push_back(igridlist[i][ig]);
			      tmpgridwgtlist.push_back(gridwgtlist[i][ig]*factori);
			    }
			  if(factorj!=0)
			    for(int ig=0; ig<igridlist[j].size(); ig++){
			      tmpigridlist.push_back(igridlist[j][ig]);
			      tmpgridwgtlist.push_back(gridwgtlist[j][ig]*factorj);
			    }

			  vector<int> tmpcomb;
			  for(int iiaxis=0; iiaxis<3; iiaxis++)
			    tmpcomb.push_back(comb[i][iiaxis]);
			  tmpcomb[iaxis] = 1;

			  corelist.push_back(tmpcorelist);
			  spulselist.push_back(tmpspulselist);
			  poslist.push_back(tmpposlist);
			  igridlist.push_back(tmpigridlist);
			  gridwgtlist.push_back(tmpgridwgtlist);
			  comb.push_back(tmpcomb);
			}// end of combine i,j

		      }// end of loop poslist j

		    }// end of loop iaxis

		    if((comb[i][0] + comb[i][1] + comb[i][2])==3) break;

		    if(poslist.size()>1){
		      corelist.erase(corelist.begin());
		      spulselist.erase(spulselist.begin());
		      poslist.erase(poslist.begin());
		      igridlist.erase(igridlist.begin());
		      gridwgtlist.erase(gridwgtlist.begin());
		      comb.erase(comb.begin());
		      i--;
		    }
		  }// end of loop poslist i

		  if(poslist.size()<1) cerr<<"something wrong...poslist.size() = "<<poslist.size()<<endl;

		  extrpl[idet].back()+=1000;

		  if((comb[0][0] + comb[0][1] + comb[0][2])<3){
		    extrpl[idet].back() += 10000;
		    kskip = true; // remove evt cannot get simulated position
		  }

		}// end of cannot get simulated position in XYZ combine
		
	      }// end of kextrapol

	    }// end of poslist>1

	    ngrid[idet].push_back(igridlist[0].size());
	    for(int ig=0; ig<igridlist[0].size(); ig++){
	      if(kgridip){
		gridip[idet].push_back(ip[igridlist[0][ig]]);
		gridwgt[idet].push_back(gridwgtlist[0][ig]);
	      }
	      if(kextrapol){
		if(exflag[igridlist[0][ig]]==1) extrpl[idet].back()++;
	      }
	    }
	    
	    for(int isig=0; isig<nsig; isig++) core[idet][isig]+=corelist[0](isig,0)*energy[inter[idet][it]];
	    for(int iseg=0; iseg<nseg; iseg++)
	      for(int i=0; i<nsig; i++)
		spulse[idet][iseg*nsig+i]+=spulselist[0](iseg*nsig+i,0)*energy[inter[idet][it]];

	  } // end of more than 1 grid

	}// end of loop interactions

	if(knoise){
	  // add noise here
	  float noise[4477],noise2[4477];
	  float tmpnoise;
	  for(int iseg=0; iseg<nseg+1; iseg++){
	    tmpnoise = gRandom->Uniform(-1,1);
	    for(int i=0; i<nsig; i++){
	      tmpnoise += -0.2*(tmpnoise+gRandom->Uniform(-10,10));
	      noise[iseg*nsig+i]=tmpnoise;
	    }
	  }

	  // smooth noise
	  int nsmooth=3;
	  for(int i=0; i<4477; i++){
	    noise2[i] = noise[i];
	    if(i<nsmooth || i>=4477-nsmooth) continue;
	    for(int j=1; j<=nsmooth; j++){
	      noise2[i]+=noise[i+j]+noise[i-j];
	    }
	    noise2[i]=noise2[i]/(2*nsmooth+1)*2;// scale noise
	  }
	
	  for(int i=0; i<nsig; i++) core[idet][i]+=noise2[i];
	  for(int i=0; i<nseg*nsig; i++) spulse[idet][i]+=noise2[i+nsig];
	}
	
      }// end of loop pdet


      // find number of fired seg
      int maxfiredseg = 0;
      int nfireddet = 0;
      int nComptondet = 0;
      bool kOneSegFired = false;
      for(int idet=0; idet<spulse.size(); idet++){

	int nfiredseg = 0;
	int nComptonseg = 0;

	for(int iseg=0; iseg<nseg; iseg++){
	  double emax = 0;
	  for(int isig=nsig/4*3; isig<nsig; isig++){
	    double etmp = spulse[idet][iseg*nsig+isig];
	    if(etmp>emax) emax=etmp; // max e of a seg
	  }

	  if(emax>threshold){
	    int counter=0;
	    int counter2=0;
	    for(int isig=nsig/4*3; isig<nsig; isig++){
	      double etmp = spulse[idet][iseg*nsig+isig];
	      if(etmp>threshold && fabs(etmp-emax)<threshold) counter++;
	      if(etmp>threshold2 && fabs(etmp-emax)<threshold) counter2++;
	    }
	    if(counter>15) nfiredseg++;
	    if(counter2>15) nComptonseg++;
	  }
	}// end of loop seg

        if(nfiredseg==1 && nComptonseg==1) kOneSegFired = true;
	
	if(nfiredseg>maxfiredseg) maxfiredseg = nfiredseg;
	if(nfiredseg>0) nfireddet++;
	if(nComptonseg>0) nComptondet++;
	
      }// end of loop det

      if(!kskip){
	category = 0;
	//if(maxfiredseg==1 && nfireddet>1) category = 1;
	if(maxfiredseg==1 && nComptondet>1) category = 1;
	if(maxfiredseg>1 && nComptondet>1 && kOneSegFired) category = 2;
      
	//if(category>0) tree->Fill(); // fill only evt with >1 seg fired
	//if(category==1) tree->Fill(); // fill only evt with max 1 seg fired in a det, >1 det fired
	if(category==1 || category==2) tree->Fill(); //

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

      pseg.clear();
      ngrid.clear();
      extrpl.clear();
      core.clear();
      spulse.clear();
      gridip.clear();
      gridwgt.clear();
      
    } // end of one event
    
  }
  cout<<"\r finish "<<ievt<<" / "<<nentries<<" events..."<<endl;
  cout<<"last event id in Geant4 "<<ievent<<endl;
  fin->Close();
  
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
    MakeData_AddPS( string(argv[1]), string(argv[2]) );
  }else if(argc>1){
    MakeData_AddPS( string(argv[1]), "rootfiles/PS/G4SimData0000.root" );
  }else{
    MakeData_AddPS( "rootfiles/noPS/G4SimData0000.root", "rootfiles/PS/G4SimData0000.root" );
  }

  return 0;
}
#endif
