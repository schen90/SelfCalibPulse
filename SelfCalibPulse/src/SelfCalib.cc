#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <mutex>

#include <TROOT.h>
#include <TGraph.h>
#include <TAxis.h>
#include "TInterpreter.h"

#include "Global.hh"
#include "AGATA.hh"
#include "TreeReaderPulse.hh"

using namespace std;

void help(){
  cout<<setw(30)<<left<<" -h"<<" : print options"<<endl
      <<setw(30)<<left<<" -config configurefile"<<" : read G4SimTree from configure file"<<endl
      <<setw(30)<<left<<" -init"<<" : initial folders for selfcalib"<<endl
      <<setw(30)<<left<<" -Fit"<<" : Fit HCs pos"<<endl
      <<setw(30)<<left<<" -loop Nfit"<<" : set iterate Nfit"<<endl
      <<endl;
  return;
}


int main(int argc, char* argv[]){
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");

#ifdef NTHREADS
  ROOT::EnableThreadSafety();
#endif

  // run options
  bool kConfig     = false;
  bool kMakeInit   = false;
  bool kFit        = false;

  string configfile;
  int NFitLoop = 0;
  string PSCPath = "PSCfiles";
  int MaxEvts = -1;

  if(argc==1){
    help();
    return 0;
  }

  for(int i=1; i<argc; i++){
    if(TString(argv[i]) == "-h"){
      help();
      return 0;
    }else if(TString(argv[i]) == "-config"){
      kConfig = true;
      configfile = string(argv[++i]);
    }else if(TString(argv[i]) == "-init"){
      kMakeInit = true;
    }else if(TString(argv[i]) == "-Fit"){
      kFit = true;
    }else if(TString(argv[i]) == "-loop"){
      NFitLoop = atoi(argv[++i]);
    }

  }

  // output Run Options
  cout<<"\e[1;31m Run Opt: ";

  if(kConfig)   cout<<"Config: "<<configfile<<"; ";

  if(kMakeInit) cout<<"ClearFolder; ";

  if(kFit) cout<<"Fit loop "<<" "<<NFitLoop<<"; ";

  cout<<"MaxMem - "<<MaxMemoryUsage<<"%; ";

#ifdef NTHREADS
  if(kFit) cout<<"NthrdsFit - "<<NTHREADS<<"; ";
#endif

  cout<<"\e[0m"<<endl;


  // clock
  time_t start, stop;
  time_t stepstart, stepstop;
  time(&start);

  
  //--ooo000ooo----ooo000ooo----ooo000ooo----ooo000ooo----ooo000ooo----ooo000ooo--
  // start selfcalib procedure
  //--ooo000ooo----ooo000ooo----ooo000ooo----ooo000ooo----ooo000ooo----ooo000ooo--

  // check configure file
  if(!kConfig){
    cerr<<"Please input configure file: -config xxxx"<<endl;
    return 1;
  }
  if(gSystem->AccessPathName(configfile.c_str())){
    cerr<<"Cannot find configure file: "<<configfile<<endl;
    return 1;
  }


  // init TreeReader
  TreeReaderPulse* treereader;
  if(kMakeInit || kFit){
    treereader = new TreeReaderPulse();
    treereader->Load(configfile);
    treereader->SetMaxMemUsage(MaxMemoryUsage); //Max Memory Usage in %
  }

  // clear folders
  if(kMakeInit){
    cout<<"Clear Folders..."<<endl;
    treereader->MakeInit();
  }

  // init AGATA
  AGATA *agata = new AGATA();
  agata->SetMaxMemUsage(MaxMemoryUsage); //Max Memory Usage in %

#ifdef NOISE
  cout<<Form("With noise")<<endl;
#else
  cout<<Form("Without noise")<<endl;
#endif

  //*****************************************//
  // Generate HitCollections
  //*****************************************//
  gROOT->ProcessLine(Form(".!rm -rf %s/iter/it[0-9]*",PSCPath.c_str()));
  if(kFit){
    time(&stepstart);
    agata->InitPSC();
    time(&stepstop);
    printf("=== Create Init PSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

    for(int ifit=0; ifit<NFitLoop; ifit++){
      cout<<"Fit repeat "<<ifit<<" times.."<<endl;

      if(ifit>0){
	time(&stepstart);
	agata->CopyPSC();
	time(&stepstop);
	printf("=== Copy PSC1 to PSC0 time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
      
	time(&stepstart);
	agata->ClearPSC1();
	time(&stepstop);
	printf("=== Clear PSC1 time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
      }
      
      time(&stepstart);
      treereader->GeneratePSC(agata);
      time(&stepstop);
      printf("=== Generate PSC1 time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      time(&stepstart);
      agata->WritePSCfiles(0);
      time(&stepstop);
      printf("=== WritePSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      gROOT->ProcessLine(Form(".!cp -fpdr %s/Det0000.root %s/iter/it/Det0000_fit%d.root",PSCPath.c_str(),PSCPath.c_str(),ifit));

    }

    gROOT->ProcessLine(Form(".!mv -f %s/iter/it %s/iter/it0",PSCPath.c_str(),PSCPath.c_str()));

    
    time(&stepstart);
    agata->WritePSCfiles(-1);
    time(&stepstop);
    printf("=== WritePSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
  }


  time(&stop);
  printf("\n============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  if(kMakeInit || kFit) delete treereader;
  delete agata;

  return 0;
}
