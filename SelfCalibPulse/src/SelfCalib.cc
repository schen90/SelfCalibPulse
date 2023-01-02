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
      <<setw(30)<<left<<" -config configurefile"<<" : read Tree from configure file"<<endl
      <<setw(30)<<left<<" -init"<<" : initial folders for selfcalib"<<endl
      <<setw(30)<<left<<" -gp i"<<" : group pulse shape for det i.  i=-1 : all det"<<endl
      <<setw(30)<<left<<" -PSA"<<" : PSA to assign initial hit pos"<<endl
      <<setw(30)<<left<<" -comb"<<" : combine Hit files for every run"<<endl
      <<setw(30)<<left<<" -Fit"<<" : Fit HCs pos"<<endl
      <<setw(30)<<left<<" -loop Ntrack Nfit"<<" : set iterate Ntrack Nfit"<<endl
      <<endl;
  return;
}

int main(int argc, char* argv[]){
  
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");

#ifdef NTHREADS
  ROOT::EnableThreadSafety();
#endif

  //gROOT->ProcessLine("gErrorIgnoreLevel = 2002;");

  // run options
  bool kConfig     = false;
  bool kMakeInit   = false;
  bool kGroupPulse = false;
  bool kPSA        = false;
  bool kComb       = false;
  bool kFit        = false;
  
  string configfile;
  int Detid = -1;
  int NFullLoop = 0;
  int NFitLoop = 0;
  string PSCPath = "PSCfiles";
  int MaxEvts = -1;
  double Diff = -1;
  
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
    }else if(TString(argv[i]) == "-gp"){
      kGroupPulse = true;
      Detid = atoi(argv[++i]);
    }else if(TString(argv[i]) == "-PSA"){
      kPSA = true;
    }else if(TString(argv[i]) == "-comb"){
      kComb = true;
    }else if(TString(argv[i]) == "-Fit"){
      kFit = true;
    }else if(TString(argv[i]) == "-loop"){
      NFullLoop = atoi(argv[++i]);
      NFitLoop = atoi(argv[++i]);
    }

  }

  // output Run Options
  cout<<"\e[1;31m Run Opt: ";

  if(kConfig)   cout<<"Config: "<<configfile<<"; ";

  if(kMakeInit) cout<<"ClearFolder; ";

  if(kGroupPulse){
    cout<<"Group";
    if(Detid>-1) cout<<" Det "<<Detid;
    else         cout<<" AllDet";
    cout<<"; ";

    if(kPSA) cout<<"PSA Radius0 "<<RADIUS0<<"; ";
  }

  if(kComb) cout<<"Combine Hitfiles; ";
  if(kFit) cout<<"Fit loop "<<NFullLoop<<" "<<NFitLoop<<"; ";
  
  cout<<"MaxMem - "<<MaxMemoryUsage<<"%; ";

#ifdef NTHREADS
  if(kGroupPulse) cout<<"NthrdsGP - "<<NTHREADS<<"; ";
#endif
#ifdef NTHREADS2
  if(kFit) cout<<"NthrdsFit - "<<NTHREADS2<<"; ";
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

  // check run opt
  if(kFit && Detid>-1){
    cerr<<"Cannot fit HC pos for one det!!!"<<endl;
    return 1;
  }
  
  // init TreeReader
  TreeReaderPulse* treereader;
  if(kMakeInit || kGroupPulse){
    treereader = new TreeReaderPulse(Detid);
    treereader->Load(configfile);
    treereader->SetMaxMemUsage(MaxMemoryUsage); //Max Memory Usage in %
  }

  // clear folders
  if(kMakeInit){
    cout<<"Clear Folders..."<<endl;
    treereader->MakeInit();
  }


  // init AGATA
  AGATA *agata = new AGATA(Detid);  
  agata->SetMaxMemUsage(MaxMemoryUsage); //Max Memory Usage in %
  agata->SetPSA(kPSA);


  //*****************************************//
  // Generate HitCollections
  //*****************************************//
  if(kGroupPulse){

    time(&stepstart);
    if(kConfig) treereader->GenerateHCs(0,agata);
    time(&stepstop);
    printf("=== InitialHCs time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

    int ndiv = 0;
    long long PSCstat[10];
    agata->CheckPSCstat(PSCstat);
    while(PSCstat[0]>MAXHITS){
      cout<<"\033[1;31m"<<"Divide "<<ndiv<<": \033[0m"<<endl;

      int DivDir = ndiv%3;
      if(DivDir<0) cout<<"divide in all direction..."<<endl;
      else         cout<<"divide in direction "<<DivDir<<endl;
      if(kConfig) agata->SetDivDir(DivDir);
      
      time(&stepstart);
      if(kConfig) treereader->GenerateHCs(1,agata);
      if(kConfig) agata->FindDevCut();
      time(&stepstop);
      printf("=== FindMaxDeviation time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      time(&stepstart);
      if(kConfig) treereader->GenerateHCs(2,agata);
      time(&stepstop);
      printf("=== Divide PSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      time(&stepstart);
      if(kConfig) agata->RemoveMotherPSC();
      if(kConfig) agata->RemoveSmallPSC(MINHITS);
      time(&stepstop);
      printf("=== Remove PSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
      
      agata->CheckPSCstat(PSCstat);
      cout<<"\033[1m"<<"PSC stats:"
	  <<"  maxnhits = "<<PSCstat[0]
	  <<" ; nPSC = "<<PSCstat[1]
	  <<" ; nEmpty = "<<PSCstat[2]
	  <<"\033[0m"<<endl<<endl;
      ndiv++;
    }

    time(&stepstart);
    if(kConfig) treereader->GenerateHCs(3,agata);
    if(kConfig) agata->CalcDevSigma();
    time(&stepstop);
    printf("=== FindDevSigma time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
    
    // remove strange PS from PSC
    int nremove = 0;
    bool kRemove = true;
    while(kRemove && nremove<0){
      cout<<"\033[1;31m"<<"Remove strange PS "<<nremove<<": \033[0m"<<endl;
      if(nremove<2) agata->SetNSigma(2.);
      else          agata->SetNSigma(3.);
      
      time(&stepstart);
      if(kConfig) agata->MakeCPulse();
      if(kConfig) treereader->GenerateHCs(4,agata);
      time(&stepstop);
      printf("=== Remove strange PS time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      kRemove = false;
      if( treereader->GetRemovePSNumber() > 0 ) kRemove = true;

      time(&stepstart);
      if(kConfig) agata->RemoveSmallPSC(MINHITS);
      time(&stepstop);
      printf("=== Remove PSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
      
      agata->CheckPSCstat(PSCstat);
      cout<<"\033[1m"<<"PSC stats:"
	  <<"  maxnhits = "<<PSCstat[0]
	  <<" ; nPSC = "<<PSCstat[1]
	  <<" ; nEmpty = "<<PSCstat[2]
	  <<"\033[0m"<<endl<<endl;

      time(&stepstart);
      if(kConfig) treereader->GenerateHCs(3,agata);
      if(kConfig) agata->CalcDevSigma();
      time(&stepstop);
      printf("=== FindDevSigma time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      nremove++;
    }

    
    // save grouped HCs and Hits
    time(&stepstart);
    agata->WritePSCfiles(Detid);
    agata->WriteHCfiles(Detid);
    agata->WriteEvtHitsfiles(Detid);
    time(&stepstop);
    printf("=== Write initial PSC&HC&Hits files time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
  }


  if(kComb){
    time(&stepstart);
    agata->Load(configfile);
    agata->CombEvtHitsfiles();
    time(&stepstop);
    printf("=== Combine Hitfiles time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
  }

  
  //*****************************************//
  // HitCollections position optimization
  //*****************************************//
  if(kFit){

    if(!kGroupPulse){
      // read data from PSCfiles and HCfiles
      time(&stepstart);
      agata->Load(configfile);
      agata->LoadHCfiles();
      agata->LoadEvtHitsconfigs();
      agata->LoadPSCfiles(0);
      time(&stepstop);
      printf("=== Load HC&EventHits time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
    }
    

    // HCs position optimize
    gROOT->ProcessLine(Form(".!rm -rf %s/iter/it[0-9]*",PSCPath.c_str()));
    for(int FullCounter=0; FullCounter<NFullLoop; FullCounter++){

      cout<<"\e[1;31m"<<"Full iteration "<<FullCounter<<" : \e[0m"<<endl;

      time(&stepstart);
      agata->Tracking(FullCounter);
      time(&stepstop);
      printf("=== Tracking time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

      time(&stepstart);
      agata->RegisterPathswithHCs();
      time(&stepstop);
      printf("=== Register Paths time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
  
      time(&stepstart);
      agata->ExecFit(NFitLoop);
      gROOT->ProcessLine(Form(".!mv -f %s/iter/it %s/iter/it%d",PSCPath.c_str(),PSCPath.c_str(),FullCounter));
      gROOT->ProcessLine(Form(".!mkdir %s/iter/it",PSCPath.c_str()));
      time(&stepstop);
      printf("=== Fit time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

    }

    /*
    // write to PSCfiles
    time(&stepstart);
    agata->WritePSCfiles();
    time(&stepstop);
    printf("=== WritePSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
    */
  }
  
  time(&stop);
  printf("\n============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  if(kMakeInit || kGroupPulse) delete treereader;
  delete agata;

  return 0;
}
