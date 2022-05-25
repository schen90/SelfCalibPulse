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
#ifdef ADDPS
    cout<<"AddPS; ";
#endif
    if(kPSA) cout<<"PSA; ";
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
  agata->SetWithPS(true);
  agata->SetGroupPos(false);

#ifdef NOISE
  double chi2slimit[3] = {2,0.3,2}; // group pulse shape with chi2s[3] w/ noise
  agata->SetMaxChi2s(chi2slimit[0],chi2slimit[1],chi2slimit[2]);
  cout<<Form("With noise Initial chi2s limits: %.1f  %.1f  %.1f",chi2slimit[0],chi2slimit[1],chi2slimit[2])<<endl;
  cout<<Form("PSCEmin = %.0f keV",PSCEMIN)<<endl;
#else
  double chi2slimit[3] = {2,0.1,1}; // group pulse shape with chi2s[3] w/o noise
  agata->SetMaxChi2s(chi2slimit[0],chi2slimit[1],chi2slimit[2]);
  cout<<Form("Without noise, Initial chi2s limits: %.1f  %.1f  %.1f",chi2slimit[0],chi2slimit[1],chi2slimit[2])<<endl;
#endif

  

  //*****************************************//
  // Generate HitCollections
  //*****************************************//
  if(kGroupPulse){

    time(&stepstart);
    if(kConfig) treereader->GenerateHCs(agata);
    time(&stepstop);
    printf("=== GenerateHCs time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));

    // remove PSCs with nhits<MINHITS, and divide PSCs with nhits>MAXHITS
    time(&stepstart);
    agata->CheckPSCs(MINHITS,MAXHITS);
    time(&stepstop);
    printf("=== CheckHCs time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
    
    // loop till no new PSC added
    treereader->SetUpdateHCs(2); //2: divide HCs
    agata->SetAddNewPSC(true);
    for(int iupdate=0; iupdate<3; iupdate++){
      bool kNewPSC = agata->IsNewPSC(); // if new PSC added
      if(!kNewPSC) break;
      cout<<"\e[1;31m Update Divided HCs "<<iupdate<<" ... \e[0m"<<endl;
      treereader->SetNewPSC(kNewPSC);
      agata->SetAddNewPSC(false);

      time(&stepstart);
      if(kConfig) treereader->GenerateHCs(agata);
      time(&stepstop);
      printf("=== Update DivideHCs %d time: %.1f seconds ===\n\n",iupdate,difftime(stepstop,stepstart));
    }
    treereader->SetNewPSC(false);
    agata->SetAddNewPSC(false);
    agata->ClearHitLevelMarker(2);

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
      agata->Tracking();
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

    // write to PSCfiles
    time(&stepstart);
    agata->WritePSCfiles();
    time(&stepstop);
    printf("=== WritePSC time: %.1f seconds ===\n\n",difftime(stepstop,stepstart));
  }
  
  time(&stop);
  printf("\n============ Elapsed time: %.1f seconds =============\n",difftime(stop,start));

  if(kMakeInit || kGroupPulse) delete treereader;
  delete agata;

  return 0;
}