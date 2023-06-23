#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <TMatrixD.h>
#include <TVector3.h>
#include <vector>
#include "TMath.h"
#include <algorithm>

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;

const float mec2_kev = 511;

#define MaxNDets 180

// transform matrix
TMatrixD Rt[MaxNDets];
TMatrixD Rt2[MaxNDets];
TMatrixD Tr[MaxNDets];

Int_t NDets; // number of detectors from LookUpTable

void LoadMatrix(string LookUpTable = "CrystalPositionLookUpTable"){
  // find input LookUpTable
  ifstream fin;
  int dummy_i;
  fin.open(LookUpTable.c_str());
  if(!fin){ cerr<<"Cannot find "<<LookUpTable<<endl; return;}
  cout<<"\e[1;32m find CrystalPosition from "<<LookUpTable<<"... \e[0m"<<endl;

  int ir;
  for(int i=0; i<MaxNDets; i++){
    ir = -1;
    fin >> ir >> dummy_i;
    if(ir<0 || ir>=MaxNDets){ ir=i-1; break;}

    Tr[ir].ResizeTo(3,1); Tr[ir].Zero();
    Rt[ir].ResizeTo(3,3); Rt[ir].Zero();
    Rt2[ir].ResizeTo(3,3); Rt2[ir].Zero();
    for(int it=0; it<3; it++) fin >> Tr[ir](it,0);
    for(int it=0; it<3; it++){
      fin >> dummy_i;
      for(int it2=0; it2<3; it2++) fin >> Rt[ir](it,it2);
      for(int it2=0; it2<3; it2++) Rt2[ir](it,it2) =  Rt[ir](it,it2);
    }
    Rt[ir].Invert();  // change to rot from world frame -> detector frame
  }
  fin.close();
  NDets = ir+1;
  cout<<"read position for "<<NDets<<" detectors"<<endl;

  return;
}


TMatrixD Lab2DetPos(int detid, TMatrixD LabPos){
  if(detid<0 || detid>=NDets){
    cerr<<"cannot find matrix for detid = "<<detid<<endl;
    return LabPos;
  }
  TMatrixD DetPos(3,1);
  DetPos = Rt[detid]*(LabPos-Tr[detid]);
  return DetPos;
}

TMatrixD Det2LabPos(Int_t detid, TMatrixD DetPos){
  if(detid<0 || detid>=NDets){
    cerr<<"cannot find matrix for detid = "<<detid<<endl;
    return DetPos;
  }
  TMatrixD LabPos(3,1);
  LabPos = Rt2[detid]*DetPos+Tr[detid];
  return LabPos;
}


Double_t ComptonEsca(Double_t *x, Double_t *par){
  Double_t Einc = par[0];
  Double_t theta = x[0]/180.*TMath::Pi();

  Double_t costheta = TMath::Cos(theta);
  Double_t tmpa = 1 + mec2_kev/Einc - costheta;
  Double_t Esca = mec2_kev / tmpa;
  return Esca;
}


Double_t ComptonTheta(Double_t *x, Double_t *par){
  Double_t Einc = par[0];
  Double_t Esca = x[0];

  Double_t costheta = 1 + mec2_kev/Einc - mec2_kev/Esca;
  Double_t theta = TMath::ACos(costheta)/TMath::Pi()*180.;
  return theta;
}


void Compton(){
  LoadMatrix("CrystalPositionLookUpTable");

  struct OBJ{
    int   EntryID;
    int   CrystalId;
    float CoreE[2];
    int   numNetCharges;
    int   seg;
    float NetCh;
    float detpos[3];
    float labpos[3];
  };
  OBJ obj;

  TChain *tree = new TChain();
  for(int itree=0; itree<54; itree++) 
    tree->AddFile(Form("run_0008/Tree_%04d.root",itree),0,"tree");

  /*
  for(int itree=0; itree<17; itree++) 
    tree->AddFile(Form("run_0009/Tree_%04d.root",itree),0,"tree");

  for(int itree=0; itree<11; itree++) 
    tree->AddFile(Form("run_0010/Tree_%04d.root",itree),0,"tree");
  */

  //for(int itree=0; itree<20; itree++) 
  //  tree->AddFile(Form("run_0016/Tree_%04d.root",itree),0,"tree");

  tree->SetBranchAddress("EntryID",       &obj.EntryID);
  tree->SetBranchAddress("CrystalId",     &obj.CrystalId);
  tree->SetBranchAddress("CoreE",          obj.CoreE);
  tree->SetBranchAddress("numNetCharges", &obj.numNetCharges);
  tree->SetBranchAddress("seg",           &obj.seg);
  tree->SetBranchAddress("NetCh",         &obj.NetCh);
  //tree->SetBranchAddress("segpos",        &obj.detpos);
  tree->SetBranchAddress("pos",           &obj.detpos);

  int nentries = tree->GetEntriesFast();
  cout<<"find "<<nentries<<" entries..."<<endl;

  // Compton scattering functions
  TF1 *fEsca = new TF1("fEsca",ComptonEsca,0,180,1);
  fEsca->SetParameter(0,1274.);

  TF1 *fTheta = new TF1("fTheta",ComptonTheta,0,1500,1);
  fTheta->SetParameter(0,1274.);

  // histogram
  TH1F *hEtot = new TH1F("hEtot","total Energy",2000,0,2000);
  TH2F *hEscaTheta = new TH2F("hEscaTheta","Scattering angle vs. Scattering Energy",180,-0.5,179.5,1024,-0.5,1023.5);
  TH1F *hThetadiff = new TH1F("hThetadiff","Theta PSA - Theta Compton",200,-100,100);

  vector<OBJ> fhits;
  for(int ientry=0; ientry<nentries; ientry++){
    if(ientry%10000==0) cout<<"\r finish "<<ientry<<" / "<<nentries<<"..."<<flush;

    tree->GetEntry(ientry);

    TMatrixD DetPos(3,1);
    for(int ix=0; ix<3; ix++) DetPos(ix,0) = obj.detpos[ix];
    TMatrixD LabPos = Det2LabPos(obj.CrystalId, DetPos);
    for(int ix=0; ix<3; ix++) obj.labpos[ix] = LabPos(ix,0);

    if(fhits.size()>0){
      if(obj.EntryID!=fhits[0].EntryID){ 

	// sort
	sort(fhits.begin(), fhits.end(), 
	     [](const OBJ lhs, const OBJ rhs){
	       return lhs.NetCh>rhs.NetCh;});

	// calculate Compton scattering 
	// total energy of the array
	float Etot = 0;
	for(int i=0; i<fhits.size(); i++) Etot += fhits[i].NetCh;
	hEtot->Fill(Etot);

	if(fabs(Etot-1274)<5 && fhits.size()==2){
	  TVector3 v0(fhits[0].labpos[0],fhits[0].labpos[1],fhits[0].labpos[2]);
	  TVector3 v1(fhits[1].labpos[0],fhits[1].labpos[1],fhits[1].labpos[2]);

	  TVector3 vsca01 = v1-v0;
	  TVector3 vsca10 = v0-v1;
	  float theta01 = v0.Angle(vsca01)/TMath::Pi()*180;
	  float theta10 = v1.Angle(vsca10)/TMath::Pi()*180;

	  float theta01diff = theta01 - fTheta->Eval(fhits[1].NetCh);
	  float theta10diff = theta10 - fTheta->Eval(fhits[0].NetCh);


	  if(theta01diff<theta10diff){
	    hEscaTheta->Fill(theta01,fhits[1].NetCh);
	    hThetadiff->Fill(theta01diff);
	  }else{
	    hEscaTheta->Fill(theta10,fhits[0].NetCh);
	    hThetadiff->Fill(theta10diff);
	  }


	  /*
	  hEscaTheta->Fill(theta01,fhits[1].NetCh);
	  hThetadiff->Fill(theta01diff);
	  */
	}

	fhits.clear();
      }
    }

    fhits.push_back(obj);
  }

  TCanvas *c = new TCanvas("c","c");
  hEscaTheta->GetXaxis()->SetRangeUser(0,180);
  hEscaTheta->GetYaxis()->SetRangeUser(125,1150);
  hEscaTheta->DrawClone("colz");
  fEsca->SetLineWidth(1);
  fEsca->Draw("same");

  TFile *fout = new TFile("htmp.root","RECREATE");
  hEtot->Write();
  hEscaTheta->Write();
  hThetadiff->Write();
  fout->Close();
  return;
}
