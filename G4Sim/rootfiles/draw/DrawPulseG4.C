#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVector3.h"
#include <vector>
#include "TInterpreter.h"

TGraph *gr;
TGraph *grgrid[9];

void DrawPulseG4(){

  gInterpreter->GenerateDictionary("vector<vector<float>>","vector");
  gInterpreter->GenerateDictionary("vector<vector<int>>","vector");

  //TCanvas *c = new TCanvas("c","c",1265,670);
  TCanvas *c1 = new TCanvas("c1","c1",1450,780);
  c1->Divide(3,3);
  c1->cd(1);
  //TFile *f = new TFile("rootfiles/extra/G4SimData0000.root");
  TFile *f = new TFile("rootfiles/PS/G4SimData0000.root");
  TTree *tree = (TTree *)f->Get("tree");

  int ievent;
  vector<int> *ndet = 0;
  vector<int> *g4seg = 0;
  vector<float> *energy = 0;
  vector<vector<float>> *posa = 0;
  vector<vector<float>> *posr = 0;

  vector<int> *pdet = 0;
  vector<float> *ecore = 0;
  vector<vector<int>> *inter = 0;

  vector<vector<int>> *pseg = 0;
  vector<vector<int>> *ngrid = 0;
  vector<vector<int>> *extrpl = 0;
  vector<vector<float>> *core = 0;
  vector<vector<float>> *spulse = 0;
  vector<vector<int>> *gridip = 0;
  vector<vector<float>> *gridwgt = 0;

  tree->SetBranchAddress("ievent",&ievent);
  tree->SetBranchAddress("ndet",&ndet);
  tree->SetBranchAddress("g4seg",&g4seg);
  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("posa",&posa);
  tree->SetBranchAddress("posr",&posr);

  tree->SetBranchAddress("pdet",&pdet);
  tree->SetBranchAddress("ecore",&ecore);
  tree->SetBranchAddress("inter",&inter);

  tree->SetBranchAddress("pseg",&pseg);
  tree->SetBranchAddress("ngrid",&ngrid);
  tree->SetBranchAddress("extrpl",&extrpl);
  tree->SetBranchAddress("core",&core);
  tree->SetBranchAddress("spulse",&spulse);
  tree->SetBranchAddress("gridip",&gridip);
  tree->SetBranchAddress("gridwgt",&gridwgt);

  int nentries = tree->GetEntriesFast();

  int simdet;
  int simseg;
  float simeng;
  float simdetpos[3];
  int simngrid = 0;
  int simextrpl;
  vector<int> simgridip;
  vector<float> simgridwgt;
  float tmpspulse[4477];
  int iety;
  int idet;
  int nselect = 1;
  int iselect = 0;
  double xmin = -0.5, xmax = 36.5;
  //double xmin = 17.5, xmax = 22.5;
  for(iety=0; iety<nentries; iety++){
    tree->GetEntry(iety);

    for(idet=0; idet<pdet->size(); idet++){
      if(inter->at(idet).size()>1) continue; // use pulse with only one interaction in a detector

      simdet = pdet->at(idet);
      int interid = inter->at(idet)[0];
      simseg = pseg->at(idet)[0]-1;
      simeng = energy->at(interid);
      for(int iaxis=0; iaxis<3; iaxis++) simdetpos[iaxis] = posr->at(interid)[iaxis];
      simngrid = ngrid->at(idet)[0];
      simextrpl = extrpl->at(idet)[0];
      for(int isig=0; isig<4356; isig++) tmpspulse[isig] = spulse->at(idet)[isig];
      for(int isig=0; isig<121; isig++) tmpspulse[isig+4356] = core->at(idet)[isig];
      simgridip.clear();
      simgridwgt.clear();
      for(int igrid=0; igrid<gridip->at(idet).size(); igrid++){
	simgridip.push_back(gridip->at(idet)[igrid]);
	simgridwgt.push_back(gridwgt->at(idet)[igrid]);
      }

      //if(simngrid==8) iselect++;
      //if(simngrid>8) iselect++;
      //if(simngrid<8) iselect++;
      //if(simextrpl==0) iselect++;
      if(simextrpl==7) iselect++;
      //if(simextrpl>100) iselect++;

      if(iselect==nselect) break;
    }

    if(iselect==nselect) break;
  }
  
  float x[4477];
  for(int iseg=0; iseg<37; iseg++)
    for(int i=0; i<121; i++)
      x[iseg*121+i] = iseg+i*1./121-0.5;

  gr = new TGraph(4477,x,tmpspulse);
  gr->SetTitle(Form("pos: %.2f, %.2f, %.2f",simdetpos[0],simdetpos[1],simdetpos[2]));

  gr->GetXaxis()->SetRangeUser(xmin,xmax);
  //gr->GetYaxis()->SetRangeUser(-0.1,1.1);
  gr->Draw("APL");

  TVector3 ivec(simdetpos[0],simdetpos[1],0);
  float simphi = ivec.Phi()/TMath::Pi()*180;
  float simr = ivec.Mag();
  float simz = simdetpos[2];
  
  cout<<"entry : "<<iety<<" idet: "<<idet<<endl
      <<"det   : "<<simdet<<endl
      <<"seg   : "<<simseg<<endl
      <<"eng   : "<<simeng<<endl
      <<"posXYZ: "<<Form("%.2f  %.2f  %.2f",simdetpos[0],simdetpos[1],simdetpos[2])<<endl
      <<"PhiRZ : "<<Form("%.2f  %.2f  %.2f",simphi,simr,simz)<<endl
      <<"ngrid : "<<simngrid<<endl
      <<"extr  : "<<simextrpl<<endl
      <<"grid  : ";
  for(int ig=0; ig<simgridip.size(); ig++) cout<<simgridip[ig]<<" ";
  cout<<endl;

  string dbfile[3] = {"pulsedb/pulseA.root","pulsedb/pulseB.root","pulsedb/pulseC.root"};
  int itype = simdet%3;
  TFile *fdb = new TFile(dbfile[itype].c_str());
  if(!fdb->IsOpen()){
    cerr<<"cannot fine dbfile "<<dbfile[itype]<<endl;
    return;
  }
  TTree *dbtree = (TTree *)fdb->Get("tree");
  Int_t dbseg;
  Double_t dbpos[3];
  Double_t dbcore[121];
  Double_t dbspulse[4356];
  dbtree->SetBranchAddress("seg",&dbseg);
  dbtree->SetBranchAddress("pos",dbpos);
  dbtree->SetBranchAddress("core",dbcore);
  dbtree->SetBranchAddress("spulse",dbspulse);
  int npoint = dbtree->GetEntriesFast();
  
  for(int ig=0; ig<simgridip.size() && ig<8; ig++){
    c1->cd(ig+2);
    int gip = simgridip[ig];
    float gwgt = simgridwgt[ig];
    dbtree->GetEntry(gip);
    
    for(int isig=0; isig<4356; isig++) tmpspulse[isig] = dbspulse[isig];
    for(int isig=0; isig<121; isig++) tmpspulse[isig+4356] = dbcore[isig];

    float dist = 0;
    for(int iaxis=0; iaxis<3; iaxis++) dist+=pow(simdetpos[iaxis]-dbpos[iaxis],2);
    dist = sqrt(dist);
    grgrid[ig] = new TGraph(4477,x,tmpspulse);
    grgrid[ig]->SetTitle(Form("grid %d,  wgt %.3f,  pos: %.2f, %.2f, %.2f, dist: %.2f",gip,gwgt,dbpos[0],dbpos[1],dbpos[2],dist));
    grgrid[ig]->GetXaxis()->SetRangeUser(xmin,xmax);
    //grgrid[ig]->GetYaxis()->SetRangeUser(-0.1,1.1);
    grgrid[ig]->Draw("APL");
  }

  fdb->Close();
  f->Close();
  return;
}
