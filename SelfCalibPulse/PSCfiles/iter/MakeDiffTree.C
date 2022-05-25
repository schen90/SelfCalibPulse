#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMath.h"
#include "stdlib.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

using namespace std;

void MakeDiffTree(){
  TChain *tree = new TChain();
  for(int i=0; i<36; i++){
    tree->AddFile("it0/Det0000_fit4.root",0,Form("tree%d",i));
  }
    
  int det, seg;
  float detpos[3];
  float cadpos[3];
  float cadpos2[3];
  int npaths;
  tree->SetBranchAddress("det",&det);
  tree->SetBranchAddress("seg",&seg);
  tree->SetBranchAddress("detpos",detpos);
  tree->SetBranchAddress("cadpos",cadpos);
  tree->SetBranchAddress("cadpos2",cadpos2);
  tree->SetBranchAddress("npaths",&npaths);

  int nentries = tree->GetEntriesFast();

  // output
  float phi0, r0, z0; // real
  float phi1, r1, z1; // init
  float phi2, r2, z2; // fit
  float dist1, dist2;
  float diffr1, diffz1, rdiffphi1, diffphi1;
  float diffr2, diffz2, rdiffphi2, diffphi2;
  TFile *fout = new TFile("DiffTree.root","RECREATE");
  TTree *outtree = new TTree("tree","diff tree");
  outtree->Branch("det",&det);
  outtree->Branch("seg",&seg);
  outtree->Branch("detpos",detpos,"detpos[3]/F");
  outtree->Branch("cadpos",cadpos,"cadpos[3]/F");
  outtree->Branch("cadpos2",cadpos2,"cadpos2[3]/F");
  outtree->Branch("npaths",&npaths);
  outtree->Branch("phi",&phi0);
  outtree->Branch("r",&r0);
  outtree->Branch("z",&z0);
  outtree->Branch("dist1",&dist1);
  outtree->Branch("dist2",&dist2);
  outtree->Branch("diffr1",&diffr1);
  outtree->Branch("diffz1",&diffz1);
  outtree->Branch("rdiffphi1",&rdiffphi1);
  outtree->Branch("diffr2",&diffr2);
  outtree->Branch("diffz2",&diffz2);
  outtree->Branch("rdiffphi2",&rdiffphi2);

  const int nbinsxy = 60;
  const double xymin = -50;
  const double xymax = 50;
  const int nbinsz = 60;
  const double zmin = -5;
  const double zmax = 95;

  TProfile2D *pxy[3];
  TProfile2D *pxz[3];
  TProfile3D *pxyz[3];
  string title[3] = {"diffz","diffr","rdiffphi"};
  for(int ix=0; ix<3; ix++){
    pxy[ix] = new TProfile2D(Form("pxy%d",ix),title[ix].c_str(),nbinsxy,xymin,xymax,nbinsxy,xymin,xymax);  
    pxz[ix] = new TProfile2D(Form("pxz%d",ix),title[ix].c_str(),nbinsxy,xymin,xymax,nbinsz,zmin,zmax);  
    pxyz[ix] = new TProfile3D(Form("pxyz%d",ix),title[ix].c_str(),nbinsxy,xymin,xymax,nbinsxy,xymin,xymax,nbinsz,zmin,zmax);  
  }
  
  for(int ientry=0; ientry<nentries; ientry++){
    if(ientry%1000==0) cout<<"\r finish "<<ientry<<" / "<<nentries<<"..."<<flush;
    tree->GetEntry(ientry);
    
    TVector3 rpos(detpos[0],detpos[1],0);
    TVector3 ipos(cadpos[0],cadpos[1],0);
    TVector3 fpos(cadpos2[0],cadpos2[1],0);
    phi0 = rpos.Phi()/TMath::Pi()*180;
    r0 = rpos.Mag();
    z0 = detpos[2];

    phi1 = ipos.Phi()/TMath::Pi()*180;
    r1 = ipos.Mag();
    z1 = cadpos[2];

    phi2 = fpos.Phi()/TMath::Pi()*180;
    r2 = fpos.Mag();
    z2 = cadpos2[2];

    diffr1 = r1 - r0;
    diffz1 = z1 - z0;
    diffphi1 = phi1 - phi0;
    if(diffphi1>180) diffphi1-=360;
    if(diffphi1<-180) diffphi1+=360;
    rdiffphi1 = (r1+r0)/2*diffphi1/180*TMath::Pi();

    diffr2 = r2 - r0;
    diffz2 = z2 - z0;
    diffphi2 = phi2 - phi0;
    if(diffphi2>180) diffphi2-=360;
    if(diffphi2<-180) diffphi2+=360;
    rdiffphi2 = (r2+r0)/2*diffphi2/180*TMath::Pi();

    dist1 = 0;
    dist2 = 0;
    for(int ix=0; ix<3; ix++){
      dist1 += pow(cadpos[ix]-detpos[ix],2);
      dist2 += pow(cadpos2[ix]-detpos[ix],2);
    }
    dist1 = sqrt(dist1);
    dist2 = sqrt(dist2);

    outtree->Fill();

    if(npaths>10){
      pxy[0]->Fill(detpos[0],detpos[1],fabs(diffz2));
      pxy[1]->Fill(detpos[0],detpos[1],fabs(diffr2));
      pxy[2]->Fill(detpos[0],detpos[1],fabs(rdiffphi2));

      pxz[0]->Fill(detpos[0],detpos[2],fabs(diffz2));
      pxz[1]->Fill(detpos[0],detpos[2],fabs(diffr2));
      pxz[2]->Fill(detpos[0],detpos[2],fabs(rdiffphi2));

      pxyz[0]->Fill(detpos[0],detpos[1],detpos[2],fabs(diffz2));
      pxyz[1]->Fill(detpos[0],detpos[1],detpos[2],fabs(diffr2));
      pxyz[2]->Fill(detpos[0],detpos[1],detpos[2],fabs(rdiffphi2));
    }
  }

  fout->cd();
  outtree->Write();
  for(int i=0; i<3; i++){
    pxy[i]->Write();
    pxz[i]->Write();
    pxyz[i]->Write();
  }
  fout->Close();
  
  return;
}
