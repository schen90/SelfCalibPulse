#include "TStyle.h"
#include "TList.h"
#include "TRandom.h"
#include "TH3.h"
#include "TProfile3D.h"

void DrawPulseDiff3D(){
  // get hist 3
  TFile *f = new TFile("PulseDiffTree.root");
  TProfile3D *h3 = (TProfile3D *)f->Get("pxyz");
  h3->Fill(0.,0.,0.,20.);

  TList *lf = h3->GetListOfFunctions();
  if (lf){
    //TF1 *tf = new TF1("TransferFunction","x*x/100",0,10.);
    TF1 *tf = new TF1("TransferFunction","x*100/(x*x*x)",0.,10.);
    lf->Add(tf);
  }
  gStyle->SetCanvasPreferGL(1);
  gStyle->SetPalette(1);
  TCanvas *c = new TCanvas("c","abc",700,700);
  h3->Draw("glcolz");

  
}
