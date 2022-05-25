{
  TFile *f = new TFile("DiffTree.root");
  TTree *tree = (TTree *)f->Get("tree");
  
  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(2,2);

  c->cd(1);
  tree->Draw("dist1>>h10(100,0,30)","npaths>150");
  tree->Draw("dist2>>h1(100,0,30)","npaths>150");
  h1->SetTitle("dist {npaths>150}");
  h10->SetLineColor(17);
  h1->Draw();
  h10->Draw("same");
  h1->Draw("same");

  c->cd(2);
  tree->Draw("diffz1>>h20(100,-15,15)","npaths>150");
  tree->Draw("diffz2>>h2(100,-15,15)","npaths>150");
  h2->SetTitle("diffz {npaths>150}");
  h20->SetLineColor(17);
  h2->Draw();
  h20->Draw("same");
  h2->Draw("same");

  c->cd(3);
  tree->Draw("diffr1>>h30(100,-15,15)","npaths>150");
  tree->Draw("diffr2>>h3(100,-15,15)","npaths>150");
  h3->SetTitle("diffr {npaths>150}");
  h30->SetLineColor(17);
  h3->Draw();
  h30->Draw("same");
  h3->Draw("same");

  c->cd(4);
  tree->Draw("rdiffphi1>>h40(100,-15,15)","npaths>150");
  tree->Draw("rdiffphi2>>h4(100,-15,15)","npaths>150");
  h4->SetTitle("rdiffphi {npaths>150}");
  h40->SetLineColor(17);
  h4->Draw();
  h40->Draw("same");
  h4->Draw("same");


}
