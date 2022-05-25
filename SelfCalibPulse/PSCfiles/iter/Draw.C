{
  TChain *tree = new TChain();
  tree->AddFile("PulseDiffTree.root",0,"tree");

  TCanvas *c = new TCanvas("c","c",1000,800);
  c->Divide(2,2);
  c->cd(1)->SetLogz();
  tree->Draw("chi2:npaths>>h1(100,0.5,300.5,100,0,150)","seg==dbseg && extrpl==0","colz");
  h1->SetTitle("chi2 : npaths {extrpl==0}");
  h1->Draw("colz");

  c->cd(3);
  TProfile *p1 = new TProfile("p1","chi2 : npaths {extrpl==0}",100,0.5,300.5,0,3000);
  tree->Draw("chi2:npaths>>p1","seg==dbseg && extrpl==0");
  p1->GetYaxis()->SetRangeUser(0,80);

  c->cd(2)->SetLogz();
  tree->Draw("chi2:npaths>>h2(100,0.5,300.5,100,0,150)","seg==dbseg && extrpl>0","colz");
  h2->SetTitle("chi2 : npaths {extrpl>0}");
  h2->Draw("colz");

  c->cd(4);
  TProfile *p2 = new TProfile("p2","chi2 : npaths {extrpl>0}",100,0.5,300.5,0,3000);
  tree->Draw("chi2:npaths>>p2","seg==dbseg && extrpl>0");
  p2->GetYaxis()->SetRangeUser(0,80);

  
}
