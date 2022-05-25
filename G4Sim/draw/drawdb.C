{
  TCanvas *c = new TCanvas("c","c",600,300);
  c->Divide(2,1);
  gStyle->SetOptStat(0);
  
  c->cd(1)->SetMargin(0.1,0.03,0.1,0.03);
  TH2D *h1 = new TH2D("h1","",100,-50,50,100,-50,50);
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetYaxis()->SetLabelSize(0.06);
  tree->Draw("pos[1]:pos[0]>>h1","","col");

  c->cd(2)->SetMargin(0.1,0.03,0.1,0.03);
  TH2D *h2 = new TH2D("h2","",100,-50,50,100,-5,95);
  h2->GetXaxis()->SetLabelSize(0.06);
  h2->GetYaxis()->SetLabelSize(0.06);
  tree->Draw("pos[2]:pos[0]>>h2","","col");

}
