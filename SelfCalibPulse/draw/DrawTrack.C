{
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c",400,600);
  c->Divide(1,2);

  c->cd(1);
  TH1D *h1 = new TH1D("h1","nhits",10,-0.5,9.5);
  tree->Draw("nhits>>h1","","");

  c->cd(2);
  THStack *h2 = new THStack("h2","correctness");
  
  TH1D *h21 = new TH1D("h21","source at 0",4,-1.5,2.5);
  tree->Draw("correct>>h21","sourcepos[2]>-1 && sourcepos[2]<1","");
  h21->SetFillColor(2);
  h21->SetFillStyle(3004);
  h2->Add(h21);

  TH1D *h22 = new TH1D("h22","source at -160",4,-1.5,2.5);
  tree->Draw("correct>>h22","sourcepos[2]>-161 && sourcepos[2]<-159","");
  h22->SetFillColor(3);
  h22->SetFillStyle(3004);
  h2->Add(h22);

  TH1D *h23 = new TH1D("h23","source at -400",4,-1.5,2.5);
  tree->Draw("correct>>h23","sourcepos[2]>-401 && sourcepos[2]<-399","");
  h23->SetFillColor(4);
  h23->SetFillStyle(3004);
  h2->Add(h23);

  h2->Draw();

  double bin0 = h21->GetBinContent(2) + h22->GetBinContent(2) + h23->GetBinContent(2);
  double bin1 = h21->GetBinContent(3) + h22->GetBinContent(3) + h23->GetBinContent(3);
  double total = bin0+bin1;
  double ratio0 = bin0/total*100;
  double ratio1 = bin1/total*100;
  cout<<"bin0 = "<<bin0<<endl;
  cout<<"bin1 = "<<bin1<<endl;
  cout<<"ratio1 = "<<ratio1<<endl;
  
  TLegend *label = new TLegend(0.15,0.65,0.4,0.85);
  label->AddEntry(h23,"source at -400","f");
  label->AddEntry(h22,"source at -160","f");
  label->AddEntry(h21,"source at 0","f");
  label->Draw("same");
}
