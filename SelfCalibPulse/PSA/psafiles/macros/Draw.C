

void Draw(){
  TFile *f0 = new TFile("htmp.root");
  TH1F *h0 = (TH1F *)f0->Get("hThetadiff");

  TFile *f1 = new TFile("outputSelector.root");
  TH1F *h1 = (TH1F *)f1->Get("mAngleDiff");

  h0->SetName("SimpleCode");
  h1->SetName("AGAPRO");
  h0->SetLineColor(2);
  h0->Scale(1.3);

  TCanvas *c = new TCanvas("c","c");
  h0->SetStats(1);

  h1->Draw("h");
  h0->Draw("hsame");

  TLegend *label = new TLegend(0.1,0.7,0.3,0.9);
  label->AddEntry(h0,"SimpleCode","l");
  label->AddEntry(h1,"AGAPRO","l");
  label->Draw();

}
