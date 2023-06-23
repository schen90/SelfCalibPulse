

void Draw(){
  TFile *f0 = new TFile("htmp.root");
  TH1F *h0 = (TH1F *)f0->Get("hThetadiff");

  TFile *f1 = new TFile("htmp_calib.root");
  TH1F *h1 = (TH1F *)f1->Get("hThetadiff");

  h0->SetName("Simulation");
  h1->SetName("Calibration");
  h0->SetLineColor(2);

  h1->Scale(h0->GetSum()/h1->GetSum());

  TCanvas *c = new TCanvas("c","c");
  h0->SetStats(1);
  h0->SetTitle("Angle PSA - Angle Compton");

  h0->Draw("h");
  h1->Draw("hsame");

  TLegend *label = new TLegend(0.1,0.7,0.3,0.9);
  label->AddEntry(h0,"Simulation","l");
  label->AddEntry(h1,"Calibration","l");
  label->Draw();

}
