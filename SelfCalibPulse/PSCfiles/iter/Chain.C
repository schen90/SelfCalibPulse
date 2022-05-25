{
  TChain *tree = new TChain();
  for(int i=0; i<36; i++){
    tree->AddFile("it0/Det0000_fit4.root",0,Form("tree%d",i));
  }
  
  /*
  TCanvas *c = new TCanvas("c","c",800,400);
  c->Divide(2,1);
  c->cd(1);
  //tree->Draw("cadpos[1]:cadpos[0]>>h1(100,-50,50,100,-50,50)","npaths>30","colz");
  tree->Draw("cadpos2[1]:cadpos2[0]>>h1(100,-50,50,100,-50,50)","npaths>30","colz");
  h1->SetStats(0);
  h1->SetTitle("HC position Y:X {npaths>30}");
  c->cd(2);
  //tree->Draw("dist>>h2(100,0,30)","npaths>30","colz");
  tree->Draw("dist2>>h2(100,0,30)","npaths>30","colz");
  h2->SetStats(0);
  h2->SetTitle("dist: calib pos - real pos {npaths>30}");
  */
  
  /*
  TCanvas *c = new TCanvas("c","c",400,800);
  c->Divide(1,2);
  c->cd(1);
  //tree->Draw("detpos[1]:detpos[0]>>h1(200,-50,50,200,-50,50)","npaths>150","colz");
  tree->Draw("cadpos2[1]:cadpos2[0]>>h1(200,-50,50,200,-50,50)","npaths>150","colz");
  c->cd(2);
  //tree->Draw("detpos[2]:detpos[0]>>h2(200,-50,50,200,-5,95)","npaths>150","colz");
  tree->Draw("cadpos2[2]:cadpos2[0]>>h2(200,-50,50,200,-5,95)","npaths>150","colz");
  */
  
  /*
  TCanvas *c = new TCanvas("c","c",400,400);
  tree->Draw("dist2>>h(100,0,30)","npaths>3");
  h->SetTitle("dist: calib_pos - real_pos {npaths>3}");
  h->Draw();
  */

  /*
  TCanvas *c = new TCanvas("c","c",400,400);
  tree->Draw("cadpos2[2]-detpos[2]>>h(100,-15,15)","npaths>3");
  h->SetTitle("calib_Z - det_Z {npaths>3}");
  h->Draw();
  */

  /*
  TCanvas *c = new TCanvas("c","c",400,800);
  c->Divide(1,2);
  c->cd(1);
  tree->Draw("dist2>>h1(100,0,30)","npaths>70");
  h1->SetTitle("dist: calib_pos - real_pos {npaths>70}");
  h1->Draw();
  c->cd(2);
  tree->Draw("cadpos2[2]-detpos[2]>>h2(100,-15,15)","npaths>70");
  h2->SetTitle("calib_Z - det_Z {npaths>70}");
  h2->Draw();
  */


  //gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c",500,800);
  c->Divide(1,2);
  c->cd(1)->SetLogz();
  tree->Draw("dist2:npaths>>h(100,0.5,300.5,100,0,30)","","colz");
  h->SetTitle("dist : npaths");
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->Draw("colz");
  //tree->Draw("cadpos2[2]-detpos[2]:npaths>>h(200,0.5,200.5,100,-15,15)","","colz");
  //h->SetTitle("calib_Z - det_Z : npaths");
  c->cd(2);
  TProfile *p = new TProfile("p","dist : npaths",100,0.5,300.5,0,30);
  p->GetYaxis()->SetRangeUser(1,11);
  p->GetXaxis()->SetLabelSize(0.05);
  p->GetYaxis()->SetLabelSize(0.05);
  tree->Draw("dist2:npaths>>p");

  
}
