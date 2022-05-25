{
  float minchi2 = 0;
  float maxchi2 = 10;
  int nbins = 100;
  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(2,2);

  c->cd(1);
  //tree->Draw(Form("chi2s[0]:dist>>h1(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1","colz");
  //tree->Draw(Form("chi2s[1]:dist>>h1(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&chi2s[0]<1&&chi2s[2]<1","colz");
  tree->Draw("dist>>h1(100,0,10)","simseg==1&&chi2s[0]<2&&chi2s[2]<1&&chi2s[1]<0.1");
  h1->SetTitle("dist");
  
  c->cd(2);
  //tree->Draw(Form("chi2s[0]:rdiffphi>>h2(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&diffr<1&&diffz<1","colz");
  //tree->Draw(Form("chi2s[1]:rdiffphi>>h2(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&chi2s[0]<1&&chi2s[2]<1","colz");
  tree->Draw("rdiffphi>>h2(100,0,10)","simseg==1&&chi2s[0]<2&&chi2s[2]<1&&chi2s[1]<0.1");
  h2->SetTitle("rdiffphi");
  
  c->cd(3);
  //tree->Draw(Form("chi2s[0]:diffr>>h3(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&rdiffphi<1&&diffz<1","colz");
  //tree->Draw(Form("chi2s[1]:diffr>>h3(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&chi2s[0]<1&&chi2s[2]<1","colz");
  tree->Draw("diffr>>h3(100,0,10)","simseg==1&&chi2s[0]<2&&chi2s[2]<1&&chi2s[1]<0.1");
  h3->SetTitle("diffr");

  c->cd(4);
  //tree->Draw(Form("chi2s[0]:diffz>>h4(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&diffr<1&&rdiffphi<1","colz");
  //tree->Draw(Form("chi2s[1]:diffz>>h4(%d,0,40,%d,%f,%f)",nbins,nbins,minchi2,maxchi2),"simseg==1&&chi2s[0]<1&&chi2s[2]<1","colz");
  tree->Draw("diffz>>h4(100,0,10)","simseg==1&&chi2s[0]<2&&chi2s[2]<1&&chi2s[1]<0.1");
  h4->SetTitle("diffz");

}
