{
  float minchi2 = 0;
  float maxchi2 = 10;
  int nbins = 100;
  TCanvas *c = new TCanvas("c","c",1200,700);
  c->Divide(4,3);
  int segs[4] = {0,1,2,5};
  //int segs[4] = {0,3,4,5};
  //int segs[4] = {0,6,12,18};

  TH2D *hr[4];
  TH2D *hz[4];
  TH2D *hphi[4];
  for(int i=0; i<4; i++){
    c->cd(1+i);
    hr[i] = new TH2D(Form("hr%d",i+1),Form("chi2s[1]:diffr{seg%d}",segs[i]),nbins,0,40,nbins,minchi2,maxchi2);
    tree->Draw(Form("chi2s[1]:diffr>>hr%d",i+1),Form("simseg==%d&&rdiffphi<1&&diffz<1&&chi2s[0]<2&&chi2s[2]<2",segs[i]),"colz");
    c->cd(5+i);
    hz[i] = new TH2D(Form("hz%d",i+1),Form("chi2s[1]:diffz{seg%d}",segs[i]),nbins,0,40,nbins,minchi2,maxchi2);
    tree->Draw(Form("chi2s[1]:diffz>>hz%d",i+1),Form("simseg==%d&&diffr<1&&rdiffphi<1&&chi2s[0]<2&&chi2s[2]<2",segs[i]),"colz");
    c->cd(9+i);
    hphi[i] = new TH2D(Form("hphi%d",i+1),Form("chi2s[1]:rdiffphi{seg%d}",segs[i]),nbins,0,40,nbins,minchi2,maxchi2);
    tree->Draw(Form("chi2s[1]:rdiffphi>>hphi%d",i+1),Form("simseg==%d&&diffr<1&&diffz<1&&chi2s[0]<2&&chi2s[2]<2",segs[i]),"colz");
  }
  

}
