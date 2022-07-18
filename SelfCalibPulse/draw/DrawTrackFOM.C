{
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c","c",1000,600);
  c->Divide(3,2);

  char spos0[3][100] = {"source at 0",
			"source at -160",
			"source at -400"};

  char spos[3][100] = {"sourcepos[2]>-1 && sourcepos[2]<1",
		       "sourcepos[2]>-161 && sourcepos[2]<-159",
		       "sourcepos[2]>-401 && sourcepos[2]<-399"};

  TH1D *h1[3][2];
  TH1D *h2[3][2];
  
  for(int i=0; i<3; i++){
    c->cd(i+1)->SetLogy();
    h1[i][0] = new TH1D(Form("h1_%d_0",i),Form("Log(FOM1/FOM2) {%s}",spos0[i]),100,-100,0);
    h1[i][1] = new TH1D(Form("h1_%d_1",i),Form("Log(FOM1/FOM2) {%s}",spos0[i]),100,-100,0);
    h1[i][0]->SetLineColor(2);
    tree->Draw(Form("log(FOM1/FOM2)>>h1_%d_0",i),Form("%s && correct==0",spos[i]),"");
    tree->Draw(Form("log(FOM1/FOM2)>>h1_%d_1",i),Form("%s && correct==1",spos[i]),"");
    h1[i][1]->Draw();
    h1[i][0]->Draw("same");

    c->cd(i+4)->SetLogy();
    h2[i][0] = new TH1D(Form("h2_%d_0",i),Form("Log(FOM1) {%s}",spos0[i]),100,0,100);
    h2[i][1] = new TH1D(Form("h2_%d_1",i),Form("Log(FOM1) {%s}",spos0[i]),100,0,100);
    h2[i][0]->SetLineColor(2);
    tree->Draw(Form("log(FOM1)>>h2_%d_0",i),Form("%s && correct==0",spos[i]),"");
    tree->Draw(Form("log(FOM1)>>h2_%d_1",i),Form("%s && correct==1",spos[i]),"");
    h2[i][1]->Draw();
    h2[i][0]->Draw("same");
  }

}
