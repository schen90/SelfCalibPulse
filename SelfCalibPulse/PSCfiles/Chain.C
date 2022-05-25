{
  TChain *tree = new TChain();
  for(int i=0; i<36; i++)
    tree->AddFile("./Det0000.root",0,Form("tree%d",i));


}
