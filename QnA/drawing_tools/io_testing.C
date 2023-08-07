void io_testing() {
  Picture hp("pic", {1000, 1000});
//   Picture hp;
//   TCanvas hp("pic", "", 1000, 1000);
//   Correlation hp;

  TFile* fileOut = new TFile("fileOut.root", "recreate");
  hp.Write("mypic");
  fileOut->Close();
}
