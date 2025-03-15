void histoPrepare() {
  gROOT->Macro("/home/oleksii/alidir/macros_on_git/qa2/exe_based/styles/mc_qa2.style.cc");

  const std::vector<std::pair<Color_t, std::string>> histos {
    {kBlue, "topLeft"},
    {kRed, "topMid"},
    {kGreen+2, "topRight"},
    {kMagenta, "bottomLeft"},
    {kCyan, "bottomMid"},
    {kBlack, "bottomRight"},
  };

  int iH{0};
  for(auto& h : histos) {
    TCanvas cc("cc", "", 1200, 800);
    cc.SetRightMargin(0.03);
    TH1F histo("histo", "", 100, 0, 100);
    histo.GetXaxis()->SetTitle("X, a.u.");
    histo.GetYaxis()->SetTitle("Entries");
    histo.GetXaxis()->CenterTitle();
    histo.GetYaxis()->CenterTitle();
    const double xlo = iH<3 ? 0.001 : 0;
    const double xhi = iH == 6 ? 100 : 99.999;
    const double ylo = iH == 0 || iH == 3 ? 0 : 0.01;
    const double yhi = iH == 0 ? 5000 : 4999.99;

    histo.GetXaxis()->SetLimits(xlo, xhi);
    histo.GetYaxis()->SetRangeUser(ylo, yhi);

    histo.SetLineColor(h.first);
    for(int i=0; i<100000; i++) {
      histo.Fill(gRandom->Gaus(50, 10));
    }
    histo.Draw();
    cc.Print((h.second + ".pdf").c_str(), "pdf");
    iH++;
  }
}
