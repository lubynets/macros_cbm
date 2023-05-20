// cuts a window from TH2 (allows to decrease amount of entries if they are unnecessary)
void TH2_cut(std::string fileName, std::string histoname) {
  const float lox = -10;
  const float hix = 10;
  const float loy = -1;
  const float hiy = 1.5;

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  TH2F* histoIn = fileIn->Get<TH2F>(histoname.c_str());

  if(histoIn == nullptr) throw std::runtime_error("histoIn == nullptr");

  const int bin_lox = histoIn->GetXaxis()->FindBin(lox);
  const int bin_hix = histoIn->GetXaxis()->FindBin(hix) - 1;
  const int bin_loy = histoIn->GetYaxis()->FindBin(loy);
  const int bin_hiy = histoIn->GetYaxis()->FindBin(hiy) - 1;

  TH2F* histoOut = new TH2F(histoIn->GetName(),
                            histoIn->GetTitle(),
                            bin_hix - bin_lox +1,
                            histoIn->GetXaxis()->GetBinLowEdge(bin_lox),
                            histoIn->GetXaxis()->GetBinUpEdge(bin_hix),
                            bin_hiy - bin_loy +1,
                            histoIn->GetYaxis()->GetBinLowEdge(bin_loy),
                            histoIn->GetYaxis()->GetBinUpEdge(bin_hiy)
                            );
  histoOut->GetXaxis()->SetTitle(histoIn->GetXaxis()->GetTitle());
  histoOut->GetYaxis()->SetTitle(histoIn->GetYaxis()->GetTitle());

  for(int iBinx=1; iBinx<=bin_hix - bin_lox +1; iBinx++) {
    for(int iBiny=1; iBiny<=bin_hiy - bin_loy +1; iBiny++) {
      histoOut->SetBinContent(iBinx, iBiny, histoIn->GetBinContent(bin_lox + iBinx - 1, bin_loy + iBiny -1));
    }
  }

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  histoOut->Write();
  fileOut->Close();
}
