void Pdfer(const std::string&& fileName) {
  TFile* file_fit = TFile::Open(fileName.c_str(), "read");

  bool is_first_canvas = true;

  std::string dirName = "Fits";

  TDirectory* dir = file_fit->Get<TDirectory>(dirName.c_str());
  TList* lok = dir->GetListOfKeys();

  for(auto&& k : *lok) {
    std::string canvasname = k->GetName();
    TCanvas* cc = file_fit -> Get<TCanvas>((dirName + "/" + canvasname).c_str());

    if(is_first_canvas)
      cc->Print("fileOut.pdf(", "pdf");
    else
      cc->Print("fileOut.pdf", "pdf");

    is_first_canvas = false;
  }

  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("fileOut.pdf)", "pdf");
}
