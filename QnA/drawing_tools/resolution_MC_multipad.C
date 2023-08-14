void resolution_MC_multipad() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style_multipad.cc" );

  std::vector<std::string> method{"sp", "ep"};
  std::vector<std::string> setup{"u12", "d12", "d3"};

  std::array<std::array<HeapPicture*, 2>, 3> hp;

  std::string filePath = "/home/oleksii/cbmdir/working/qna/aXmass/root/resolutions/psi_rp/";

  MultiPicture* mpic = new MultiPicture(3, 2);
  for(int i=0; i<3; i++) {
    for(int j=0; j<2; j++) {
      TFile* fileIn = TFile::Open((filePath + "res.psd." + method.at(j) + "." + setup.at(i) + ".root").c_str());
      hp.at(i).at(j) = fileIn->Get<HeapPicture>("heap_picture");
      mpic->SetPicture(i, j, hp.at(i).at(j));
    }
  }

  mpic->SetAspectRatio(1500, 1000);
  mpic->SetLeftMargin(0.05);
  mpic->SetRightMargin(0.01);
  mpic->SetBottomMargin(0.05);
  mpic->SetTopMargin(0.01);
  mpic->SetIsEqualizeXaxisRanges();
  mpic->SetIsEqualizeYaxisRanges();

  mpic->Draw();

  mpic->GetCanvas()->SaveAs("cc.pdf");
  mpic->GetCanvas()->SaveAs("cc.root");

}
