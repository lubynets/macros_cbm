void resolution_multipad() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style_multipad.cc" );

  std::vector<std::string> method{"sp", "ep"};
  std::vector<std::string> setup{"u12", "d12", "d3"};

  MultiPicture* mpic = new MultiPicture(3, 2);
  for(int i=0; i<3; i++) {
    for(int j=0; j<2; j++) {
      TFile* fileIn = TFile::Open(("/home/oleksii/cbmdir/flow_drawing_tools/input/res.psd." + method.at(j) + "." + setup.at(i) + ".root").c_str());
      HeapPicture* hp = fileIn->Get<HeapPicture>("heap_picture");
      mpic->SetPicture(i, j, hp);
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
  mpic->GetCanvas()->Draw();
  mpic->GetCanvas()->SaveAs("cc.pdf");
  mpic->GetCanvas()->SaveAs("cc.root");

}
