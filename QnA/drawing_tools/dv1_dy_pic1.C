void dv1_dy_pic1() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style_multipad.cc" );
  gStyle->SetTitleSize(36, "XY");
  gStyle->SetLabelSize(36, "XY");
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.6, "Y");

  std::vector<std::string> pbeam{"12", "3.3"};

  std::array<std::array<HeapPicture*, 1>, 2> hp;

  std::string filePath = "/home/oleksii/cbmdir/working/qna/aXmass/root/dv1dy/";

  MultiPicture* mpic = new MultiPicture(2, 1);
  for(int i=0; i<2; i++) {
    TFile* fileIn = TFile::Open((filePath + "slope.dcmqgsm." + pbeam.at(i) + "agev.3122.root").c_str());
    hp.at(i).at(0) = fileIn->Get<HeapPicture>("heap_picture_uQ_R1_sub4");
    mpic->SetPicture(i, 0, hp.at(i).at(0));
  }

  mpic->SetAspectRatio(1800, 1000);
  mpic->SetLeftMargin(0.06);
  mpic->SetRightMargin(0.01);
  mpic->SetBottomMargin(0.10);
  mpic->SetTopMargin(0.01);
  mpic->SetIsEqualizeXaxisRanges();
  mpic->SetIsEqualizeYaxisRanges();
//   mpic->GetPicture(1, 0)->ClearTexts({0, 2, 4});

  mpic->Draw();

  mpic->GetCanvas()->SaveAs("cc.png");
  mpic->GetCanvas()->SaveAs("cc.pdf");
  mpic->GetCanvas()->SaveAs("cc.root");
  mpic->GetCanvas()->SaveAs("cc.C");
}
