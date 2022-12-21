void Q_parper() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/covariances_scol.lambda.root";

  MultiCorrelation multicor;
  multicor.SetPalette({kBlue, kRed});
  multicor.SetMarkers({kFullCircle, kFullCircle});

  multicor.AddCorrelation(fileName, {"Qperp.all"}, "Qperp");

  HeapPicture pic("picture", {1000, 1000});
  pic.AddText({0.2, 0.89, "Au+Au"}, 0.025);
  pic.AddText({0.2, 0.86, "DCM-QGSM-SMM"}, 0.025);
  pic.AddText({0.2, 0.83, "12A GeV/c"}, 0.025);

//   auto leg1 = new TLegend();

  for(auto& obj : multicor.GetCorrelations()) {
    obj->SetCalculateSystematicsFromVariation(false);
    pic.AddDrawable(obj);
  }

  pic.SetAxisTitles( {"b, fm", "#LT Q_{#perp} #GT"} );
  pic.CustomizeXRange();
  pic.CustomizeYRange();
  pic.SetAutoLegend(false);
//   pic.AddLegend(leg1);
//   pic.CustomizeLegend(leg1);
  pic.Draw();

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  fileOut->cd();
  pic.GetCanvas()->Write();
  fileOut->Close();

  pic.GetCanvas()->Print("fileOut.pdf", "pdf");
}
