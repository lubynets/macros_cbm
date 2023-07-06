void resolution_presubevent4() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3";
  std::string evegen = "urqmd";  std::string pbeam = "12";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/" + pbeam + "agev/v1andR1.stf." + evegen + "." + pbeam + "agev.root";

  std::vector<std::string> correls{"psd1", "psd2", "psd3"};

  std::string sub4th = "sts_pipos";

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::string fileOutName = "fileOut";

  MultiCorrelation multicor;
  multicor.SetIsFillSysErrors(false);
  multicor.SetPalette( {kBlue, kBlue, kBlue, kBlue, kRed, kRed, kRed, kRed, kGreen+2, kGreen+2, kGreen+2, kGreen+2, kMagenta, kMagenta, kMagenta, kMagenta} );
  multicor.SetMarkers( {-1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare} );

  for(int i=0; i<correls.size(); i++) {
    for(auto& comp : components) {
      multicor.AddCorrelation(fileName, {"R1R1/r1r1." + correls.at(i) + "_RECENTERED_" + sub4th + "." + comp}, correls.at(i) + "_" + sub4th);
      multicor.AddCorrelation(fileName, {"QQ/qq." + correls.at(i) + "_RECENTERED_" + sub4th + "." + comp}, correls.at(i) + "_" + sub4th);
    }
  }

  multicor.Scale(-1);

  for(auto& comp : components) {
    multicor.AddCorrelation(fileName, {"R1R1/r1r1.psd1_RECENTERED_psd3_RECENTERED." + comp}, "psd1_psd3");
    multicor.AddCorrelation(fileName, {"QQ/qq.psd1_RECENTERED_psd3_RECENTERED." + comp}, "psd1_psd3");
  }

  multicor.Rebin({{"SimEventHeader_centrality_impactpar", {0, 10, 20, 30, 40, 50, 60, 70}}});

  multicor.SlightShiftXAxis(0.);
  HeapPicture pic("picture", {1000, 1000});
  if(evegen == "dcmqgsm") {
    if(pbeam == "12") pic.AddText({0.18, 0.92, "5M Au+Au"}, 0.025);
    else              pic.AddText({0.18, 0.92, "5.2M Au+Au"}, 0.025);
    pic.AddText({0.18, 0.89, "DCM-QGSM-SMM"}, 0.025);
  }
  if(evegen == "urqmd") {
    pic.AddText({0.18, 0.92, "2M Au+Au"}, 0.025);
    pic.AddText({0.18, 0.89, "UrQMD"}, 0.025);
  }
  pic.AddText({0.18, 0.86, (pbeam + "A GeV/c").c_str()}, 0.025);

  auto leg1 = new TLegend();
  leg1->SetBorderSize(1);

  for(auto& obj : multicor.GetCorrelations()) {
    obj->SetCalculateSystematicsFromVariation(false);
    pic.AddDrawable(obj);
  }

  TLegendEntry* entry = new TLegendEntry();
  entry->SetLineWidth(3);
  entry->SetLineColor(kBlack);
  entry->SetLineStyle(1);
  leg1->AddEntry(entry, "#LTQ#Psi#GT#LTQ#Psi#GTxx", "L");
  entry = new TLegendEntry();
  entry->SetLineStyle(2);
  leg1->AddEntry(entry, "#LTQ#Psi#GT#LTQ#Psi#GTyy", "L");
  entry = new TLegendEntry();
  entry->SetMarkerSize(2);
  entry->SetMarkerColor(kBlack);
  entry->SetMarkerStyle(kFullSquare);
  leg1->AddEntry(entry, "#LTQQ#GTxx", "P");
  entry = new TLegendEntry();
  entry->SetMarkerStyle(kOpenSquare);
  leg1->AddEntry(entry, "#LTQQ#GTyy", "P");
  leg1->AddEntry(multicor.GetCorrelations().at(0)->GetPoints(), multicor.GetCorrelations().at(0)->GetTitle().c_str(), "L");
  leg1->AddEntry(multicor.GetCorrelations().at(4)->GetPoints(), multicor.GetCorrelations().at(4)->GetTitle().c_str(), "L");
  leg1->AddEntry(multicor.GetCorrelations().at(8)->GetPoints(), multicor.GetCorrelations().at(8)->GetTitle().c_str(), "L");
  leg1->AddEntry(multicor.GetCorrelations().at(12)->GetPoints(), multicor.GetCorrelations().at(12)->GetTitle().c_str(), "L");

  pic.SetAxisTitles( {"Centrality, %", ""} );
  pic.CustomizeXRange();
  pic.CustomizeYRange();
  pic.AddLegend(leg1);
  pic.CustomizeLegend(leg1);
  pic.Draw();

  if(is_write_rootfile) {
    TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
    fileOut->cd();
    pic.GetCanvas()->Write();
    fileOut->Close();
  }

  pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
}
