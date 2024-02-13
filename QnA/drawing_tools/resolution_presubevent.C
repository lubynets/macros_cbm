void resolution_presubevent() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3";
//   std::string evegen = "urqmd";  std::string pbeam = "12";

//   bool is_write_rootfile = false;
  bool is_write_rootfile = true;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/" + pbeam + "agev/v1andR1.stf." + evegen + "." + pbeam + "agev.root";

  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
//   std::vector<std::string> correls{"etacut_1_charged", "etacut_2_charged", "etacut_3_charged"};
//   std::vector<std::string> correls{"etacut_1_all", "etacut_2_all", "etacut_3_all"};

  std::string step;
  bool average_comp;

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::string fileOutName;

  if(correls.at(0)[0] == 'p') { step = "_RECENTERED"; fileOutName = "factoriz_3sub.psd"; average_comp = false; }
  if(correls.at(0)[0] == 'e') { step = "_PLAIN"; fileOutName = "factoriz_3sub.eta"; average_comp = true; }

  MultiCorrelation multicor;
  multicor.SetIsFillSysErrors(false);
  if(!average_comp) {
    multicor.SetPalette( {kBlue, kBlue, kBlue, kBlue, kRed, kRed, kRed, kRed, kGreen+2, kGreen+2, kGreen+2, kGreen+2} );
    multicor.SetMarkers( {-1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare} );
  } else {
    multicor.SetPalette( {kBlue, kBlue, kRed, kRed, kGreen+2, kGreen+2} );
    multicor.SetMarkers( {-1, kFullSquare, -1, kFullSquare, -1, kFullSquare} );
  }

  for(int i=0; i<correls.size(); i++) {
    for(int j=i; j<correls.size(); j++) {
      if(i == j) continue;
      if(!average_comp) {
        for(auto& comp : components) {
          multicor.AddCorrelation(fileName, {"R1R1/r1r1." + correls.at(i) + step + "_" + correls.at(j) + step + "." + comp}, correls.at(i) + "_" + correls.at(j));
          multicor.AddCorrelation(fileName, {"QQ/qq." + correls.at(i) + step + "_" + correls.at(j) + step + "." + comp}, correls.at(i) + "_" + correls.at(j));
        }
      } else {
        multicor.AddCorrelation(fileName, {"R1R1/r1r1." + correls.at(i) + step + "_" + correls.at(j) + step + "." + components.at(0),
                                           "R1R1/r1r1." + correls.at(i) + step + "_" + correls.at(j) + step + "." + components.at(1)}, correls.at(i) + "_" + correls.at(j));
        multicor.AddCorrelation(fileName, {"QQ/qq." + correls.at(i) + step + "_" + correls.at(j) + step + "." + components.at(0),
                                           "QQ/qq." + correls.at(i) + step + "_" + correls.at(j) + step + "." + components.at(1)}, correls.at(i) + "_" + correls.at(j));
      }
    }
  }

  multicor.Rebin({{"SimEventHeader_centrality_impactpar", {0, 10, 20, 30, 40, 50, 60, 70}}});

  multicor.SlightShiftXAxis(0.);
  HeapPicture pic("picture", {1000, 1000});

  const float text_size = 20;
  const int text_font = 63;
//   if(evegen == "dcmqgsm") {
//     if(pbeam == "12") pic.AddText("5M Au+Au", {0.04, 0.96}, text_size, text_font);
//     else              pic.AddText("5.2M Au+Au", {0.04, 0.96}, text_size, text_font);
//     pic.AddText("DCM-QGSM-SMM", {0.04, 0.92}, text_size, text_font);
//   }
//   if(evegen == "urqmd") {
//     pic.AddText("2M Au+Au", {0.04, 0.96}, text_size, text_font);
//     pic.AddText("UrQMD", {0.04, 0.92}, text_size, text_font);
//   }
//   pic.AddText(pbeam + "A GeV/c", {0.04, 0.88}, text_size, text_font);

  auto leg1 = new TLegend();
  leg1->SetBorderSize(1);

  for(auto& obj : multicor.GetCorrelations()) {
    obj->SetCalculateSystematicsFromVariation(false);
    pic.AddDrawable(obj);
  }

  TLegendEntry* entry = new TLegendEntry();
  if(!average_comp) {
    entry->SetLineWidth(3);
    entry->SetLineColor(kBlack);
    entry->SetLineStyle(1);
    leg1->AddEntry(entry, "factorized product, xx", "L");
    entry = new TLegendEntry();
    entry->SetLineStyle(2);
    leg1->AddEntry(entry, "factorized product, yy", "L");
    entry = new TLegendEntry();
    entry->SetMarkerSize(2);
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(kFullSquare);
    leg1->AddEntry(entry, "plain product, xx", "P");
    entry = new TLegendEntry();
    entry->SetMarkerStyle(kOpenSquare);
    leg1->AddEntry(entry, "plain product, yy", "P");
    leg1->AddEntry(multicor.GetCorrelations().at(0)->GetPoints(), ("A_B: " + multicor.GetCorrelations().at(0)->GetTitle()).c_str(), "L");
    leg1->AddEntry(multicor.GetCorrelations().at(4)->GetPoints(), ("A_B: " + multicor.GetCorrelations().at(4)->GetTitle()).c_str(), "L");
    leg1->AddEntry(multicor.GetCorrelations().at(8)->GetPoints(), ("A_B: " + multicor.GetCorrelations().at(8)->GetTitle()).c_str(), "L");
  } else {
    entry->SetLineWidth(3);
    entry->SetLineColor(kBlack);
    entry->SetLineStyle(1);
    leg1->AddEntry(entry, "factorized product, x&y ave", "L");
    entry = new TLegendEntry();
    entry->SetMarkerSize(2);
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(kFullSquare);
    leg1->AddEntry(entry, "plain product, x&y ave", "P");
    leg1->AddEntry(multicor.GetCorrelations().at(0)->GetPoints(), ("A_B: " + multicor.GetCorrelations().at(0)->GetTitle()).c_str(), "L");
    leg1->AddEntry(multicor.GetCorrelations().at(2)->GetPoints(), ("A_B: " + multicor.GetCorrelations().at(2)->GetTitle()).c_str(), "L");
    leg1->AddEntry(multicor.GetCorrelations().at(4)->GetPoints(), ("A_B: " + multicor.GetCorrelations().at(4)->GetTitle()).c_str(), "L");
  }

  pic.SetAxisTitles( {"Centrality, %", "#LTQ_{1}^{A}Q_{1}^{B}#GT"} );
  pic.CustomizeXRange();
  pic.CustomizeYRange();
  pic.AddLegend(leg1);
  pic.SetIsCustomizeLegend();
  pic.Draw();

  if(is_write_rootfile) {
    TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
    fileOut->cd();
    pic.GetCanvas()->Write();
    pic.Write("heap_picture");
    fileOut->Close();
  }

  pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
}
