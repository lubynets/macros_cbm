void resolution_presubevent() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".root";

  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
//   std::vector<std::string> correls{"etacut_1_charged", "etacut_2_charged", "etacut_3_charged"};
//   std::vector<std::string> correls{"etacut_1_all", "etacut_2_all", "etacut_3_all"};

  std::string step;
  bool average_comp;

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::string fileOutName = "fileOut";

  if(correls.at(0)[0] == 'p') { step = "_RECENTERED"; average_comp = false; }
  if(correls.at(0)[0] == 'e') { step = "_PLAIN"; average_comp = true; }

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
  multicor.SlightShiftXAxis(0.);
  HeapPicture pic("picture", {1000, 1000});
  if(evegen == "dcmqgsm") {
    pic.AddText({0.18, 0.92, "5M Au+Au"}, 0.025);
    pic.AddText({0.18, 0.89, "DCM-QGSM-SMM"}, 0.025);
  }
  if(evegen == "urqmd") {
    pic.AddText({0.18, 0.92, "2M Au+Au"}, 0.025);
    pic.AddText({0.18, 0.89, "UrQMD"}, 0.025);
  }
  pic.AddText({0.18, 0.86, "12A GeV/c"}, 0.025);

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
  } else {
    entry->SetLineWidth(3);
    entry->SetLineColor(kBlack);
    entry->SetLineStyle(1);
    leg1->AddEntry(entry, "#LTQ#Psi#GT#LTQ#Psi#GT x&y ave", "L");
    entry = new TLegendEntry();
    entry->SetMarkerSize(2);
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(kFullSquare);
    leg1->AddEntry(entry, "#LTQQ#GT x&y ave", "P");
    leg1->AddEntry(multicor.GetCorrelations().at(0)->GetPoints(), multicor.GetCorrelations().at(0)->GetTitle().c_str(), "L");
    leg1->AddEntry(multicor.GetCorrelations().at(2)->GetPoints(), multicor.GetCorrelations().at(2)->GetTitle().c_str(), "L");
    leg1->AddEntry(multicor.GetCorrelations().at(4)->GetPoints(), multicor.GetCorrelations().at(4)->GetTitle().c_str(), "L");
  }

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
