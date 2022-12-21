void resolution_presubevent() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

//   std::string evegen = "dcmqgsm";
  std::string evegen = "urqmd";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".root";

  std::vector<std::string> correls{"psd1_RECENTERED", "psd2_RECENTERED", "psd3_RECENTERED"};
//   std::vector<std::string> correls{"spec1_prim_PLAIN", "spec2_prim_PLAIN", "spec3_prim_PLAIN"};

  int n_substr;
  if(correls.at(0)[0] == 'p') n_substr = 4;
  if(correls.at(0)[0] == 's') n_substr = 5;

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::string fileOutName = "fileOut";
  bool is_first_canvas = true;

  MultiCorrelation multicor_AA, multicor_AB;
  multicor_AA.SetPalette( {kBlue, kBlue, kBlue, kBlue, kRed, kRed, kRed, kRed, kGreen+2, kGreen+2, kGreen+2, kGreen+2} );
  multicor_AB.SetPalette( {kBlue, kBlue, kBlue, kBlue, kRed, kRed, kRed, kRed, kGreen+2, kGreen+2, kGreen+2, kGreen+2} );
  multicor_AA.SetMarkers( {-1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare} );
  multicor_AB.SetMarkers( {-1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare, -1, kFullSquare, -2, kOpenSquare} );

  for(int i=0; i<correls.size(); i++) {
    for(int j=i; j<correls.size(); j++) {
      for(auto& comp : components) {
        if(i == j) {
          multicor_AA.AddCorrelation(fileName, {"R1R1/r1r1." + correls.at(i) + "_" + correls.at(j) + "." + comp}, correls.at(i).substr(0, n_substr) + "_" + correls.at(j).substr(0, n_substr));
          multicor_AA.AddCorrelation(fileName, {"QQ/qq." + correls.at(i) + "_" + correls.at(j) + "." + comp}, correls.at(i).substr(0, n_substr) + "_" + correls.at(j).substr(0, n_substr));
        }
        else {
          multicor_AB.AddCorrelation(fileName, {"R1R1/r1r1." + correls.at(i) + "_" + correls.at(j) + "." + comp}, correls.at(i).substr(0, n_substr) + "_" + correls.at(j).substr(0, n_substr));
          multicor_AB.AddCorrelation(fileName, {"QQ/qq." + correls.at(i) + "_" + correls.at(j) + "." + comp}, correls.at(i).substr(0, n_substr) + "_" + correls.at(j).substr(0, n_substr));
        }
      }
    }
  }

  for(int i=0; i<2; i++) {
    MultiCorrelation multicor;
    if(i == 0) multicor = multicor_AB;
    if(i == 1) multicor = multicor_AA;
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

    pic.SetAxisTitles( {"Centrality, %", ""} );
    pic.CustomizeXRange();
    pic.CustomizeYRange();
    pic.AddLegend(leg1);
    pic.CustomizeLegend(leg1);
    pic.Draw();

    if(is_first_canvas) {
      pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
      is_first_canvas = false;
    }
    else {
      pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
    }
  }

  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + ".pdf)").c_str(), "pdf");
}
