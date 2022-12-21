void psd_cross() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );
  bool is_first_canvas = true;

//   TFile* fileOut = TFile::Open("fileOut.root", "recreate");

//   std::string evegen = "dcmqgsm";
  std::string evegen = "urqmd";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".root";

  std::vector<std::string> correls{"psd1_RECENTERED", "psd2_RECENTERED", "psd3_RECENTERED"};
//   std::vector<std::string> correls{"spec1_prim_PLAIN", "spec2_prim_PLAIN", "spec3_prim_PLAIN"}; // TODO not calculated yet

  std::string fileOutName;
  if(correls.at(0)[0] == 'p') fileOutName = "cross.psd";
  if(correls.at(0)[0] == 's') fileOutName = "cross.spec_prim";

  for(int i=0; i<correls.size(); i++) {
    MultiCorrelation multicor;
    multicor.SetPalette( {kBlue, kRed, kGreen+2} );
    multicor.SetMarkers( {kFullSquare, kFullSquare, kFullSquare} );

    for(int j=0; j<correls.size(); j++) {
      multicor.AddCorrelation(fileName, {"R1/res_cross." + correls.at(i) + "." + correls.at(j) + ".x1y1"}, "mc_" + correls.at(i) + "." + correls.at(j) + "_x1y1");
    }

    multicor.SlightShiftXAxis(1);

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
    pic.AddText({0.18, 0.83, ("x1: " + correls.at(i)).c_str()}, 0.025);

    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);

    for(auto& obj : multicor.GetCorrelations()) {
      obj->SetCalculateSystematicsFromVariation(false);
      pic.AddDrawable(obj);
    }

    leg1->SetHeader("y1:");
    leg1->AddEntry(multicor.GetCorrelations().at(0)->GetPoints(), correls.at(0).c_str(), "P");
    leg1->AddEntry(multicor.GetCorrelations().at(1)->GetPoints(), correls.at(1).c_str(), "P");
    leg1->AddEntry(multicor.GetCorrelations().at(2)->GetPoints(), correls.at(2).c_str(), "P");

    pic.SetAxisTitles( {"Centrality, %", "Q_{x}Q_{y}"} );
    pic.CustomizeXRange();
    pic.CustomizeYRange();
    pic.AddLegend(leg1);
    pic.CustomizeLegend(leg1);
    pic.Draw();

//     fileOut->cd();
//     pic.GetCanvas()->Write();

    if(is_first_canvas)
      pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
    else
      pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
    is_first_canvas = false;
  }

  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + ".pdf)").c_str(), "pdf");

//   fileOut->Close();
}
