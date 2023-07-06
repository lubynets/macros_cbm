void resolution_4sub() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12"; std::string cuts = "lc1";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3"; std::string cuts = "oc1";
  std::string evegen = "urqmd";  std::string pbeam = "12"; std::string cuts = "lc1";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/vR." + evegen + "." + pbeam + "agev." + cuts + ".3122.root";
  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
  std::vector<std::string> fourth{"sts_pipos"};

  bool average_comp = false;
  std::string step = "_RECENTERED";

  std::vector<std::string> components{"x1x1", "y1y1"};

  std::string fileOutName = "res_4sub.psd";
  TFile* fileOut{nullptr};
  if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

  MultiCorrelation multicor_mc;
  multicor_mc.SetIsFillSysErrors(false);

  if(!average_comp) {
    multicor_mc.SetPalette( {kBlue, kBlue, kRed, kRed, kGreen+2, kGreen+2} );
    multicor_mc.SetMarkers( {-1, -2, -1, -2, -1, -2} );
  } else {
    multicor_mc.SetPalette( {kBlue, kRed, kGreen+2} );
    multicor_mc.SetMarkers( {-1, -1, -1} );
  }

  for(auto& corr : correls) {
    if(!average_comp) {
      for(auto& comp : components) {
        multicor_mc.AddCorrelation(fileName, {"R1/MC/res_MC." + corr + "." + comp}, "mc_" + corr + step);
      }
    } else {
      multicor_mc.AddCorrelation(fileName, {"R1/MC/res_MC." + corr + "." + components.at(0),
                                            "R1/MC/res_MC." + corr + "." + components.at(1)}, "mc_" + corr + step);
    }
  }

  bool is_first_canvas = true;
  for(auto& s4 : fourth) {
    MultiCorrelation multicor_sub4;
    multicor_sub4.SetIsFillSysErrors(false);
    if(!average_comp) {
      multicor_sub4.SetPalette( {kBlue, kBlue, kRed, kRed, kGreen+2, kGreen+2} );
      multicor_sub4.SetMarkers( {kFullSquare, kOpenSquare, kFullSquare, kOpenSquare, kFullSquare, kOpenSquare} );
    } else {
      multicor_sub4.SetPalette( {kBlue, kRed, kGreen+2} );
      multicor_sub4.SetMarkers( {kFullSquare, kFullSquare, kFullSquare} );
    }

    for(auto& corr : correls) {
      if(!average_comp) {
        for(auto& comp : components) {
          multicor_sub4.AddCorrelation(fileName, {"R1/sub4/res.sub4." + s4 + "." + corr + "." + comp}, "sub4_" + corr + step);
        }
      } else {
        multicor_sub4.AddCorrelation(fileName, {"R1/sub4/res.sub4." + s4 + "." + corr + "." + components.at(0),
                                                "R1/sub4/res.sub4." + s4 + "." + corr + "." + components.at(1)}, "sub4_" + corr + step);
      }
    }
    multicor_sub4.SlightShiftXAxis(0.3);

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
    pic.AddText({0.18, 0.83, "MC: R^{A}_{x} = 2#LTQ^{A}_{x}Q^{#Psi}_{x}#GT"}, 0.02);

    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);

    for(auto& obj : multicor_mc.GetCorrelations()) {
      obj->SetCalculateSystematicsFromVariation(false);
      pic.AddDrawable(obj);
    }
    for(auto& obj : multicor_sub4.GetCorrelations()) {
      obj->SetCalculateSystematicsFromVariation(false);
      pic.AddDrawable(obj);
    }

    TGraph* grx = new TGraph();
    grx->SetMarkerStyle(kFullSquare);
    grx->SetLineStyle(1);
    grx->SetMarkerColor(kBlack);
    grx->SetLineColor(kBlack);
    TGraph* gry = new TGraph();
    gry->SetMarkerStyle(kOpenSquare);
    gry->SetLineStyle(2);
    gry->SetMarkerColor(kBlack);
    gry->SetLineColor(kBlack);

    leg1->AddEntry(grx, "     MC        ", "L");
    leg1->AddEntry(grx, ("     4-sub_" + s4).c_str(), "P");

    if(!average_comp) {
      leg1->AddEntry(grx, ("     " + components.at(0)).c_str(), "LP");
      leg1->AddEntry(gry, ("     " + components.at(1)).c_str(), "LP");
      leg1->AddEntry(multicor_sub4.GetCorrelations().at(0)->GetPoints(), ("     " + correls.at(0)).c_str(), "P");
      leg1->AddEntry(multicor_sub4.GetCorrelations().at(2)->GetPoints(), ("     " + correls.at(1)).c_str(), "P");
      leg1->AddEntry(multicor_sub4.GetCorrelations().at(4)->GetPoints(), ("     " + correls.at(2)).c_str(), "P");
    } else {
      leg1->AddEntry(grx, "     x&y ave", "LP");
      leg1->AddEntry(multicor_sub4.GetCorrelations().at(0)->GetPoints(), ("     " + correls.at(0)).c_str(), "P");
      leg1->AddEntry(multicor_sub4.GetCorrelations().at(1)->GetPoints(), ("     " + correls.at(1)).c_str(), "P");
      leg1->AddEntry(multicor_sub4.GetCorrelations().at(2)->GetPoints(), ("     " + correls.at(2)).c_str(), "P");
    }

    pic.SetAxisTitles( {"Centrality, %", "R_{1}"} );
    pic.CustomizeXRange();
    pic.CustomizeYRangeWithLimits(0, 0.3);
    pic.AddLegend(leg1);
    pic.CustomizeLegend(leg1);
  //   pic.SetGridX();
  //   pic.SetGridY();
    pic.Draw();

    if(is_write_rootfile) {
      fileOut->cd();
      pic.GetCanvas()->Write();
    }

    if(is_first_canvas) pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
    else pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
    is_first_canvas = false;
  }
  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

  if(is_write_rootfile) fileOut->Close();
}
