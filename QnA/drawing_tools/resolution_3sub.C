void resolution_3sub() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;
    
  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".nopass.root";

//   std::vector<std::string> correls{"psd1", "psd2", "psd3"};
//   std::vector<std::string> correls{"etacut_1_charged", "etacut_2_charged", "etacut_3_charged"};
  std::vector<std::string> correls{"etacut_1_all", "etacut_2_all", "etacut_3_all"};

  std::string step;
  bool average_comp;

  std::vector<std::string> components{"x1x1", "y1y1"};

  std::string fileOutName;
  if(correls.at(0)[0] == 'p') { step = "_RECENTERED"; fileOutName = "res_3sub.psd"; average_comp = false; }
  if(correls.at(0)[0] == 'e') { step = "_PLAIN"; fileOutName = "res_3sub.etacut"; average_comp = true; }

  MultiCorrelation multicor_mc, multicor_sub3;
  multicor_mc.SetIsFillSysErrors(false);
  multicor_sub3.SetIsFillSysErrors(false);
  if(!average_comp) {
    multicor_mc.SetPalette( {kBlue, kBlue, kRed, kRed, kGreen+2, kGreen+2} );
    multicor_mc.SetMarkers( {-1, -2, -1, -2, -1, -2} );
    multicor_sub3.SetPalette( {kBlue, kBlue, kRed, kRed, kGreen+2, kGreen+2} );
    multicor_sub3.SetMarkers( {kFullSquare, kOpenSquare, kFullSquare, kOpenSquare, kFullSquare, kOpenSquare} );
  } else {
    multicor_mc.SetPalette( {kBlue, kRed, kGreen+2} );
    multicor_mc.SetMarkers( {-1, -1, -1} );
    multicor_sub3.SetPalette( {kBlue, kRed, kGreen+2} );
    multicor_sub3.SetMarkers( {kFullSquare, kFullSquare, kFullSquare} );
  }

  for(auto& corr : correls) {
    if(!average_comp) {
      for(auto& comp : components) {
        multicor_mc.AddCorrelation(fileName, {"R1/res.mc." + corr + step + "." + comp}, "mc_" + corr + step);
        multicor_sub3.AddCorrelation(fileName, {"R1/res.sub3." + corr + step + "." + comp}, "mc_" + corr + step);
      }
    } else {
      multicor_mc.AddCorrelation(fileName, {"R1/res.mc." + corr + step + "." + components.at(0),
                                            "R1/res.mc." + corr + step + "." + components.at(1)}, "mc_" + corr + step);
      multicor_sub3.AddCorrelation(fileName, {"R1/res.sub3." + corr + step + "." + components.at(0),
                                              "R1/res.sub3." + corr + step + "." + components.at(1)}, "mc_" + corr + step);
    }
  }

//   multicor_mc.SlightShiftXAxis(0.);
            
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
  pic.AddText({0.18, 0.83, "MC: R^{A}_{x} = 2#LTQ^{A}_{x}Q^{#Psi}_{x}#GT"}, 0.02);
  
  auto leg1 = new TLegend();
  leg1->SetBorderSize(1);
  
  for(auto& obj : multicor_mc.GetCorrelations()) {
    obj->SetCalculateSystematicsFromVariation(false);
    pic.AddDrawable(obj);
  }
  for(auto& obj : multicor_sub3.GetCorrelations()) {
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
  leg1->AddEntry(grx, "     3-sub", "P");
  
  if(!average_comp) {
    leg1->AddEntry(grx, ("     " + components.at(0)).c_str(), "LP");
    leg1->AddEntry(gry, ("     " + components.at(1)).c_str(), "LP");
    leg1->AddEntry(multicor_sub3.GetCorrelations().at(0)->GetPoints(), ("     " + correls.at(0)).c_str(), "P");
    leg1->AddEntry(multicor_sub3.GetCorrelations().at(2)->GetPoints(), ("     " + correls.at(1)).c_str(), "P");
    leg1->AddEntry(multicor_sub3.GetCorrelations().at(4)->GetPoints(), ("     " + correls.at(2)).c_str(), "P");
  } else {
    leg1->AddEntry(grx, "     x&y ave", "LP");
    leg1->AddEntry(multicor_sub3.GetCorrelations().at(0)->GetPoints(), ("     " + correls.at(0)).c_str(), "P");
    leg1->AddEntry(multicor_sub3.GetCorrelations().at(1)->GetPoints(), ("     " + correls.at(1)).c_str(), "P");
    leg1->AddEntry(multicor_sub3.GetCorrelations().at(2)->GetPoints(), ("     " + correls.at(2)).c_str(), "P");
  }
  
  pic.SetAxisTitles( {"Centrality, %", "R_{1}"} );
  pic.CustomizeXRange();
  pic.CustomizeYRange();  
  pic.AddLegend(leg1);
  pic.CustomizeLegend(leg1);
//   pic.SetGridX();
//   pic.SetGridY();
  pic.Draw();
//   pic.Save("fileOut", "png");
  
  if(is_write_rootfile) {
    TFile* fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");
    fileOut->cd();
    pic.GetCanvas()->Write();
    fileOut->Close();
  }
    
  pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
}
