void rectwires()
{
  gROOT->Macro( "/home/user/cbmdir/flow_drawing_tools/example/style.cc" );
  
  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";
  
  std::string pbeam = "12";
//   std::string pbeam = "3.3";
  
//   std::string fileName = "/home/user/cbmdir/working/qna/resolutions/psd." + evegen + "." + pbeam + "agev.root";
  std::string fileName = "/home/user/cbmdir/working/qna/resolutions/res." + evegen + "." + pbeam + "agev.root";
    
//   std::vector<std::string> steps{"PLAIN", "RECENTERED", "TWIST", "RESCALED"};
  std::vector<std::string> steps{"RECENTERED", "RESCALED"};
  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
  std::vector<std::string> components{"x1x1", "y1y1"};
  
  bool is_first_canvas = true;
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  for(auto& corr : correls)
    for(auto& comp : components)
    {
      MultiCorrelation multicor_mc;
//       multicor_mc.SetPalette( {kRed, kGreen+2, kBlue, kBlack } );
//       multicor_mc.SetMarkers( {-1, -2, kOpenSquare, kOpenCircle, -1, -1, -1, -1} );
      multicor_mc.SetPalette( {kRed, kBlack } );
      multicor_mc.SetMarkers( {-1, kOpenCircle} );
      
//       for(auto& step : steps)
//         multicor_mc.AddCorrelation(fileName, {"Qpsi/" + corr + "_" + step + ".Q_psi_PLAIN." + comp}, "mc_" + corr + "_" + step + "_" + comp);
      
      for(auto& step : steps)
        multicor_mc.AddCorrelation(fileName, {corr + "_" + step + "." + comp}, "mc_" + corr + "_" + step + "_" + comp);
      
      multicor_mc.Scale(2);
      
      HeapPicture pic(corr + "_" + comp, {1000, 1000});
      
      pic.AddText({0.2, 0.90, corr.c_str()}, 0.025);
      pic.AddText({0.2, 0.87, comp.c_str()}, 0.025);
      pic.AddText({0.6, 0.90, "DCM-QGSM-SMM"}, 0.025);
      pic.AddText({0.6, 0.87, "12A GeV/c"}, 0.025);
//       pic.AddText({0.6, 0.90, "UrQMD"}, 0.025);
//       pic.AddText({0.6, 0.87, "3.3A GeV/c"}, 0.025);
      
      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);
      
      for(auto& obj : multicor_mc.GetCorrelations()) {
        obj->SetCalculateSystematicsFromVariation(false);
        pic.AddDrawable(obj);
      }
      
//       leg1->AddEntry(multicor_mc.GetCorrelations().at(0)->GetPoints(), steps.at(0).c_str(), "L");
//       leg1->AddEntry(multicor_mc.GetCorrelations().at(1)->GetPoints(), steps.at(1).c_str(), "L");
//       leg1->AddEntry(multicor_mc.GetCorrelations().at(2)->GetPoints(), steps.at(2).c_str(), "P");      
//       leg1->AddEntry(multicor_mc.GetCorrelations().at(3)->GetPoints(), steps.at(3).c_str(), "P");
     
      leg1->AddEntry(multicor_mc.GetCorrelations().at(0)->GetPoints(), steps.at(0).c_str(), "L");
      leg1->AddEntry(multicor_mc.GetCorrelations().at(1)->GetPoints(), steps.at(1).c_str(), "P");
      
      pic.SetAxisTitles( {"Centrality, %", "Q_{PSD}Q_{#Psi}"} );
      
      pic.CustomizeXRange();
      pic.CustomizeYRange();  
      pic.AddLegend(leg1);
      pic.CustomizeLegend(leg1);
      pic.SetGridX();
      pic.SetGridY();
      pic.Draw();  
      
      fileOut->cd();
      pic.GetCanvas()->Write();
      
      if(is_first_canvas)
        pic.GetCanvas()->Print("fileOut.pdf(", "pdf");
      else
        pic.GetCanvas()->Print("fileOut.pdf", "pdf");
      is_first_canvas = false;
    }
        
  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print("fileOut.pdf)", "pdf");
  
  fileOut->Close();
}