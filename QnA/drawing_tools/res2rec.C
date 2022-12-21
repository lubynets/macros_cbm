void res2rec()
{
  gROOT->Macro( "/home/user/cbmdir/flow_drawing_tools/example/style.cc" );
  
  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";
  
  std::string pbeam = "12";
//   std::string pbeam = "3.3";
  
//   std::string fileIn = "/home/user/cbmdir/working/qna/resolutions/psd." + evegen + "." + pbeam + "agev.root";
  std::string fileIn = "/home/user/cbmdir/working/qna/resolutions/res." + evegen + "." + pbeam + "agev.root";
    
  std::vector<std::string> steps{"RECENTERED", "RESCALED"};
  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
  std::vector<std::string> components{"x1x1", "y1y1"};
  
  bool is_first_canvas = true;
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  for(auto& corr : correls)
    for(auto& comp : components)
    {
//       auto corr_rece = new Correlation(fileIn, {"Qpsi/" + corr + "_RECENTERED.Q_psi_PLAIN." + comp}, "mc_" + corr + "_RECENTERED_" + comp);
//       auto corr_resc = new Correlation(fileIn, {"Qpsi/" + corr + "_RESCALED.Q_psi_PLAIN." + comp}, "mc_" + corr + "_RESCALED_" + comp);
      
      auto corr_rece = new Correlation(fileIn, {corr + "_RECENTERED." + comp}, "mc_" + corr + "_RECENTERED_" + comp);
      auto corr_resc = new Correlation(fileIn, {corr + "_RESCALED." + comp}, "mc_" + corr + "_RESCALED_" + comp);
      
      Ratio<Correlation> pic(corr + "_" + comp, {1000, 1000});
      
      pic.AddText({0.2, 0.90, corr.c_str()}, 0.025);
      pic.AddText({0.2, 0.87, comp.c_str()}, 0.025);
      pic.AddText({0.6, 0.90, "DCM-QGSM-SMM"}, 0.025);
      pic.AddText({0.6, 0.87, "12A GeV/c"}, 0.025);
//       pic.AddText({0.6, 0.90, "UrQMD"}, 0.025);
//       pic.AddText({0.6, 0.87, "3.3A GeV/c"}, 0.025);

      pic.SetReference(corr_rece);
      pic.AddObject(corr_resc);
      
      pic.SetAxisTitles( {"Centrality, %", "Q_{PSD}Q_{#Psi}", "ratio"} );
      
      pic.CustomizeXRange();
      pic.CustomizeYRange();  
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