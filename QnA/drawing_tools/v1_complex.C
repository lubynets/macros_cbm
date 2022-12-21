#include "v1_complex.h"

void v1_complex()
{
  gROOT->Macro( "/home/user/cbmdir/flow_drawing_tools/example/style.cc" );
  
  bool is_first_canvas = true;  
  
  std::string fileIn = "/home/user/cbmdir/working/qna/correlations/derivatives/v1.dcmqgsm.12agev.defcuts.3122.set2.sgnl_1.root";
  
  std::vector<std::string> correls{"psd1", "psd2", "psd3"};
  std::vector<std::string> components{"x1x1", "y1y1"};
  
  std::vector<std::string> subevents{"psd1", "psd2", "psd3"};
  std::vector<std::string> methods{"sub3", "sub4_sts_p", "sub4_sts_pipos"};
  
  const float C_lo = 10;
  const float C_hi = 20;
  const float pT_lo = 0.4;
  const float pT_hi = 1.0;
//   const float y_lo = 0.92179005;
//   const float y_hi = 2.7217901;
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  for(auto& me : methods) {
    for(auto& se : subevents) {
      
      MultiCorrelation usimQ;
      
      usimQ.AddCorrelation(fileIn, {"v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN." + components.at(0),
                                    "v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN." + components.at(1)}, "usimPsi");
      
      usimQ.SetPalette({kRed});
      usimQ.SetMarkers({-1});
      
      usimQ.Rebin({{"AnaEventHeader_centrality_tracks", 1, C_lo, C_hi}});
      usimQ.Rebin({{"SimParticles_pT", 1, pT_lo, pT_hi}});
      usimQ.Project({"SimParticles_rapidity"});
      
      MultiCorrelation urecQ;
      
      urecQ.AddCorrelation(fileIn, {"v1/uPsi/v1.u_rec_RESCALED.Q_psi_PLAIN." + components.at(0),
                                    "v1/uPsi/v1.u_rec_RESCALED.Q_psi_PLAIN." + components.at(1)}, "uPsi");      
      
      urecQ.AddCorrelation(fileIn, {"v1/uQ_R1_MC/v1.u_rec_RESCALED." + se + "_RECENTERED.res_MC." + components.at(0),
                                    "v1/uQ_R1_MC/v1.u_rec_RESCALED." + se + "_RECENTERED.res_MC." + components.at(1)}, "uQ_MC");
           
      
      urecQ.AddCorrelation(fileIn, {"v1/uQ_R1_" + me + "/v1.u_rec_RESCALED." + se + "_RECENTERED.res_" + me + "." + components.at(0),
                                    "v1/uQ_R1_" + me + "/v1.u_rec_RESCALED." + se + "_RECENTERED.res_" + me + "." + components.at(1)}, "uQ_res");
      
      urecQ.SetPalette({kGreen+2, kBlue, kRed});
      urecQ.SetMarkers({kFullTriangleUp, kFullCircle, kOpenSquare});
      
      urecQ.Rebin({{"AnaEventHeader_centrality_tracks", 1, C_lo, C_hi}});
      urecQ.Rebin({{"ReconstructedParticles_pT", 1, pT_lo, pT_hi}});
      urecQ.Project({"ReconstructedParticles_rapidity"});
  
      HeapPicture pic(me + "_" + se, {1000, 1000});
      
      pic.AddText({0.2, 0.90, "#Lambda"}, 0.035);
      pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
      pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
      pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);
      pic.AddText({0.2, 0.78, ("p_{T}: " + to_string_with_precision(pT_lo, 1) + " - " + to_string_with_precision(pT_hi, 1) + " GeV/c").c_str()}, 0.025);
      pic.AddText({0.2, 0.75, ("centrality: " + to_string_with_precision(C_lo, 0) + " - " + to_string_with_precision(C_hi, 0) + " %").c_str()}, 0.025);
      pic.AddText({0.2, 0.72, me.c_str()}, 0.025);
      pic.AddText({0.2, 0.69, se.c_str()}, 0.025);
      
      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);
      
      for(auto& obj : usimQ.GetCorrelations()) {
        obj->SetCalculateSystematicsFromVariation(false);
        Graph* gr = new Graph(obj);
        gr->ShiftXaxis(-1.6217901);
        pic.AddDrawable(gr);
      }
            
      for(auto& obj : urecQ.GetCorrelations()) {
        obj->SetCalculateSystematicsFromVariation(false);
        Graph* gr = new Graph(obj);
        gr->ShiftXaxis(-1.6217901);
        pic.AddDrawable(gr);
      }
      
      leg1->AddEntry( usimQ.GetCorrelations().at(0)->GetPoints(), "MC", "L" );
      leg1->AddEntry( urecQ.GetCorrelations().at(0)->GetPoints(), "#Psi^{RP}", "P" );
      leg1->AddEntry( urecQ.GetCorrelations().at(1)->GetPoints(), "#Psi^{PSD}(R_{1}^{true})", "P" );
      leg1->AddEntry( urecQ.GetCorrelations().at(2)->GetPoints(), "#Psi^{PSD}(R_{1}^{reco})", "P" );
            
      pic.SetAxisTitles( {"y", "v_{1}"} );
      pic.CustomizeXRange();
      pic.CustomizeYRange();
      pic.AddLegend(leg1);
      pic.CustomizeLegend(leg1);
      pic.Draw();
      
      fileOut->cd();
      
      pic.GetCanvas()->Write();
      
      if(is_first_canvas)
        pic.GetCanvas()->Print("fileOut.pdf(", "pdf");
      else
        pic.GetCanvas()->Print("fileOut.pdf", "pdf");
      is_first_canvas = false;
    }
  }

  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print("fileOut.pdf)", "pdf");
  
  fileOut->Close();    
}