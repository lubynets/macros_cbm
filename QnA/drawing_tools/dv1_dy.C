#include "lambda.h"

void dv1_dy() {
  gROOT->Macro( "/home/user/cbmdir/flow_drawing_tools/example/style.cc" );
  bool is_first_canvas = true;
  
//   std::string fileMcName = "/home/user/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.lightcuts1.310.set4.root";
//   std::string fileRecName = "/home/user/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.apr20.recpid.lightcuts1.310.set4.root";
  
  std::string fileMcName = "/home/user/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.3.3agev.apr20.recpid.lightcuts1.3122.set4.root";
  std::string fileRecName = "/home/user/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.3.3agev.apr20.recpid.lightcuts1.3122.set4.root";
  
//   std::string fileMcName = "/home/user/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.lightcuts1.3122.set4.all.root";
//   std::string fileRecName = "/home/user/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.apr20.recpid.lightcuts1.3122.set4.root";
  
//   std::string fileMcName = "/home/user/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.defaultcuts.3312.set4.root";
//   std::string fileRecName = "/home/user/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.apr20.recpid.lightcuts1.310.set4.root";
  
//   std::string number = "1";
//   std::string number = "2";
//   std::string number = "3";
//   std::string number = "4";
//   std::string number = "5";
//   std::string number = "23";
//   std::string number = "451";
//   std::string number = "2345";
//   std::string number = "12345";
  
//   std::string fileMcName = "/home/user/cbmdir/working/qna/fairreview/3122_12agev/" + number + "/v1andR1.root";
//   std::string fileRecName = "/home/user/cbmdir/working/qna/fairreview/3122_12agev/" + number + "/out.fitter.root";
  
  std::string particle = "#Lambda";
//   std::string particle = "K^{0}_{S}";
//   std::string particle = "#Xi^{-}";
  
  const float C_lo = 0;
  const float C_hi = 40;
  const float pT_lo = 0.;
  const float pT_hi = 1.;
  
  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");
  std::string component_1 = "x1x1";
  std::string component_2 = "y1y1";
  std::string harmonic = "1";
  
  TFile* fileMc = TFile::Open(fileMcName.c_str(), "open");
  auto* dc_pre = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileMc->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN.x1x1");
//   auto dc_value = dc_pre->Rebin({"AnaEventHeader_centrality_tracks", {0,20,40,70}});
  auto dc_value = dc_pre->Rebin({"AnaEventHeader_centrality_tracks", {C_lo, C_hi}});
  auto* dc = &dc_value;
  

  
  
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }
  
//   IntegrateSelectAxis();
//   IntegrateSliceAxis();

  SetSliceAxisBinEdges({pT_lo, pT_hi});

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
//   TFile* fileOut = TFile::Open(("fileOut.root").c_str(), "recreate");
  
for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){
  
    auto v1_sim = DoubleDifferentialCorrelation( fileMcName.c_str(),
                                                {("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN." + component_1).c_str(),
                                                 ("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN." + component_2).c_str() } );
    v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
//     v1_sim.Scale(2);
    v1_sim.SetMarker(-1);
    v1_sim.SetPalette({kBlue});
    v1_sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_sim.Calculate();
    v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);  
   
    
    auto v1_rec_fit = DoubleDifferentialCorrelation( fileRecName.c_str(),
                                                {("parameters/signal." + component_1).c_str(),
                                                 ("parameters/signal." + component_2).c_str() } );
    
//     auto v1_rec_fit = DoubleDifferentialCorrelation( fileMcName.c_str(),
//                                                 {("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_1).c_str(),
//                                                  ("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_2).c_str() } );
    
    
    v1_rec_fit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
//     v1_rec_fit.Scale(2);
    v1_rec_fit.SetMarker(kFullSquare);
    v1_rec_fit.SetPalette({kBlue});
    v1_rec_fit.Rebin({{axes.at(kSelect).reco_fit_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_rec_fit.SetProjectionAxis({axes.at(kProjection).reco_fit_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_rec_fit.SetSliceAxis({axes.at(kSlice).reco_fit_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_rec_fit.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_rec_fit.Calculate();
    v1_rec_fit.ShiftProjectionAxis(axes.at(kProjection).shift_);
    
    
    HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});
    
    pic.AddText({0.2, 0.90, particle.c_str()}, 0.035);
    pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
    pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
//     pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);
    pic.AddText({0.2, 0.81, "3.3A GeV/c"}, 0.025);

    pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + 
                            " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
    
    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);
    
    auto* entry = leg1->AddEntry("", "MC input", "L");
    entry->SetLineColor(kBlack);
    entry->SetLineWidth(2);
        
    entry = leg1->AddEntry("", "REC, invmassfit", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerStyle(kFullSquare);
    
    leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());
    
    for( auto obj : v1_sim.GetProjections() ){
//       obj->SetErrorOption("E3");
      pic.AddDrawable( obj );
    }    

    for( auto obj : v1_rec_fit.GetProjections() ){
      pic.AddDrawable( obj );
      leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
    }

    TF1* f1 = new TF1("f1", "pol1", axes.at(kProjection).bin_edges_.front()-0.01, axes.at(kProjection).bin_edges_.back()+0.01);
//     TF1* f1 = new TF1("f1", "[1]*(x-[0])", axes.at(kProjection).bin_edges_.front()-0.01, axes.at(kProjection).bin_edges_.back()+0.01);
//     f1->FixParameter(0, 0.98);
    v1_rec_fit.GetProjections().at(0)->GetPoints()->Fit(f1, "0");
    pic.AddFunction(f1);
    
    pic.AddText({0.7, 0.2, ("slope = " + to_string_with_precision(f1->GetParameter(1), 6)).c_str()}, 0.025);
    pic.AddText({0.7, 0.17, ("slopeerror = " + to_string_with_precision(f1->GetParError(1), 6)).c_str()}, 0.025);

    pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), ("v_{" + harmonic + "}").c_str()});
    
//     pic.SetXRange({-0.05, 0.95});
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
  
  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print("fileOut.pdf)", "pdf");
  
  fileOut->Close();  
  
}
