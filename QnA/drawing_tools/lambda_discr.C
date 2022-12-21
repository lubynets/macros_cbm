#include "lambda.h"

void lambda_discr() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );
  bool is_first_canvas = true;
  
  std::string fileMcName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.lightcuts1.310.set4.root";
  std::string fileRecName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.apr20.recpid.lightcuts1.310.set4_1.root";
  
//   std::string fileMcName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.3.3agev.apr20.recpid.lightcuts1.3122.set4.root";
//   std::string fileRecName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.3.3agev.apr20.recpid.lightcuts1.3122.set4.root";
  
//   std::string fileMcName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.lightcuts1.3122.set4.sgnl.root";
//   std::string fileRecName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.apr20.recpid.lightcuts1.3122.set4_1.root";
  
//   std::string fileMcName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/v1andR1.dcmqgsm.apr20.recpid.defaultcuts.3312.set4.root";
  
  std::string particle = "#Lambda";
//   std::string particle = "K^{0}_{S}";
//   std::string particle = "#Xi^{-}";
  
  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");
  std::string component_1 = "x1x1";
  std::string component_2 = "y1y1";
  std::string harmonic = "1";
  
//   SetAxis("rapidity", "select");
//   SetAxis("pT", "projection");
//   SetAxis("centrality", "slice");
//   std::string component_1 = "x2x2";
//   std::string component_2 = "y2y2";
//   std::string harmonic = "2";

  TFile* fileMc = TFile::Open(fileMcName.c_str(), "open");
  auto* dc_pre = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileMc->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN.x1x1");
  auto dc_value = dc_pre->Rebin({"AnaEventHeader_centrality_tracks", {0,10,20,40,70}});
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

//   SetSelectAxisBinEdges({10, 40});
//   SetSliceAxisBinEdges({0.2, 1.0});
//   SetSliceAxisBinEdges({0.2, 1.4}); // 3122
//   SetSliceAxisBinEdges({0, 0.9}); // 310
//   IntegrateSliceAxis();
  SetSliceAxisBinEdges({0, 0.3, 0.6, 0.9});
  
//   SetProjectionAxisBinEdges({
//     -0.5-axes.at(kProjection).shift_,
//     -0.1-axes.at(kProjection).shift_,
//     0.1-axes.at(kProjection).shift_,
//     0.5-axes.at(kProjection).shift_,
//     0.7-axes.at(kProjection).shift_,
//   });

//   SetProjectionAxisBinEdges({-0.5-axes.at(kProjection).shift_,
//                              -0.3-axes.at(kProjection).shift_,
//                              -0.1-axes.at(kProjection).shift_,
//                               0.7-axes.at(kProjection).shift_});
  
  
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){
  
    auto v1_sim = DoubleDifferentialCorrelation( fileMcName.c_str(),
                                                {("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN." + component_1).c_str(),
                                                 ("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN." + component_2).c_str() } );
    v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
//     v1_sim.Scale(2);
    v1_sim.SetMarker(-1);
    v1_sim.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed,
                       kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed});
    v1_sim.SetBiasPalette(false);
    v1_sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_sim.Calculate();
    v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);  
  

    auto v1_rec_nofit = DoubleDifferentialCorrelation( fileMcName.c_str(),
                                                {("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_sgnl_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_1).c_str(),
                                                 ("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_sgnl_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_2).c_str() } ); // 310, 3312
//     auto v1_rec_nofit = DoubleDifferentialCorrelation( fileMcName.c_str(),
//                                                 {("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_1).c_str(),
//                                                  ("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_2).c_str() } );         // 3122
    
    v1_rec_nofit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
//     v1_rec_nofit.Scale(2);
    v1_rec_nofit.SetMarker(kOpenSquare);
    v1_rec_nofit.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed,
                             kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed});
    v1_rec_nofit.SetBiasPalette(false);
    v1_rec_nofit.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_rec_nofit.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_rec_nofit.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_rec_nofit.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_rec_nofit.Calculate();
    v1_rec_nofit.ShiftProjectionAxis(axes.at(kProjection).shift_);
    
    
    auto v1_rec_fit = DoubleDifferentialCorrelation( fileRecName.c_str(),
                                                {("parameters/signal." + component_1).c_str(),
                                                 ("parameters/signal." + component_2).c_str() } ); // 3122, 310
    
//     auto v1_rec_fit = DoubleDifferentialCorrelation( fileMcName.c_str(),
//                                                 {("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_1).c_str(),
//                                                  ("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd2_RECENTERED.res_sub4_sts_pipos." + component_2).c_str() } ); // 3312
    
    
    v1_rec_fit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
//     v1_rec_fit.Scale(2);
    v1_rec_fit.SetMarker(kFullSquare);
    v1_rec_fit.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed,
                           kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed});
    v1_rec_fit.SetBiasPalette(false);
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
    pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);
//     pic.AddText({0.2, 0.81, "3.3A GeV/c"}, 0.025);

    pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + 
                            " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
    
    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);
    
    auto* entry = leg1->AddEntry("", "MC input", "L");
    entry->SetLineColor(kBlack);
    entry->SetLineWidth(2);
    
    entry = leg1->AddEntry("", "REC, MC-match", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerStyle(kOpenSquare);
    
    entry = leg1->AddEntry("", "REC", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerStyle(kFullSquare);
    
    leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());
    
    for( auto obj : v1_sim.GetProjections() ){
//       obj->SetErrorOption("E3");
      pic.AddDrawable( obj );
    }    
    for( auto obj : v1_rec_nofit.GetProjections() ){
      pic.AddDrawable( obj );
//       leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
    }
    for( auto obj : v1_rec_fit.GetProjections() ){
      pic.AddDrawable( obj );
      leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
    }    


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
