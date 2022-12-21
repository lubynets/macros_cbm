//
// Created by mikhail on 5/6/21.
//

#include "lambda.h"

// void lambda(const std::string& filename){
void lambda(){
  gROOT->Macro( "/home/user/cbmdir/flow_drawing_tools/example/style.cc" );
//   auto leg1 = new TLegend( 0.65, 0.75, 0.9, 0.945 );
  bool is_first_canvas = true;
  
//   SetFileName("/home/user/cbmdir/working/qna/correlations/derivatives/v1.dcmqgsm.12agev.defcuts.3122.set2.sgnl_1.root");
  SetFileName("/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.12agev.defcuts.3122.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.urqmd.12agev.defcuts.3122.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.3.3agev.defcuts.3122.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.master.root");
  std::string particle = "#Lambda";
  
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.12agev.defcuts.3312.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.urqmd.12agev.defcuts.3312.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.3.3agev.defcuts.3312.set1.sgnl_1.root");
//   std::string particle = "#Xi^{-}";
  
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.12agev.defcuts.310.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.urqmd.12agev.defcuts.310.set1.sgnl_1.root");
//   SetFileName("/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.3.3agev.defcuts.310.set1.sgnl_1.root");
//   std::string particle = "K^{0}_{S}";
  
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
    
  TFile* fileIn = TFile::Open(fileName_.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>>("sim/u_sim_PLAIN.Q_psi_PLAIN.x1x1");
//   auto* dc = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/usimPsi/v1.u_sim_PLAIN.Q_psi_PLAIN.x1x1");
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }
  
//   SetSelectAxisBinEdges({0, 70});  
//   IntegrateSelectAxis();

//   IntegrateSliceAxis();
//   SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});
//   SetSliceAxisBinEdges({0, 10, 20, 40, 70});

//   SetProjectionAxisBinEdges({0.1-axes.at(kProjection).shift_,
//                              0.3-axes.at(kProjection).shift_,
//                              0.5-axes.at(kProjection).shift_,
//                              0.7-axes.at(kProjection).shift_,
//                              0.9-axes.at(kProjection).shift_,
//                              1.1-axes.at(kProjection).shift_});
      
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
    
  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){

    auto v1_rec = DoubleDifferentialCorrelation( fileName_.c_str(),
                                                {("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN." + component_1).c_str(),
                                                 ("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN." + component_2).c_str() } );
//     auto v1_rec = DoubleDifferentialCorrelation( fileName_.c_str(),
//                                                 {("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd1_RECENTERED.res_sub4_sts_pipos." + component_1).c_str(),
//                                                  ("v1/uQ_R1_sub4_sts_pipos/v1.u_rec_RESCALED.psd1_RECENTERED.res_sub4_sts_pipos." + component_2).c_str() } );
    
    v1_rec.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1_rec.Scale(2);
    v1_rec.SetMarker(kFullSquare);
//     v1_rec.SetPalette( {kRed, kGreen+2, kBlue, kViolet, kBlack } );
//     v1_rec.SetPalette( {kOrange+1, kGreen+2, kBlue, kRed, kViolet, kAzure-4, kBlack, kMagenta } );
    v1_rec.SetPalette({kOrange+1, kSpring, kGreen+2, kAzure-4, kRed, kViolet, kBlue});
    v1_rec.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_rec.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_rec.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_rec.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_rec.Calculate();
    v1_rec.ShiftProjectionAxis(axes.at(kProjection).shift_);

    auto v1_sim = DoubleDifferentialCorrelation( fileName_.c_str(),
                                                {("sim/u_sim_PLAIN.Q_psi_PLAIN." + component_1).c_str(),
                                                 ("sim/u_sim_PLAIN.Q_psi_PLAIN." + component_2).c_str() } );
    v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1_sim.Scale(2);
    v1_sim.SetMarker(-1);
//     v1_sim.SetPalette( {kOrange+1, kGreen+2, kBlue, kRed, kViolet, kAzure-4, kBlack, kMagenta } );
//     v1_sim.SetPalette( {kRed-2, kGreen-2, kBlue-2, kViolet-2, kBlack-2 } );
    v1_sim.SetPalette({kOrange+1, kSpring, kGreen+2, kAzure-4, kRed, kViolet, kBlue});
    v1_sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_sim.Calculate();
    v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);

    HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});
    
//     pic.AddText({0.2, 0.90, "Au+Au, DCM-QGSM-SMM, 12A GeV/c"}, 0.025);
    pic.AddText({0.2, 0.90, particle.c_str()}, 0.035);
    pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
    pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
//     pic.AddText({0.2, 0.84, "UrQMD"}, 0.025);
    pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);
//     pic.AddText({0.2, 0.81, "3.3A GeV/c"}, 0.025);

    pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + 
                            " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
//     pic.AddText({0.2, 0.75, "#Psi(PSD_{1})"}, 0.025);
    
    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);
    leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());
    
    for( auto obj : v1_rec.GetProjections() ){
      pic.AddDrawable( obj );
      leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
    }
    for( auto obj : v1_sim.GetProjections() ){
//       obj->SetErrorOption("E3");
      pic.AddDrawable( obj );
    }

    pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), ("v_{" + harmonic + "}").c_str()});
    

//     pic.SetXRange({-0.05, 0.65});
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
  
  
//   pic.Save("./v1_au123_urqmd_pi_neg_pT_y","png");
//   pic.Save("./v1_au123_urqmd_pi_neg_pT_y","root");
}


