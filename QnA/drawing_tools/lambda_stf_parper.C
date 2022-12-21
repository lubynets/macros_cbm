#include "lambda.h"

bool float_equality(float a, float b) {
  if (std::fabs(a-b)<1e-5)
    return true;
  else
    return false;
}

float C_via_b(float b) {
  if (float_equality(b, 0.0))  return 0.;
  if (float_equality(b, 4.2))  return 10.;
  if (float_equality(b, 6.0))  return 20.;
  if (float_equality(b, 7.3))  return 30.;
  if (float_equality(b, 8.4))  return 40.;
  if (float_equality(b, 9.4))  return 50.;
  if (float_equality(b, 10.3)) return 60.;
  if (float_equality(b, 11.2)) return 70.;
  return -1;
}

void lambda_stf_parper() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );
  bool is_first_canvas = true;
    
  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/v1andR1.stf.root";
  
//   std::string particle = "lambda";
//   std::string particle = "pipos";
  std::string particle = "pineg";
  
  std::string subevent = "Q_spec_PLAIN";
  
  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");
  std::string component_1 = "x1x1";
  std::string component_2 = "y1y1";
  std::string harmonic = "1";
  
  axes.at(0).sim_name_ = "SimParticles_pT";
  axes.at(1).sim_name_ = "SimParticles_rapidity";
  axes.at(2).sim_name_ = "SimEventHeader_b";
  axes.at(0).reco_name_ = "SimParticles_pT";
  axes.at(1).reco_name_ = "SimParticles_rapidity";
  axes.at(2).reco_name_ = "SimEventHeader_b";
  
  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/lambda/uPsi/v1.uPsi.x1x1");
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }
  
//   SetSliceAxisBinEdges({0, 0.4, 0.8, 1.2, 1.6});    // lambda
//   SetSliceAxisBinEdges({0, 0.6, 1.0, 1.4, 2.0});
  SetSliceAxisBinEdges({0, 0.4, 0.6, 1.0, 1.4, 2.0});   // pipos pineg
  
//   SetProjectionAxisBinEdges({-1.0-axes.at(kProjection).shift_,
//                              -0.6-axes.at(kProjection).shift_,
//                              -0.2-axes.at(kProjection).shift_,
//                               0.2-axes.at(kProjection).shift_,
//                               0.6-axes.at(kProjection).shift_,
//                               1.0-axes.at(kProjection).shift_});
  
  TF1 oneline("oneline", "1", -100, 100);
  TF1 minusoneline("minusoneline", "-1", -100, 100);  
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
    
  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){

    auto v1_R_MC = DoubleDifferentialCorrelation( fileName.c_str(),
                                                {("v1/" + particle + "/uQ_R1/v1.uQ_R1." + subevent + "." + component_1).c_str(),
                                                 ("v1/" + particle + "/uQ_R1/v1.uQ_R1." + subevent + "." + component_2).c_str()} );
    v1_R_MC.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1_R_MC.Scale(2);
    v1_R_MC.SetMarker(kFullSquare);
    v1_R_MC.SetPalette({kOrange+1, kSpring, kGreen+2, kAzure-4, kRed, kViolet, kBlue,
                        kOrange+1, kSpring, kGreen+2, kAzure-4, kRed, kViolet, kBlue});
    v1_R_MC.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                   {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_R_MC.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_R_MC.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_R_MC.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_R_MC.Calculate();
    v1_R_MC.ShiftProjectionAxis(axes.at(kProjection).shift_);

    auto v1_PsiRP = DoubleDifferentialCorrelation( fileName.c_str(),
                                                {("v1/" + particle + "/uPsi/v1.uPsi." + component_1).c_str(),
                                                 ("v1/" + particle + "/uPsi/v1.uPsi." + component_2).c_str() } );
    v1_PsiRP.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1_PsiRP.Scale(2);
    v1_PsiRP.SetMarker(-1);
    v1_PsiRP.SetPalette({kOrange+1, kSpring, kGreen+2, kAzure-4, kRed, kViolet, kBlue,
                         kOrange+1, kSpring, kGreen+2, kAzure-4, kRed, kViolet, kBlue});
    v1_PsiRP.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_PsiRP.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_PsiRP.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_PsiRP.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_PsiRP.Calculate();
    v1_PsiRP.ShiftProjectionAxis(axes.at(kProjection).shift_);

    HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});
    
    pic.AddText({0.2, 0.90, particle.c_str()}, 0.025);
//     pic.AddText({0.2, 0.90, "Au+Au, DCM-QGSM-SMM, 12A GeV/c"}, 0.025);
    pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
    pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
    pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);

//     pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + 
//                             " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
    pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(C_via_b(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_), axes.at(kSelect).precision_) + 
                            " - " + to_string_with_precision(C_via_b(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_), axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
    pic.AddText({0.2, 0.75, subevent.c_str()}, 0.025);
    
    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);
    leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());
    
    if(!is_draw_difference) {    
      for( auto obj : v1_R_MC.GetProjections() ){
        pic.AddDrawable( obj );
        leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
      }
      for( auto obj : v1_PsiRP.GetProjections() ){
        pic.AddDrawable( obj );
      }
      pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), ("v_{" + harmonic + "}").c_str()});
    }
    else {    
      pic.DrawZeroLine(false);
      int i{0};
      for( auto obj_R_MC : v1_R_MC.GetProjections() ) {
        auto obj_psiRP = v1_PsiRP.GetProjections().at(i);
        GraphSubtractor gr_sub;
        gr_sub.SetMinuend(obj_R_MC);
        gr_sub.SetSubtrahend(obj_psiRP);
        gr_sub.Calculate();
        auto obj_R_MC_psiRP = gr_sub.GetResult();
        pic.AddDrawable(obj_R_MC_psiRP);
        
        leg1->AddEntry( obj_R_MC->GetPoints(), obj_R_MC->GetTitle().c_str(), "P" );
        
        i++;
      }
      pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(),
                        ("(v_{" + harmonic + "}{#Psi^{" + subevent + "}(R_{1}^{true})} - v_{" + harmonic + "}{#Psi^{RP}})  /"+
                        "  #sqrt{#sigma^{2}_{v_{" + harmonic + "}(#Psi^{" + subevent + "}(R_{1}^{true}))} - #sigma^{2}_{v_{" + harmonic + "}(#Psi^{RP})}}").c_str()});
    }    
    
    pic.CustomizeXRange();
    pic.CustomizeYRange();
    pic.AddLegend(leg1);
    pic.CustomizeLegend(leg1);
    pic.Draw();
    
    if(is_draw_difference) {
      oneline.Draw("same");
      minusoneline.Draw("same");
    }    
    
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
