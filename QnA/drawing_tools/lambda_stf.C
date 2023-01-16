#include "lambda.h"

void lambda_stf() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".root";
  
  std::vector<std::string> particles{"lambda", "kshort", "pipos", "pineg"};
  std::vector<std::string> subevents{"psd1", "psd2", "psd3", "spec1_prim", "spec2_prim", "spec3_prim"};
  std::string step;
  bool average_comp;

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");

  std::string uQ_R1 = "uQ_R1"; std::string y_axis_title = "v_{1}";
//   std::string uQ_R1 = "uQ_R1_even"; std::string y_axis_title = "v_{1}^{even}";

  std::vector<std::string> components{"x1x1", "y1y1"};
  
  axes.at(0).sim_name_ = "SimParticles_pT";
  axes.at(1).sim_name_ = "SimParticles_rapidity";
  axes.at(2).sim_name_ = "SimEventHeader_centrality_impactpar";
  axes.at(0).reco_name_ = "SimParticles_pT";
  axes.at(1).reco_name_ = "SimParticles_rapidity";
  axes.at(2).reco_name_ = "SimEventHeader_centrality_impactpar";
  
  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/lambda/uPsi/v1.uPsi.x1x1");
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }
  
//   SetProjectionAxisBinEdges({-1.0-axes.at(kProjection).shift_,
//                              -0.6-axes.at(kProjection).shift_,
//                              -0.2-axes.at(kProjection).shift_,
//                               0.2-axes.at(kProjection).shift_,
//                               0.6-axes.at(kProjection).shift_,
//                               1.0-axes.at(kProjection).shift_});

  for(auto& particle : particles) {
    if(evegen == "dcmqgsm") {
      if(particle == "lambda") SetSliceAxisBinEdges({0, 0.4, 0.8, 1.2, 1.6});
      if(particle == "kshort") SetSliceAxisBinEdges({0, 0.4, 0.8, 1.6});
      if(particle == "pipos" || particle == "pineg") SetSliceAxisBinEdges({0, 0.4, 0.6, 1.0, 1.4, 2.0});
    }

    if(evegen == "urqmd") {
      if(particle == "lambda") SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});
      if(particle == "kshort") SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});
      if(particle == "pipos" || particle == "pineg") SetSliceAxisBinEdges({0, 1.0, 1.4, 2.0});
    }

    for(auto& subevent : subevents) {
      if(subevent[0] == 'p') {
        step = "_RECENTERED";
        average_comp = false;
      }
      if(subevent[0] == 's') {
        step = "_PLAIN";
        average_comp = true;
      }

      bool is_first_canvas = true;
      std::string fileOutName;
      if(uQ_R1 == "uQ_R1") fileOutName = "v1_res." + particle + "." + subevent;
      if(uQ_R1 == "uQ_R1_even") fileOutName = "v1_res_even." + particle + "." + subevent;
      //       TFile* fileOut = TFile::Open("fileOut.root", "recreate");

      for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){

        std::vector<DoubleDifferentialCorrelation> v1_R_MC;
        if(average_comp) {
          v1_R_MC.resize(1);
          v1_R_MC.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/v1." + uQ_R1 + "." + subevent + step + ".x1x1").c_str(),
                                                                            ("v1/" + particle + "/uQ_R1/v1." + uQ_R1 + "." + subevent + step + ".y1y1").c_str()} );
          v1_R_MC.at(0).SetMarker(kFullSquare);
        } else {
          v1_R_MC.resize(2);
          v1_R_MC.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/v1." + uQ_R1 + "." + subevent + step + ".x1x1").c_str()} );
          v1_R_MC.at(1) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/v1." + uQ_R1 + "." + subevent + step + ".y1y1").c_str()} );
          v1_R_MC.at(0).SetMarker(kFullSquare);
          v1_R_MC.at(1).SetMarker(kOpenSquare);
        }

        for(auto& vc : v1_R_MC) {
          vc.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
          vc.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                         kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
          vc.SetBiasPalette(false);
          vc.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
          vc.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
          vc.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
          vc.ShiftSliceAxis(axes.at(kSlice).shift_);
          vc.Calculate();
          vc.ShiftProjectionAxis(axes.at(kProjection).shift_);
        }

        v1_R_MC.at(0).SlightShiftProjectionAxis(0.025);
        if(!average_comp) {
          v1_R_MC.at(1).SlightShiftProjectionAxis(0.025, 0.0125);
        }

        auto v1_PsiRP = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uPsi/v1.uPsi.x1x1").c_str(),
                                                                          ("v1/" + particle + "/uPsi/v1.uPsi.y1y1").c_str()} );
        v1_PsiRP.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_PsiRP.SetMarker(-1);
        v1_PsiRP.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                            kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
        v1_PsiRP.SetBiasPalette(false);
        v1_PsiRP.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                        {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
        v1_PsiRP.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_PsiRP.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_PsiRP.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_PsiRP.Calculate();
        v1_PsiRP.ShiftProjectionAxis(axes.at(kProjection).shift_);

        HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});
        pic.SetRelErrorThreshold(0.2);

        pic.AddText({0.2, 0.90, particle.c_str()}, 0.025);
        if(evegen == "dcmqgsm") {
          pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
          pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
        }
        if(evegen == "urqmd") {
          pic.AddText({0.2, 0.87, "2M Au+Au"}, 0.025);
          pic.AddText({0.2, 0.84, "UrQMD"}, 0.025);
        }
        pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);
        pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                                " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
        pic.AddText({0.2, 0.75, subevent.c_str()}, 0.025);

        auto leg1 = new TLegend();
        leg1->SetBorderSize(1);
        leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

        TLegendEntry* entry;
        entry = leg1->AddEntry("", "#Psi^{RP}, x&y averaged", "L");
        entry->SetMarkerSize(2);
        if(average_comp) {
          entry = leg1->AddEntry("", "#Psi(R_{1}^{true}), x&y averaged", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kFullSquare);
        } else {
          entry = leg1->AddEntry("", "#Psi(R_{1}^{true}), X", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kFullSquare);
          entry = leg1->AddEntry("", "#Psi(R_{1}^{true}), Y", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kOpenSquare);
        }

        for(int i=0; i<v1_R_MC.size(); i++) {
          for( auto obj : v1_R_MC.at(i).GetProjections() ){
            pic.AddDrawable( obj );
            if(i==0) leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
          }
        }
        for( auto obj : v1_PsiRP.GetProjections() ){
          pic.AddDrawable( obj );
        }
        pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title});

        pic.CustomizeXRange();
        pic.CustomizeYRange();
        pic.AddLegend(leg1);
        pic.CustomizeLegend(leg1);
        pic.Draw();

//           fileOut->cd();
//           pic.GetCanvas()->Write();

        if(is_first_canvas)
          pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
        else
          pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
        is_first_canvas = false;

//         if(average_comp) break;
      }

      TCanvas emptycanvas("", "", 1000, 1000);
      emptycanvas.Print((fileOutName + ".pdf)").c_str(), "pdf");

//       fileOut->Close();
    }
  }
}
