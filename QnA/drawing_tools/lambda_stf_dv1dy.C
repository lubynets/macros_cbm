#include "lambda.h"

void lambda_stf_dv1dy() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

//   std::string evegen = "dcmqgsm";
  std::string evegen = "urqmd";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/dv1dy.stf." + evegen + ".root";

  std::vector<std::string> particles{
                                     "lambda",
                                     "kshort",
                                     "xi",
                                     "pipos",
                                     "pineg"
                                    };
  std::vector<std::string> subevents{"psd1", "psd2", "psd3",
                                     "etacut_1_charged", "etacut_2_charged", "etacut_3_charged",
                                     "etacut_1_all", "etacut_2_all", "etacut_3_all"};
  std::string step;
  bool average_comp;
  float y_lo, y_hi;

//   SetAxis("centrality", "select");
  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title;

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::vector<std::string> fitcoeffs{"slope", "intercept"};

  axes.at(kSlice).sim_name_ = "SimParticles_pT";
  axes.at(kProjection).sim_name_ = "SimEventHeader_centrality_impactpar";
  axes.at(kSlice).reco_name_ = "SimParticles_pT";
  axes.at(kProjection).reco_name_ = "SimEventHeader_centrality_impactpar";

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");

  for(auto& particle : particles) {

    auto* dc = (Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>>((particle + "/v1sim_intercept.psi.ave").c_str());
    if(dc==nullptr) {
      throw std::runtime_error("DataContainer is nullptr");
    };
    for(auto& ax : axes) {
      if(ax.name_ == "rapidity") continue;
      Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
      ax.bin_edges_.clear();
      for(int i=0; i<=qnaxis.size(); i++) {
        ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
      }
    }
    delete dc;

    for(auto& subevent : subevents) {
      if(subevent[0] == 'p') {
        step = "_RECENTERED";
        average_comp = false;
      }
      if(subevent[0] == 'e') {
        step = "_PLAIN";
        average_comp = true;
      }

      bool is_first_canvas = true;
      std::string fileOutName = "dv1dy." + particle + "." + subevent;
      //       TFile* fileOut = TFile::Open("fileOut.root", "recreate");


      for(auto fc : fitcoeffs) {

        if(evegen == "dcmqgsm") {
          if(particle == "lambda") {
            if(fc == "slope")     { y_lo = -0.2; y_hi = 0.4; }
            if(fc == "intercept") { y_lo = -0.1; y_hi = 0.05; }
          }
          if(particle == "kshort") {
            if(fc == "slope")     { y_lo = -0.1; y_hi = 0.1; }
            if(fc == "intercept") { y_lo = -0.05; y_hi = 0.02; }
          }
          if(particle == "pipos" || particle == "pineg") {
            if(fc == "slope")     { y_lo = -0.2; y_hi = 0.4; }
            if(fc == "intercept") { y_lo = -0.1; y_hi = 0.1; }
          }
          if(particle == "xi") {
            if(fc == "slope")     { y_lo = -0.2; y_hi = 0.4; }
            if(fc == "intercept") { y_lo = -0.1; y_hi = 0.1; }
          }
        }
        if(evegen == "urqmd") {
          if(particle == "lambda") {
            if(fc == "slope")     { y_lo = -0.15; y_hi = 0.1; }
            if(fc == "intercept") { y_lo = -0.1; y_hi = 0.1; }
          }
          if(particle == "kshort") {
            if(fc == "slope")     { y_lo = -0.1; y_hi = 0.3; }
            if(fc == "intercept") { y_lo = -0.1; y_hi = 0.1; }
          }
          if(particle == "pipos" || particle == "pineg") {
            if(fc == "slope")     { y_lo = -0.2; y_hi = 0.1; }
            if(fc == "intercept") { y_lo = -0.15; y_hi = 0.05; }
          }
          if(particle == "xi") {
            if(fc == "slope")     { y_lo = -0.2; y_hi = 0.1; }
            if(fc == "intercept") { y_lo = -0.1; y_hi = 0.1; }
          }
        }


        if(fc == "slope") y_axis_title = "dv_{1}/dy";
        if(fc == "intercept") y_axis_title = "v_{1}|_{y=0}";

        std::vector<DoubleDifferentialCorrelation> v1_R_MC;
        if(average_comp) {
          v1_R_MC.resize(1);
          v1_R_MC.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + subevent + ".ave").c_str()} );
          v1_R_MC.at(0).SetMarker(kFullSquare);
        } else {
          v1_R_MC.resize(2);
          v1_R_MC.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + subevent + ".x1x1").c_str()} );
          v1_R_MC.at(1) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + subevent + ".y1y1").c_str()} );
          v1_R_MC.at(0).SetMarker(kFullSquare);
          v1_R_MC.at(1).SetMarker(kOpenSquare);
        }

        for(auto& vc : v1_R_MC) {
          vc.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
          vc.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                         kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
          vc.SetBiasPalette(false);
          vc.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
          vc.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
          vc.ShiftSliceAxis(axes.at(kSlice).shift_);
          vc.Calculate();
          vc.ShiftProjectionAxis(axes.at(kProjection).shift_);
        }

        v1_R_MC.at(0).SlightShiftProjectionAxis(1);
        if(!average_comp) {
          v1_R_MC.at(1).SlightShiftProjectionAxis(1, 0.5);
        }

//         auto v1_R_MC = DoubleDifferentialCorrelation( fileName.c_str(), correlnames );
//         v1_R_MC.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
//         v1_R_MC.SetMarker(kFullSquare);
//         v1_R_MC.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
//                             kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
//         v1_R_MC.SetBiasPalette(false);
//         v1_R_MC.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
//         v1_R_MC.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
//         v1_R_MC.ShiftSliceAxis(axes.at(kSlice).shift_);
//         v1_R_MC.Calculate();
//         v1_R_MC.ShiftProjectionAxis(axes.at(kProjection).shift_);
//         v1_R_MC.SlightShiftProjectionAxis(1);
//
        auto v1_PsiRP = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1sim_" + fc + ".psi.ave").c_str()} );
        v1_PsiRP.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_PsiRP.SetMarker(-1);
        v1_PsiRP.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                              kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
        v1_PsiRP.SetBiasPalette(false);
        v1_PsiRP.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_PsiRP.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_PsiRP.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_PsiRP.Calculate();
        v1_PsiRP.ShiftProjectionAxis(axes.at(kProjection).shift_);

        HeapPicture pic(fc, {1000, 1000});
//         pic.SetRelErrorThreshold(0.2);

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
        pic.AddText({0.2, 0.78, subevent.c_str()}, 0.025);

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
//         pic.CustomizeYRange();
        pic.SetYRange({y_lo, y_hi});
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
      }

      TCanvas emptycanvas("", "", 1000, 1000);
      emptycanvas.Print((fileOutName + ".pdf)").c_str(), "pdf");

//       fileOut->Close();
    }
  }
}
