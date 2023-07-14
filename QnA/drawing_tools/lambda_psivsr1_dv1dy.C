#include "lambda.h"

void lambda_psivsr1_dv1dy() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3";
//   std::string evegen = "urqmd";   std::string pbeam = "12";

  std::string particle = "#Lambda"; std::string pdg = "3122"; std::string cuts = "lc1";
//   std::string particle = "K^{0}_{S}"; std::string pdg = "310"; std::string cuts = "oc1";
//   std::string particle = "#Xi^{-}"; std::string pdg = "3312"; std::string cuts = "dc";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileMcName = "/home/oleksii/cbmdir/working/qna/aXmass/vR.dv1dy." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";

  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title;

//   std::vector<std::string> components{"x1x1", "y1y1"};
  std::vector<std::string> components{"ave"};

  std::vector<std::string> fitcoeffs{"slope", "intercept"};

  struct Inputs {
    std::string dirname_;
    std::vector<std::string> subevents_;
    std::string resname_;
  };

  std::vector<std::string> subevents{"psd1", "psd2", "psd3",
                                     /*"etacut_1_charged", "etacut_2_charged", "etacut_3_charged",
                                     "etacut_1_all", "etacut_2_all", "etacut_3_all"*/};

  TFile* fileMc = TFile::Open(fileMcName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>*)fileMc->Get<Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>>("v1/usimPsi/slope/v1.u_sim.Q_psi.x1x1");
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  for(auto& ax : axes) {
    if(ax.name_ == "rapidity") continue;
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  TFile* fileOut;

  for(auto& fc : fitcoeffs) {

    bool is_first_canvas{true};
    std::string fileOutName = fc + "." + evegen + "." + pbeam + "agev." + pdg;
    if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

    if(fc == "slope") y_axis_title = "dv_{1}/dy";
    if(fc == "intercept") y_axis_title = "v_{1}|_{y=0}";

    for(auto& se : subevents) {
      for(auto& co : components) {
        auto v1_ref = DoubleDifferentialCorrelation( fileMcName.c_str(), {("v1/usimPsi/" + fc + "/v1.u_sim.Q_psi." + co).c_str()} );
        v1_ref.SetErrorType(error_mode);
        v1_ref.SetMeanType(mean_mode);
        v1_ref.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_ref.SetMarker(-1);
        v1_ref.SetIsFillLine();
        v1_ref.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                            kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
        v1_ref.SetBiasPalette(false);
        v1_ref.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_ref.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_ref.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_ref.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_ref.Calculate();

        auto v1_est = DoubleDifferentialCorrelation( fileMcName.c_str(), {("v1/usimQ_R1_MC/" + fc + "/v1.u_sim." + se + "_res_MC." + co).c_str()} );
        v1_est.SetErrorType(error_mode);
        v1_est.SetMeanType(mean_mode);
        v1_est.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_est.SetMarker(kFullSquare);
        v1_est.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                           kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
        v1_est.SetBiasPalette(false);
        v1_est.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_est.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_est.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_est.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_est.Calculate();

        HeapPicture pic(fc, {1000, 1000});

        pic.AddText({0.2, 0.90, particle.c_str()}, 0.03);
        if(evegen == "dcmqgsm") {
          if(pbeam == "12") pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
          else              pic.AddText({0.2, 0.87, "5.2M Au+Au"}, 0.025);
          pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
        }
        if(evegen == "urqmd") {
          pic.AddText({0.2, 0.87, "2M Au+Au"}, 0.025);
          pic.AddText({0.2, 0.84, "UrQMD"}, 0.025);
        }
        pic.AddText({0.2, 0.81, (pbeam + "A GeV/c").c_str()}, 0.025);
        pic.AddText({0.2, 0.78, se.c_str()}, 0.025);

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

        for( auto obj : v1_ref.GetProjections() ){
          pic.AddDrawable( obj );
        }
        for( auto obj : v1_est.GetProjections() ){
          pic.AddDrawable( obj );
          leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
        }


        pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title.c_str()});

    //     pic.SetXRange({-0.05, 0.95});
        pic.CustomizeXRange();
        pic.CustomizeYRange();
        pic.AddLegend(leg1);
        pic.CustomizeLegend(leg1);
        pic.Draw();

        if(is_write_rootfile) {
          fileOut->cd();
          pic.GetCanvas()->Write();
        }

        if(is_first_canvas)
          pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
        else
          pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
        is_first_canvas = false;

      }// components
    }// subevents
    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

    if(is_write_rootfile) fileOut->Close();
  }// fitcoeffs
}
