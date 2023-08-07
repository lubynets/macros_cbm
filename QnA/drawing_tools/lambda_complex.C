#include "lambda.h"

void lambda_complex() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3"; axes.at(1).shift_ = -0.985344;
//   std::string evegen = "urqmd";   std::string pbeam = "12";

//   std::string particle = "#Lambda"; std::string pdg = "3122";
  std::string particle = "K^{0}_{S}"; std::string pdg = "310";
// //   std::string particle = "#Xi^{-}"; std::string pdg = "3312";

  std::string cuts = "lc1";
  if(pbeam == "3.3") cuts = "oc1";
  if(pdg == "3312") cuts = "dc";

  bool is_imf = true;
  if(pdg=="3312" || (pdg=="3122" && pbeam=="3.3")) is_imf = false;

//   bool is_write_rootfile = false;
  bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileMcName = "/home/oleksii/cbmdir/working/qna/aXmass/vR." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";
  std::string fileRecName = "/home/oleksii/cbmdir/working/qna/aXmass/of." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");
  std::string component_1 = "x1x1";
  std::string component_2 = "y1y1";

  if(!is_imf) {
    fileRecName = fileMcName;
    axes.at(0).reco_fit_name_ = axes.at(0).reco_name_;
    axes.at(1).reco_fit_name_ = axes.at(1).reco_name_;
    axes.at(2).reco_fit_name_ = axes.at(2).reco_name_;
  }

  struct Inputs {
    std::string dirname_;
    std::vector<std::string> subevents_;
    std::string resname_;
  };

  std::vector<Inputs> inputs {
//     {"uPsi", {"Q_psi"}, ""},
//     {"uQ_R1_MC", {"psd1", "psd2", "psd3"}, "_res_MC"},
//     {"uQ_R1_sub3", {"psd1", "psd2", "psd3"}, "_res_sub3"},
    {"uQ_R1_sub4", {/*"psd1", */"psd2"/*, "psd3"*/}, "_res_sub4_sts_pipos"},
  };

  TFile* fileMc = TFile::Open(fileMcName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileMc->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/usimPsi/v1.u_sim.Q_psi.x1x1");
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

//   SetSelectAxisBinEdges({10, 30});
//   SetSliceAxisBinEdges({0.2, 1.0});
//   IntegrateSliceAxis();
//   SetProjectionAxisBinEdges({-0.5-axes.at(kProjection).shift_,
//                              -0.1-axes.at(kProjection).shift_,
//                               0.1-axes.at(kProjection).shift_,
//                               0.5-axes.at(kProjection).shift_,
//                               0.7-axes.at(kProjection).shift_});

  TFile* fileOut;

  for(auto& ip : inputs) {
    for(auto& se : ip.subevents_) {

      bool is_first_canvas = true;
      std::string fileOutName = "v1." + evegen + "." + pbeam + "agev." + pdg + "." + se + ip.resname_;

      if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

      for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

        auto v1_sim = DoubleDifferentialCorrelation( fileMcName.c_str(),
                                                    {("v1/usimPsi/v1.u_sim.Q_psi." + component_1).c_str(),
                                                     ("v1/usimPsi/v1.u_sim.Q_psi." + component_2).c_str() } );
        v1_sim.SetErrorType(error_mode);
        v1_sim.SetMeanType(mean_mode);
        v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_sim.SetMarker(-1);
        v1_sim.SetIsFillLine();
        v1_sim.SetPalette(palette1);
        v1_sim.SetBiasPalette(false);
        v1_sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                      {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
        v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_sim.Calculate();

        std::string rec_nofit_name = "v1/" + ip.dirname_ + "/v1.u_rec_sgnl." + se + ip.resname_ + ".";
        auto v1_rec_nofit = DoubleDifferentialCorrelation( fileMcName.c_str(),
                                                    {(rec_nofit_name + component_1).c_str(),
                                                     (rec_nofit_name + component_2).c_str() } );
        v1_rec_nofit.SetErrorType(error_mode);
        v1_rec_nofit.SetMeanType(mean_mode);
        v1_rec_nofit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_rec_nofit.SetMarker(kOpenSquare);
        v1_rec_nofit.SetPalette(palette1);
        v1_rec_nofit.SetBiasPalette(false);
        v1_rec_nofit.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                            {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
        v1_rec_nofit.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_rec_nofit.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_rec_nofit.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_rec_nofit.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_rec_nofit.SlightShiftProjectionAxis(0.025, 0.0125);
        v1_rec_nofit.Calculate();

        std::string rec_fit_name;
        if(is_imf) rec_fit_name = "Pars/" + ip.dirname_ + "/" + se + ip.resname_ + "/signal.";
        else       rec_fit_name = "v1/" + ip.dirname_ + "/v1.u_rec." + se + ip.resname_ + ".";

        auto v1_rec_fit = DoubleDifferentialCorrelation( fileRecName.c_str(),
                                                    {(rec_fit_name + component_1).c_str(),
                                                     (rec_fit_name + component_2).c_str() } );
        v1_rec_fit.SetErrorType(error_mode);
        v1_rec_fit.SetMeanType(mean_mode);
        v1_rec_fit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_rec_fit.SetMarker(kFullSquare);
        v1_rec_fit.SetPalette(palette1);
        v1_rec_fit.SetBiasPalette(false);
        v1_rec_fit.Rebin({{axes.at(kSelect).reco_fit_name_.c_str(),
                          {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
        v1_rec_fit.SetProjectionAxis({axes.at(kProjection).reco_fit_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_rec_fit.SetSliceAxis({axes.at(kSlice).reco_fit_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_rec_fit.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_rec_fit.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_rec_fit.SlightShiftProjectionAxis(0.025);
        v1_rec_fit.Calculate();

        HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});

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
        pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                                " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
        pic.AddText({0.2, 0.75, se.c_str()}, 0.025);
        if(ip.resname_ != "") pic.AddText({0.2, 0.72, ip.resname_.substr(1, ip.resname_.size()).c_str()}, 0.025);

        auto leg1 = new TLegend();
        leg1->SetBorderSize(1);

        auto* entry = leg1->AddEntry("", "MC input", "F");
        entry->SetFillColorAlpha(kBlack, 0.2);
        entry->SetLineColor(kWhite);
        entry->SetFillStyle(1000);

//         entry = leg1->AddEntry("", "REC, MC-match", "P");
//         entry->SetMarkerSize(2);
//         entry->SetMarkerStyle(kOpenSquare);

        entry = leg1->AddEntry("", "Reconstructed", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kFullSquare);

        leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

        for( auto obj : v1_sim.GetProjections() ){
          pic.AddDrawable( obj );
        }
//         for( auto obj : v1_rec_nofit.GetProjections() ){
//           pic.AddDrawable( obj );
//         }
        for( auto obj : v1_rec_fit.GetProjections() ){
          pic.AddDrawable( obj );
          leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
        }


        pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), "v_{1}"});

    //     pic.SetXRange({-0.05, 0.95});
        pic.CustomizeXRange();
        pic.CustomizeYRangeWithLimits(-0.3, 0.6);
        pic.AddLegend(leg1);
        pic.CustomizeLegend(leg1);
//         pic.SetGridXY();
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

      }

      TCanvas emptycanvas("", "", 1000, 1000);
      emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

      if(is_write_rootfile) fileOut->Close();
    }
  }
}
