#include "lambda.h"

void lambda_rtpshort() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3"; axes.at(1).shift_ = -0.985344;
  std::string evegen = "urqmd";   std::string pbeam = "12";

  std::string particle = "#Lambda"; std::string pdg = "3122";
//   std::string particle = "K^{0}_{S}"; std::string pdg = "310";
// //   std::string particle = "#Xi^{-}"; std::string pdg = "3312";

  std::string cuts = "lc1";
  if(pbeam == "3.3") cuts = "oc1";
  if(pdg == "3312") cuts = "dc";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileInName = "/home/oleksii/cbmdir/working/qna/aXmass/vR." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";
  std::string fileDv1DyName = "/home/oleksii/cbmdir/working/qna/aXmass/vR.dv1dy." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");
  std::string component_1 = "x1x1";
  std::string component_2 = "y1y1";

  struct Inputs {
    std::string dirname_;
    std::vector<std::string> subevents_;
    std::string resname_;
  };

  std::vector<Inputs> inputs {
    {"uPsi", {"Q_psi"}, ""},
//     {"uQ_R1_MC", {"psd1", "psd2", "psd3"}, "_res_MC"},
//     {"uQ_R1_sub3", {"psd1", "psd2", "psd3"}, "_res_sub3"},
//     {"uQ_R1_sub4", {"psd1", "psd2", "psd3"}, "_res_sub4_sts_pipos"},
  };

  TFile* fileMc = TFile::Open(fileInName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileMc->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("v1/usimPsi/v1.u_sim.Q_psi.x1x1");
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  TFile* fileDv1Dy = TFile::Open(fileDv1DyName.c_str(), "open");
  if(fileDv1Dy == nullptr)  throw std::runtime_error("fileDv1Dy == nullptr");

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
      std::string fileOutName = "v1.rtp." + evegen + "." + pbeam + "agev." + pdg + "." + se + ip.resname_;

      if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

      auto* sim_slope = (QnDcSD*)fileDv1Dy->Get<QnDcSD>("v1/usimPsi/slope/v1.u_sim.Q_psi.ave");
      auto* sim_offset = (QnDcSD*)fileDv1Dy->Get<QnDcSD>("v1/usimPsi/intercept/v1.u_sim.Q_psi.ave");
      if(sim_slope==nullptr || sim_offset==nullptr) throw std::runtime_error("sim_slope==nullptr || sim_offset==nullptr");

      auto* rec_slope = (QnDcSD*)fileDv1Dy->Get<QnDcSD>(("v1/" + ip.dirname_ + "/slope/v1.u_rec_sgnl." + se + ip.resname_ + ".ave").c_str());
      auto* rec_offset = (QnDcSD*)fileDv1Dy->Get<QnDcSD>(("v1/" + ip.dirname_ + "/intercept/v1.u_rec_sgnl." + se + ip.resname_ + ".ave").c_str());
      if(rec_slope==nullptr || rec_offset==nullptr) throw std::runtime_error("rec_slope==nullptr || rec_offset==nullptr");

      bool select_then_slice;
      if(sim_slope->GetAxes().at(0).Name() == axes.at(kSelect).sim_name_.c_str() &&
         sim_slope->GetAxes().at(1).Name() == axes.at(kSlice).sim_name_.c_str()) {
        select_then_slice = true;
      } else if (sim_slope->GetAxes().at(1).Name() == axes.at(kSelect).sim_name_.c_str() &&
                 sim_slope->GetAxes().at(0).Name() == axes.at(kSlice).sim_name_.c_str()) {
        select_then_slice  = false;
      } else throw std::runtime_error("Something wrong with Select and Slice axes in dv1dy");

      for(unsigned long iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

        std::string sim_name = "v1/usimPsi/v1.u_sim.Q_psi.";
        auto v1_sim = DoubleDifferentialCorrelation( fileInName.c_str(),
                                                    {(sim_name + component_1).c_str(),
                                                     (sim_name + component_2).c_str() } );
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
        auto v1_rec = DoubleDifferentialCorrelation( fileInName.c_str(),
                                                    {(rec_nofit_name + component_1).c_str(),
                                                     (rec_nofit_name + component_2).c_str() } );
        v1_rec.SetErrorType(error_mode);
        v1_rec.SetMeanType(mean_mode);
        v1_rec.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_rec.SetMarker(kFullSquare);
        v1_rec.SetPalette(palette1);
        v1_rec.SetBiasPalette(false);
        v1_rec.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                      {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
        v1_rec.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_rec.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_rec.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_rec.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_rec.SlightShiftProjectionAxis(0.025);
        v1_rec.Calculate();

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

        std::string ref_title = "MC-true tracks";
        std::string est_title = "Reconstructed tracks";

        TLegendEntry* entry;
        entry = leg1->AddEntry("", ref_title.c_str(), "F");
        entry->SetFillColorAlpha(kBlack, 0.2);
        entry->SetLineColor(kWhite);
        entry->SetFillStyle(1000);

        entry = leg1->AddEntry("", est_title.c_str(), "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kFullSquare);

        leg1->AddEntry("", (axes.at(kSlice).title_+axes.at(kSlice).unit_ + ":").c_str(), "");

        for( auto obj : v1_sim.GetProjections() ){
          pic.AddDrawable( obj );
        }
        for( auto obj : v1_rec.GetProjections() ){
          pic.AddDrawable( obj );
          leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
        }
        for(unsigned long iSlice=0; iSlice<v1_sim.GetProjections().size(); iSlice++) {
          float k, b;
          if(select_then_slice) {
            k = sim_slope->At({iEdge, iSlice}).Mean();
            b = sim_offset->At({iEdge, iSlice}).Mean();
          } else {
            k = sim_slope->At({iSlice, iEdge}).Mean();
            b = sim_offset->At({iSlice, iEdge}).Mean();
          }
          TF1* fsim = new TF1("fsim", "[0]+[1]*x", -10, 10);
          fsim->SetParameters(b, k);
          fsim->SetLineStyle(2);
          fsim->SetLineColor(palette1.at(iSlice));
          pic.AddFunction(fsim);

          if(select_then_slice) {
            k = rec_slope->At({iEdge, iSlice}).Mean();
            b = rec_offset->At({iEdge, iSlice}).Mean();
          } else {
            k = rec_slope->At({iSlice, iEdge}).Mean();
            b = rec_offset->At({iSlice, iEdge}).Mean();
          }
          TF1* frec = new TF1("frec", "[0]+[1]*x", -10, 10);
          frec->SetParameters(b, k);
          frec->SetLineStyle(1);
          frec->SetLineColor(palette1.at(iSlice));
          pic.AddFunction(frec);
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
