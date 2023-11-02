#include "lambda.h"

void lambda_complex(int iSetup=1, int iPdg=1) {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen, pbeam, pdg, particle;

  if(iSetup==1) {evegen = "dcmqgsm"; pbeam = "12"; }
  if(iSetup==2) {evegen = "dcmqgsm"; pbeam = "3.3"; axes.at(1).shift_ = -0.985344;}
  if(iSetup==3) {evegen = "urqmd";   pbeam = "12"; }

  if(iPdg==1) {pdg = "3122"; particle = "#Lambda";  }
  if(iPdg==2) {pdg = "310";  particle = "K^{0}_{S}";}
  if(iPdg==3) {pdg = "3312"; particle = "#Xi^{-}";  }

  std::string cuts = "lc1";
  if(pbeam == "3.3") cuts = "oc1";
  if(pdg == "3312") cuts = "dc";

  bool is_imf = true;
  if(pdg=="3312" || (pdg=="3122" && pbeam=="3.3")) is_imf = false;

//   pdg += "_finept";

  bool draw_fit{true};
//   bool draw_fit{false};
//   std::string pol = "pol1";
  std::string pol = "pol3";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileMcName = "/home/oleksii/cbmdir/working/qna/aXmass/vR." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";
  std::string fileRecName = "/home/oleksii/cbmdir/working/qna/aXmass/of." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";
  std::string fileMcName_Dv1Dy = "/home/oleksii/cbmdir/working/qna/aXmass/vR.dv1dy_" + pol + "." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";
  std::string fileRecName_Dv1Dy = "/home/oleksii/cbmdir/working/qna/aXmass/of.dv1dy_" + pol + "." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");
  std::string component_1 = "x1x1";
  std::string component_2 = "y1y1";

  if(!is_imf) {
    fileRecName = fileMcName;
    fileRecName_Dv1Dy = fileMcName_Dv1Dy;
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
    {"uPsi", {"Q_psi"}, ""},
    {"uQ_R1_MC", {"psd1", "psd2", "psd3"}, "_res_MC"},
//     {"uQ_R1_sub3", {"psd1", "psd2", "psd3"}, "_res_sub3"},
    {"uQ_R1_sub4", {"psd1", "psd2", "psd3"}, "_res_sub4_sts_pipos"},
  };

  TFile* fileMc = TFile::Open(fileMcName.c_str(), "open");
  if(fileMc == nullptr) throw std::runtime_error("fileMc == nullptr");
  auto* dc = (QnDcSCa*)fileMc->Get<QnDcSCa>("v1/usimPsi/v1.u_sim.Q_psi.x1x1");
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  TFile* fileRec = TFile::Open(fileRecName.c_str(), "open");
  if(fileRec == nullptr) throw std::runtime_error("fileRec == nullptr");

  TFile* fileMc_Dv1Dy = TFile::Open(fileMcName_Dv1Dy.c_str(), "open");
  if(draw_fit && fileMc_Dv1Dy == nullptr) throw std::runtime_error("fileMc_Dv1Dy == nullptr");

  TFile* fileRec_Dv1Dy = TFile::Open(fileRecName_Dv1Dy.c_str(), "open");
  if(draw_fit && fileRec_Dv1Dy == nullptr) throw std::runtime_error("fileRec_Dv1Dy == nullptr");


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

      auto* sim_slope = (QnDcSD*)fileMc_Dv1Dy->Get<QnDcSD>("v1/usimPsi/slope/v1.u_sim.Q_psi.ave");
      auto* sim_offset = (QnDcSD*)fileMc_Dv1Dy->Get<QnDcSD>("v1/usimPsi/intercept/v1.u_sim.Q_psi.ave");
      auto* sim_third = (QnDcSD*)fileMc_Dv1Dy->Get<QnDcSD>("v1/usimPsi/third/v1.u_sim.Q_psi.ave");
      if(draw_fit && (sim_slope==nullptr || sim_offset==nullptr)) throw std::runtime_error("sim_slope==nullptr || sim_offset==nullptr");
      if(draw_fit && pol == "pol3" && sim_third == nullptr) throw std::runtime_error("sim_third == nullptr");

      QnDcSD* rec_slope{nullptr};
      QnDcSD* rec_offset{nullptr};
      QnDcSD* rec_third{nullptr};

      if(is_imf) {
        rec_slope = (QnDcSD*)fileRec_Dv1Dy->Get<QnDcSD>(("Pars/" + ip.dirname_ + "/" + se + ip.resname_ + "/slope/signal.ave").c_str());
        rec_offset = (QnDcSD*)fileRec_Dv1Dy->Get<QnDcSD>(("Pars/" + ip.dirname_ + "/" + se + ip.resname_ + "/intercept/signal.ave").c_str());
        rec_third = (QnDcSD*)fileRec_Dv1Dy->Get<QnDcSD>(("Pars/" + ip.dirname_ + "/" + se + ip.resname_ + "/third/signal.ave").c_str());
      } else {
        rec_slope = (QnDcSD*)fileRec_Dv1Dy->Get<QnDcSD>(("v1/" + ip.dirname_ + "/slope/v1.u_rec." + se + ip.resname_ + ".ave").c_str());
        rec_offset = (QnDcSD*)fileRec_Dv1Dy->Get<QnDcSD>(("v1/" + ip.dirname_ + "/intercept/v1.u_rec." + se + ip.resname_ + ".ave").c_str());
        rec_third = (QnDcSD*)fileRec_Dv1Dy->Get<QnDcSD>(("v1/" + ip.dirname_ + "/third/v1.u_rec." + se + ip.resname_ + ".ave").c_str());
      }
      if(draw_fit && (rec_slope==nullptr || rec_offset==nullptr)) throw std::runtime_error("rec_slope==nullptr || rec_offset==nullptr");
      if(draw_fit && pol == "pol3" && rec_third == nullptr) throw std::runtime_error("rec_third == nullptr");

      bool select_then_slice;
      if(sim_slope->GetAxes().at(0).Name() == axes.at(kSelect).sim_name_.c_str() &&
         sim_slope->GetAxes().at(1).Name() == axes.at(kSlice).sim_name_.c_str()) {
        select_then_slice = true;
      } else if (sim_slope->GetAxes().at(1).Name() == axes.at(kSelect).sim_name_.c_str() &&
                 sim_slope->GetAxes().at(0).Name() == axes.at(kSlice).sim_name_.c_str()) {
        select_then_slice  = false;
      } else throw std::runtime_error("Something wrong with Select and Slice axes in dv1dy");

      for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

        auto v1_sim = DoubleDifferentialCorrelation( fileMcName.c_str(),
                                                    {("v1/usimPsi/v1.u_sim.Q_psi." + component_1).c_str(),
                                                     ("v1/usimPsi/v1.u_sim.Q_psi." + component_2).c_str() } );
        v1_sim.SetErrorType(error_mode);
        v1_sim.SetMeanType(mean_mode);
        v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_sim.SetMarker(-1);
        v1_sim.SetIsFillLine();
        v1_sim.SetPalette(Helper::palette1);
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
        v1_rec_nofit.SetPalette(Helper::palette1);
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
        v1_rec_fit.SetPalette(Helper::palette1);
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

        if(draw_fit) {
          for(unsigned long iSlice=0; iSlice<v1_sim.GetProjections().size(); iSlice++) {
            float k, b, p3;
            const unsigned long i1 = select_then_slice ? iEdge : iSlice;
            const unsigned long i2 = select_then_slice ? iSlice : iEdge;
            k = sim_slope->At({i1, i2}).Mean();
            b = sim_offset->At({i1, i2}).Mean();
            p3 = pol=="pol3" ? sim_third->At({i1, i2}).Mean() : 0;
            TF1* fsim = new TF1("fsim", "[0]+[1]*x+[2]*x*x*x", -10, 10);
            fsim->SetParameters(b, k, p3);
            fsim->SetLineStyle(2);
            fsim->SetLineColor(Helper::palette1.at(iSlice));
            pic.AddFunction(fsim);

            k = rec_slope->At({i1, i2}).Mean();
            b = rec_offset->At({i1, i2}).Mean();
            p3 = pol=="pol3" ? rec_third->At({i1, i2}).Mean() : 0;
            TF1* frec = new TF1("frec", "[0]+[1]*x+[2]*x*x*x", -10, 10);
            frec->SetParameters(b, k, p3);
            frec->SetLineStyle(1);
            frec->SetLineColor(Helper::palette1.at(iSlice));
            pic.AddFunction(frec);
          }
        }

        pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), "v_{1}"});

    //     pic.SetXRange({-0.05, 0.95});
        pic.CustomizeXRange();
        pic.CustomizeYRangeWithLimits(-0.3, 0.6);
        pic.AddLegend(leg1);
        pic.SetIsCustomizeLegend();
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
