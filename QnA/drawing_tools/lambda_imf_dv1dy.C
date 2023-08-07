#include "lambda.h"

void lambda_imf_dv1dy() {
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

  bool is_imf = true;
  if(pdg=="3312" || (pdg=="3122" && pbeam=="3.3")) is_imf = false;

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileMcName = "/home/oleksii/cbmdir/working/qna/aXmass/vR.dv1dy." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";
  std::string fileRecName = "/home/oleksii/cbmdir/working/qna/aXmass/of.dv1dy." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + ".root";

  DrawOption drawOption = kPlain;
//   DrawOption drawOption = kDifference;
//   DrawOption drawOption = kChi2;
//   DrawOption drawOption = kRatio;

  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");

  if(!is_imf) {
    fileRecName = fileMcName;
    axes.at(0).reco_fit_name_ = axes.at(0).reco_name_;
    axes.at(1).reco_fit_name_ = axes.at(1).reco_name_;
    axes.at(2).reco_fit_name_ = axes.at(2).reco_name_;
  }

  std::string y_axis_title;

//   std::vector<std::string> components{"x1x1", "y1y1"};
  std::vector<std::string> components{"ave"};

  std::vector<std::string> fitcoeffs{"slope", "intercept"};

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
    std::string fileOutName;
    if(drawOption == kPlain) fileOutName = "";
    if(drawOption == kChi2) fileOutName = "chi2.";
    if(drawOption == kDifference) fileOutName = "diff.";
    if(drawOption == kRatio) fileOutName = "ratio.";
    fileOutName += fc + ".imf." + evegen + "." + pbeam + "agev." + pdg;
    if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

    if(fc == "slope") y_axis_title = "dv_{1}/dy";
    if(fc == "intercept") y_axis_title = "v_{1}|_{y=0}";
    if(drawOption == kChi2) y_axis_title = "#frac{est. - ref.}{#sigma#{}{est. - ref.}} of " + y_axis_title;
    if(drawOption == kRatio) y_axis_title = "est. / ref. of " + y_axis_title;
    if(drawOption == kDifference) y_axis_title = "est. - ref. of " + y_axis_title;

    for(auto& ip : inputs) {
      for(auto& se : ip.subevents_) {
        for(auto& co : components) {

          std::string rec_nofit_name = "v1/" + ip.dirname_ + "/" + fc + "/v1.u_rec_sgnl." + se + ip.resname_ + "." + co;
          auto v1_rec_nofit = DoubleDifferentialCorrelation( fileMcName.c_str(), {rec_nofit_name.c_str()} );
          v1_rec_nofit.SetErrorType(error_mode);
          v1_rec_nofit.SetMeanType(mean_mode);
          v1_rec_nofit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
          v1_rec_nofit.SetMarker(-1);
          v1_rec_nofit.SetIsFillLine();
          v1_rec_nofit.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                                   kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
          v1_rec_nofit.SetBiasPalette(false);
          v1_rec_nofit.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
          v1_rec_nofit.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
          v1_rec_nofit.ShiftSliceAxis(axes.at(kSlice).shift_);
          v1_rec_nofit.ShiftProjectionAxis(axes.at(kProjection).shift_);
          v1_rec_nofit.Calculate();

          std::string rec_fit_name;
          if(is_imf) rec_fit_name = "Pars/" + ip.dirname_ + "/" + se + ip.resname_ + "/" + fc  + "/signal." + co;
          else       rec_fit_name = "v1/" + ip.dirname_ + "/" + fc + "/v1.u_rec." + se + ip.resname_ + "." + co;
          auto v1_rec_fit = DoubleDifferentialCorrelation( fileRecName.c_str(), {rec_fit_name.c_str()} );
          v1_rec_fit.RenameAxis(axes.at(kProjection).reco_fit_name_, axes.at(kProjection).reco_name_);
          v1_rec_fit.RenameAxis(axes.at(kSlice).reco_fit_name_, axes.at(kSlice).reco_name_);
          v1_rec_fit.SetErrorType(error_mode);
          v1_rec_fit.SetMeanType(mean_mode);
          v1_rec_fit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
          v1_rec_fit.SetMarker(kFullSquare);
          v1_rec_fit.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                                 kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
          v1_rec_fit.SetBiasPalette(false);
          v1_rec_fit.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
          v1_rec_fit.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
          v1_rec_fit.ShiftSliceAxis(axes.at(kSlice).shift_);
          v1_rec_fit.ShiftProjectionAxis(axes.at(kProjection).shift_);
          v1_rec_fit.SlightShiftProjectionAxis(1);
          v1_rec_fit.Calculate();

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
          if(ip.resname_ != "") pic.AddText({0.2, 0.75, ip.resname_.substr(1, ip.resname_.size()).c_str()}, 0.025);

          auto leg1 = new TLegend();
          leg1->SetBorderSize(1);

          std::string ref_title = "REC, MC-match";
          std::string est_title;
          if(is_imf) est_title = "REC, inv. mass fit";
          else       est_title = "REC, all candidates";

          TLegendEntry* entry;
          if(drawOption == kPlain) {
            entry = leg1->AddEntry("", ref_title.c_str(), "F");
            entry->SetFillColorAlpha(kBlack, 0.2);
            entry->SetLineColor(kWhite);
            entry->SetFillStyle(1000);

            entry = leg1->AddEntry("", est_title.c_str(), "P");
            entry->SetMarkerSize(2);
            entry->SetMarkerStyle(kFullSquare);
          } else {
            leg1->AddEntry("", ("ref: " + ref_title).c_str(), "");
            leg1->AddEntry("", ("est: " + est_title).c_str(), "");
          }

          leg1->AddEntry("", (axes.at(kSlice).title_+axes.at(kSlice).unit_ + ":").c_str(), "");

          if(drawOption != kPlain && drawOption != kDifference) {
            pic.DrawZeroLine(false);
            pic.AddHorizontalLine(1);
          }
          if(drawOption == kChi2) pic.AddHorizontalLine(-1);

          DoubleDifferentialCorrelation vplot;
          if(drawOption == kPlain) vplot = v1_rec_fit;
          if(drawOption == kDifference || drawOption == kChi2) vplot = Minus(v1_rec_fit, v1_rec_nofit);
          if(drawOption == kChi2) vplot.DivideValueByError();
          if(drawOption == kRatio) vplot = Divide(v1_rec_fit, v1_rec_nofit);

          if(drawOption == kPlain) {
            for( auto obj : v1_rec_nofit.GetProjections() ){
              pic.AddDrawable( obj );
            }
          }

          for( auto obj : vplot.GetProjections() ){
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
    }// inputs
    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

    if(is_write_rootfile) fileOut->Close();
  }// fitcoeffs
}
