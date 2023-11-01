#include "lambda.h"

void lambda_complex_dv1dy_syst() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );
//   gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.01);
//   gStyle->SetEndErrorSize(5);

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3"; axes.at(1).shift_ = -0.985344;
  std::string evegen = "urqmd";   std::string pbeam = "12";

//   std::string particle = "#Lambda"; std::string pdg = "3122";
  std::string particle = "K^{0}_{S}"; std::string pdg = "310";
// //   std::string particle = "#Xi^{-}"; std::string pdg = "3312";

  //   std::string is_fine_pt = "";
  std::string is_fine_pt = "_finept";

  std::string cuts = "lc1";
  if(pbeam == "3.3") cuts = "oc1";
  if(pdg == "3312") cuts = "dc";

  bool is_imf = true;
  if(pdg=="3312" || (pdg=="3122" && pbeam=="3.3")) is_imf = false;

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileMcName = "/home/oleksii/cbmdir/working/qna/aXmass/vR.dv1dy." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + is_fine_pt + ".root";
  std::string fileRecName = "/home/oleksii/cbmdir/working/qna/aXmass/of.dv1dy." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + is_fine_pt + ".root";

//   SetAxis("centrality", "projection");
//   SetAxis("pT", "slice");
  SetAxis("centrality", "slice");
  SetAxis("pT", "projection");

  if(!is_imf) {
    fileRecName = fileMcName;
    axes.at(0).reco_fit_name_ = axes.at(0).reco_name_;
    axes.at(1).reco_fit_name_ = axes.at(1).reco_name_;
    axes.at(2).reco_fit_name_ = axes.at(2).reco_name_;
  }

  std::string y_axis_title;

  std::vector<std::string> components{"x1x1", "y1y1"};
//   std::vector<std::string> components{"ave"};

  std::vector<std::string> fitcoeffs{"slope", "intercept"};

  struct Inputs {
    std::string dirname_;
    std::vector<std::string> subevents_;
    std::string resname_;
  };

  std::vector<Inputs> inputs {
    {"uPsi", {"Q_psi"}, ""},
    {"uQ_R1_MC", {"psd1", "psd2", "psd3"}, "_res_MC"},
    {"uQ_R1_sub3", {"psd1", "psd2", "psd3"}, "_res_sub3"},
    {"uQ_R1_sub4", {"psd1", "psd2", "psd3"}, "_res_sub4_sts_pipos"},
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

  float slightprojshift;
  if(axes.at(kSlice).name_ == "centrality") {
    if(pbeam == "12") SetSliceAxisBinEdges({0, 5, 20, 40, 70});
    else              SetSliceAxisBinEdges({0, 20, 40, 70});
    slightprojshift = 0.015;
  } else {
    slightprojshift = 1;
  }

  TFile* fileOut;

  for(auto& fc : fitcoeffs) {

    bool is_first_canvas{true};
    std::string fileOutName = fc + "." + evegen + "." + pbeam + "agev." + pdg;
    if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

    if(fc == "slope") y_axis_title = "dv_{1}/dy";
    if(fc == "intercept") y_axis_title = "v_{1}|_{y=0}";

    for(auto& ip : inputs) {

      std::vector<std::string> v1_sim_names;
      for(auto& co : components) {
        v1_sim_names.emplace_back("v1/usimPsi/" + fc + "/v1.u_sim.Q_psi." + co);
      }
      auto v1_sim = DoubleDifferentialCorrelation( fileMcName.c_str(), v1_sim_names );
      v1_sim.SetErrorType(error_mode);
      v1_sim.SetMeanType(mean_mode);
      v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_sim.SetMarker(-1);
      v1_sim.SetIsFillLine();
      v1_sim.SetPalette(Helper::palette1);
      v1_sim.SetBiasPalette(false);
      v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);
      v1_sim.SetSliceAxisPrecision(axes.at(kSlice).precision_);
      v1_sim.Calculate();

      std::vector<std::string> v1_rec_nofit_names;
      for(auto& se : ip.subevents_) {
        for(auto& co : components) {
          v1_rec_nofit_names.emplace_back("v1/" + ip.dirname_ + "/" + fc + "/v1.u_rec_sgnl." + se + ip.resname_ + "." + co);
        }
      }
      auto v1_rec_nofit = DoubleDifferentialCorrelation( fileMcName.c_str(), v1_rec_nofit_names );
      v1_rec_nofit.SetErrorType(error_mode);
      v1_rec_nofit.SetMeanType(mean_mode);
      v1_rec_nofit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_rec_nofit.SetMarker(kOpenSquare);
      v1_rec_nofit.SetPalette(Helper::palette1);
      v1_rec_nofit.SetBiasPalette(false);
      v1_rec_nofit.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_rec_nofit.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_rec_nofit.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_rec_nofit.ShiftProjectionAxis(axes.at(kProjection).shift_);
      v1_rec_nofit.SlightShiftProjectionAxis(slightprojshift, slightprojshift/2);
      v1_rec_nofit.SetSliceAxisPrecision(axes.at(kSlice).precision_);
      v1_rec_nofit.Calculate();

      std::vector<std::string> v1_rec_fit_names;
      for(auto& se : ip.subevents_) {
        for(auto& co : components) {
          if(is_imf) v1_rec_fit_names.emplace_back("Pars/" + ip.dirname_ + "/" + se + ip.resname_ + "/" + fc  + "/signal." + co);
          else       v1_rec_fit_names.emplace_back("v1/" + ip.dirname_ + "/" + fc + "/v1.u_rec." + se + ip.resname_ + "." + co);
        }
      }
      auto v1_rec_fit = DoubleDifferentialCorrelation( fileRecName.c_str(), v1_rec_fit_names );
      v1_rec_fit.SetErrorType(error_mode);
      v1_rec_fit.SetMeanType(mean_mode);
      v1_rec_fit.SetCalculateSystematicsFromVariation();
      v1_rec_fit.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_rec_fit.SetMarker(kFullSquare);
      v1_rec_fit.SetPalette(Helper::palette1);
      v1_rec_fit.SetBiasPalette(false);
      v1_rec_fit.SetProjectionAxis({axes.at(kProjection).reco_fit_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_rec_fit.SetSliceAxis({axes.at(kSlice).reco_fit_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_rec_fit.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_rec_fit.ShiftProjectionAxis(axes.at(kProjection).shift_);
      v1_rec_fit.SlightShiftProjectionAxis(slightprojshift);
      v1_rec_fit.SetSliceAxisPrecision(axes.at(kSlice).precision_);
      v1_rec_fit.Calculate();

      HeapPicture pic(fc, {1000, 1000});

      const float text_size = 24;
      const int text_font = 63;
      float text_X = 0.04;
      float text_Y = 0.96;
      if(fc == "slope") {
        pic.AddText(particle.c_str(), {text_X, text_Y}, text_size+8, text_font);
        if(evegen == "dcmqgsm") {
          if(pbeam == "12") pic.AddText("5M Au+Au", {text_X, text_Y - 0.04}, text_size, text_font);
          else              pic.AddText("5.2M Au+Au", {text_X, text_Y - 0.04}, text_size, text_font);
          pic.AddText("DCM-QGSM-SMM", {text_X, text_Y - 0.08}, text_size, text_font);
        }
        if(evegen == "urqmd") {
          pic.AddText("2M Au+Au", {text_X, text_Y - 0.04}, text_size, text_font);
          pic.AddText("UrQMD", {text_X, text_Y - 0.08}, text_size, text_font);
        }
        pic.AddText((pbeam + "A GeV/c").c_str(), {text_X, text_Y - 0.12}, text_size, text_font);
        if(ip.resname_ != "") pic.AddText(ip.resname_.substr(1, ip.resname_.size()).c_str(), {text_X, text_Y - 0.16}, text_size, text_font);
      }

      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);

      auto* entry = leg1->AddEntry("", "MC input", "F");
      entry->SetFillColorAlpha(kBlack, 0.2);
      entry->SetLineColor(kWhite);
      entry->SetFillStyle(1000);

//           entry = leg1->AddEntry("", "REC, MC-match", "P");
//           entry->SetMarkerSize(2);
//           entry->SetMarkerStyle(kOpenSquare);

      entry = leg1->AddEntry("", "Reconstructed", "P");
      entry->SetMarkerSize(2);
      entry->SetMarkerStyle(kFullSquare);

      leg1->AddEntry("", (axes.at(kSlice).title_+axes.at(kSlice).unit_ + ":").c_str(), "");

      for( auto obj : v1_sim.GetProjections() ){
        pic.AddDrawable( obj );
      }
//           for( auto obj : v1_rec_nofit.GetProjections() ){
//             pic.AddDrawable( obj );
//           }
      for( auto obj : v1_rec_fit.GetProjections() ){
        pic.AddDrawable( obj );
        leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
      }

      pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title.c_str()});

  //     pic.SetXRange({-0.05, 0.95});
      pic.CustomizeXRange();
      pic.CustomizeYRange();
//       if(fc == "slope") pic.SetYRange({-0.079, 0.46});
//       else              pic.SetYRange({-0.075, 0.019});
      if(fc == "slope") pic.AddLegend(leg1);
      pic.SetIsCustomizeLegend();
      pic.Draw();

      if(is_write_rootfile) {
        fileOut->cd();
        pic.Write(("heap_picture_" + ip.dirname_).c_str());
        pic.GetCanvas()->Write(("canvas_" + ip.dirname_).c_str());
      }

      if(is_first_canvas)
        pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
      else
        pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
      is_first_canvas = false;

    }// inputs
    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

    if(is_write_rootfile) fileOut->Close();
  }// fitcoeffs
}
