#include "lambda.h"

void lambda_stf_dv1dy_mod() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm"; std::string pbeam = "12";
//   std::string evegen = "dcmqgsm"; std::string pbeam = "3.3";
//   std::string evegen = "urqmd";   std::string pbeam = "12";

  DrawOption drawOption = kPlain;
//   DrawOption drawOption = kDifference;
//   DrawOption drawOption = kChi2;
//   DrawOption drawOption = kRatio;

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
//   Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::BOOTSTRAP};

//   Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  enum Mode {
    kModel,
    kDetector
  };
  Mode mode = kModel;
//   Mode mode = kDetector;

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen  + "/" + pbeam + "agev/v1andR1.dv1dy.stf." + evegen + "." + pbeam + "agev.root";
//   std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen  + "/" + pbeam + "agev/dv1dy.rebinned.stf." + evegen + "." + pbeam + "agev.root";

  std::vector<std::string> particles{
                                     "lambda",
//                                      "kshort",
//                                      "xi",
//                                      "pipos",
//                                      "pineg"
                                    };
  std::vector<std::string> subevents{
                                     /*"etacut_1_charged", */"etacut_2_charged"/*, "etacut_3_charged",*/
//                                      "etacut_1_all", "etacut_2_all", "etacut_3_all",
//                                      "psd1", "psd2", "psd3"
                                    };

  std::vector<std::string> subevents_mod{
                                     "etacut_1_charged", "etacut_2_charged", "etacut_3_charged",
//                                      "etacut_1_all", "etacut_2_all", "etacut_3_all"
                                    };
  std::vector<std::string> subevents_det{"psd1", "psd2", "psd3"};


  std::string step;
  bool average_comp;
  float y_lo, y_hi;

  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title;

  std::vector<std::string> fitcoeffs{"slope", "intercept"};
  if(drawOption == kRatio) fitcoeffs.pop_back();
  if(drawOption == kDifference) fitcoeffs.erase(fitcoeffs.begin());

  axes.at(kSlice).sim_name_ = "SimParticles_pT";
  axes.at(kProjection).sim_name_ = "SimEventHeader_centrality_impactpar";
  axes.at(kSlice).reco_name_ = "SimParticles_pT";
  axes.at(kProjection).reco_name_ = "SimEventHeader_centrality_impactpar";

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  if(fileIn==nullptr) throw std::runtime_error("fileIn is nullptr");
  TFile* fileOut{nullptr};

  for(auto& particle : particles) {

    std::string greek_particle = particle;
    if(particle == "lambda") greek_particle = "#Lambda";

    std::string fileOutName;
    if(drawOption == kPlain) fileOutName = "dv1dy." + particle;
    if(drawOption == kChi2) fileOutName = "chi2.dv1dy." + particle;
    if(drawOption == kDifference) fileOutName = "diff.dv1dy." + particle;
    if(drawOption == kRatio) fileOutName = "ratio.dv1dy." + particle;

    std::ofstream fileOutText;
    fileOutText.open((fileOutName + ".txt").c_str());
    fileOutText << evegen << "\n";

    std::cout << particle << "\n";
    auto* dc = (Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>>(("v1/" + particle + "/uPsi/slope/v1.uPsi.ave").c_str());
    if(dc==nullptr) throw std::runtime_error("DataContainer is nullptr");
    for(auto& ax : axes) {
      if(ax.name_ == "rapidity") continue;
      Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
      ax.bin_edges_.clear();
      for(int i=0; i<=qnaxis.size(); i++) {
        ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
      }
    }
    delete dc;

    std::vector<std::string> subevents_loop;
    if(mode == kModel)    subevents_loop = subevents;
    if(mode == kDetector) subevents_loop = subevents_mod;

    int ise{0};
    for(auto& subevent : subevents_loop) {
      if(subevent[0] == 'e') {
        step = "PLAIN";
        average_comp = true;
      }
      if(subevent[0] == 'p') {
        step = "RECENTERED";
        average_comp = true;
      }

      bool is_first_canvas = true;

      std::string fileOutNameSubE;
      if(mode == kModel)    fileOutNameSubE = fileOutName + ".sub2psi." + subevent;
      if(mode == kDetector) fileOutNameSubE = fileOutName + ".sub2sub." + subevent;
      if(is_write_rootfile) fileOut = TFile::Open((fileOutNameSubE + ".root").c_str(), "recreate");

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
        if(drawOption == kChi2) y_axis_title = "#chi^{2} of " + y_axis_title;
        if(drawOption == kRatio) y_axis_title = "REC / MC of " + y_axis_title;

        DoubleDifferentialCorrelation v1_ref;
        if(mode == kModel)    v1_ref = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uPsi/" + fc + "/v1.uPsi.ave").c_str()} );
        if(mode == kDetector) v1_ref = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevents_mod.at(ise) + "_PLAIN.ave").c_str()} );
        v1_ref.SetErrorType(error_mode);
        v1_ref.SetMeanType(mean_mode);
        v1_ref.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        v1_ref.SetMarker(-1);
        v1_ref.SetIsFillLine();
        v1_ref.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                           kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
        v1_ref.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
        v1_ref.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
        v1_ref.ShiftSliceAxis(axes.at(kSlice).shift_);
        v1_ref.ShiftProjectionAxis(axes.at(kProjection).shift_);
        v1_ref.Calculate();

        std::vector<DoubleDifferentialCorrelation> v1_est;
        if(average_comp) {
          v1_est.resize(1);
          if(mode == kModel) v1_est.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevent + "_" + step + ".ave").c_str()} );
          if(mode == kDetector) v1_est.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevents_det.at(ise) + "_RECENTERED.ave").c_str()} );
          v1_est.at(0).SetMarker(kFullSquare);
          v1_est.at(0).SlightShiftProjectionAxis(1);
        } else {
          v1_est.resize(2);
          if(mode == kModel) {
            v1_est.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevent + "_" + step + ".x1x1").c_str()} );
            v1_est.at(1) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevent + "_" + step + ".y1y1").c_str()} );
          }
          if(mode == kDetector) {
            v1_est.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevents_det.at(ise) + "_RECENTERED.x1x1").c_str()} );
            v1_est.at(1) = DoubleDifferentialCorrelation( fileName.c_str(), {("v1/" + particle + "/uQ_R1/" + fc + "/v1.uQ_R1." + subevents_det.at(ise) + "_RECENTERED.y1y1").c_str()} );
          }
          v1_est.at(0).SetMarker(kFullSquare);
          v1_est.at(1).SetMarker(kOpenSquare);
          v1_est.at(1).SlightShiftProjectionAxis(1, 0.5);
        }

        for(auto& vc : v1_est) {
          vc.SetErrorType(error_mode);
          vc.SetMeanType(mean_mode);
          vc.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
          vc.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                          kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
          vc.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
          vc.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
          vc.ShiftSliceAxis(axes.at(kSlice).shift_);
          vc.ShiftProjectionAxis(axes.at(kProjection).shift_);
          vc.Calculate();
        }

        HeapPicture pic(fc, {1000, 1000});
        pic.AddText({0.2, 0.90, greek_particle.c_str()}, 0.025);
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
//         leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

        TLegendEntry* entry;
        if(mode == kModel)    entry = leg1->AddEntry("", "REF: #Psi^{RP}, x&y averaged", "F");
        if(mode == kDetector) entry = leg1->AddEntry("", "REF: eta-cut, R_{MC}, x&y averaged", "F");
        entry->SetMarkerSize(2);
        entry->SetFillColorAlpha(kBlack, 0.2);
        entry->SetLineColor(kWhite);
        entry->SetFillStyle(1000);
        if(average_comp) {
          if(mode == kModel)    entry = leg1->AddEntry("", "EST: #Psi^{rec}, R_{MC}, x&y averaged", "P");
          if(mode == kDetector) entry = leg1->AddEntry("", "EST: PSD, R_{MC}, x&y averaged", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kFullSquare);
        } else {
          if(mode == kModel)    entry = leg1->AddEntry("", "EST: #Psi^{rec}, R_{MC}, X", "P");
          if(mode == kDetector) entry = leg1->AddEntry("", "EST: PSD, R_{MC}, x&y averaged", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kFullSquare);
          if(mode == kModel)    entry = leg1->AddEntry("", "EST: #Psi^{rec}, R_{MC}, Y", "P");
          if(mode == kDetector) entry = leg1->AddEntry("", "EST: PSD, R_{MC}, x&y averaged", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kOpenSquare);
        }
        leg1->AddEntry("", (axes.at(kSlice).title_+axes.at(kSlice).unit_ + ":").c_str(), "");

        if(drawOption != kPlain && drawOption != kDifference) {
          pic.DrawZeroLine(false);
          pic.AddHorizontalLine(1);
        }
        if(drawOption == kChi2) pic.AddHorizontalLine(-1);

        std::vector<DoubleDifferentialCorrelation> vplot;
        vplot.resize(v1_est.size());
        for(int iv=0; iv<v1_est.size(); iv++) { // iv: X, Y, ave
          if(drawOption == kPlain) vplot.at(iv) = v1_est.at(iv);
          if(drawOption == kDifference || drawOption == kChi2) vplot.at(iv) = Minus(v1_est.at(iv), v1_ref);
          if(drawOption == kChi2) vplot.at(iv).DivideValueByError();
          if(drawOption == kRatio) vplot.at(iv) = Divide(v1_est.at(iv), v1_ref);

          int j{0}; // j: slice axis
          for( auto obj : vplot.at(iv).GetProjections() ) {
            std::vector<float> vec_values = obj->GetPointsValues();
            std::vector<float> vec_errors = obj->GetPointsErrors();
            pic.AddDrawable(obj);
            if(iv==0) leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
            fileOutText << subevent << "\t" << fc << "\t" << "pt" << j+1 << "\t" << "AVE" << "\t";
            for(int iv=0; iv<vec_values.size(); iv++) {
              fileOutText << vec_values.at(iv) << "\t" << vec_errors.at(iv) << "\t";
            }
            fileOutText << "\n";

            j++;
          } // j
        } // iv
        if(drawOption == kPlain) {
          for( auto obj : v1_ref.GetProjections() ){
            pic.AddDrawable( obj );
          }
        }
        pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title});

        pic.CustomizeXRange();
        if(drawOption != kPlain) {
          pic.CustomizeYRange();
        } else {
          pic.SetYRange({y_lo, y_hi});
        }
        pic.AddLegend(leg1);
        pic.SetIsCustomizeLegend();
        pic.Draw();

        if(is_write_rootfile) {
          fileOut->cd();
          pic.GetCanvas()->Write();
//           pic.GetCanvas()->SaveAs("cc.C");
        }

        if(is_first_canvas)
          pic.GetCanvas()->Print((fileOutNameSubE + ".pdf(").c_str(), "pdf");
        else
          pic.GetCanvas()->Print((fileOutNameSubE + ".pdf").c_str(), "pdf");
        is_first_canvas = false;
      }

      TCanvas emptycanvas("", "", 1000, 1000);
      emptycanvas.Print((fileOutNameSubE + ".pdf]").c_str(), "pdf");

      if(is_write_rootfile) fileOut->Close();
      ise++;
    }
    fileOutText.close();
  }
}
