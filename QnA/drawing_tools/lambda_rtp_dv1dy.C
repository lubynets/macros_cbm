#include "lambda.h"

void lambda_rtp_dv1dy() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  enum DrawOption {
    kPlain,
    kDifference,
    kChi2,
    kRatio
  };
  DrawOption drawOption = kPlain;
//   DrawOption drawOption = kDifference;
//   DrawOption drawOption = kChi2;
//   DrawOption drawOption = kRatio;

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

//   bool average_comp{true}; std::vector<std::string> components{"ave"}; std::string avestatus = ".ave";
  bool average_comp{false}; std::vector<std::string> components{"x1x1", "y1y1"}; std::string avestatus = "";

  std::vector<std::string> particles{"lambda", "kshort"};
  std::vector<std::string> steps{"PLAIN", "RECENTERED", "TWIST", "RESCALED"};
  std::vector<std::string> weightstatus{"now", "wei"};

  std::string fileName = "/home/oleksii/cbmdir/working/qna/rectrackspsi/" + evegen + "/dv1dy.rebinned.rtp." + evegen + avestatus + ".root";

  float y_lo, y_hi;

  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title;

  std::vector<std::string> fitcoeffs{"slope", "intercept"};
  if(drawOption == kRatio) fitcoeffs.pop_back();
  if(drawOption == kDifference) fitcoeffs.erase(fitcoeffs.begin());

  axes.at(kSlice).sim_name_ = "SimParticles_pT";
  axes.at(kProjection).sim_name_ = "RecEventHeader_centrality_tracks";
  axes.at(kSlice).reco_name_ = "ReconstructedParticles_pT";
  axes.at(kProjection).reco_name_ = "RecEventHeader_centrality_tracks";

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  TFile* fileOut{nullptr};

  for(auto& particle : particles) {

    std::string fileOutName;
    if(drawOption == kPlain) fileOutName = "dv1dy." + particle;
    if(drawOption == kChi2) fileOutName = "chi2.dv1dy." + particle;
    if(drawOption == kDifference) fileOutName = "diff.dv1dy." + particle;
    if(drawOption == kRatio) fileOutName = "ratio.dv1dy." + particle;

    std::ofstream fileOutText;
    fileOutText.open((fileOutName + ".txt").c_str());
    fileOutText << evegen << "\n";
    fileOutText << "particle\t" << particle << "\n";

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

    bool is_first_canvas = true;
    if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

    for(auto& step : steps) {
      if(step != "PLAIN" && step!= "RESCALED") continue;
      for(auto& ws : weightstatus) {

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

          std::vector<DoubleDifferentialCorrelation> v1_est;
          if(average_comp) {
            v1_est.resize(1);
            v1_est.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + ws + "." + step + ".ave").c_str()} );
            v1_est.at(0).SetMarker(kFullSquare);
            v1_est.at(0).SlightShiftProjectionAxis(1);
          } else {
            v1_est.resize(2);
            v1_est.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + ws + "." + step + ".x1x1").c_str()} );
            v1_est.at(1) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + ws + "." + step + ".y1y1").c_str()} );
            v1_est.at(0).SetMarker(kFullSquare);
            v1_est.at(1).SetMarker(kOpenSquare);
            v1_est.at(1).SlightShiftProjectionAxis(1, 0.5);
          }

          for(auto& vc : v1_est) {
            vc.Scale(2.);
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

          auto v1_ref = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1sim_" + fc + ".psi.ave").c_str()} );
          v1_ref.Scale(2.);
          v1_ref.RenameAxis("SimParticles_pT", "ReconstructedParticles_pT");
          v1_ref.SetErrorType(error_mode);
          v1_ref.SetMeanType(mean_mode);
          v1_ref.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
          v1_ref.SetMarker(-1);
          v1_ref.SetIsFillLine();
          v1_ref.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                             kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
          v1_ref.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
          v1_ref.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
          v1_ref.ShiftSliceAxis(axes.at(kSlice).shift_);
          v1_ref.ShiftProjectionAxis(axes.at(kProjection).shift_);
          v1_ref.Calculate();

          HeapPicture pic(fc, {1000, 1000});
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
          pic.AddText({0.2, 0.78, step.c_str()}, 0.025);
          pic.AddText({0.2, 0.75, ws.c_str()}, 0.025);

          auto leg1 = new TLegend();
          leg1->SetBorderSize(1);
          leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

          TLegendEntry* entry;
          entry = leg1->AddEntry("", "REF: sim-tracks", "L");
          entry->SetMarkerSize(2);
          entry = leg1->AddEntry("", "EST: reco-tracks", "P");
          entry->SetMarkerSize(2);
          entry->SetMarkerStyle(kFullSquare);

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
              std::string comp;
              if(vplot.size() == 2) {
                if(iv==0) comp = "X";
                if(iv==1) comp = "Y";
              } else {
                comp = "AVE";
              }
              std::vector<float> vec_values = obj->GetPointsValues();
              std::vector<float> vec_errors = obj->GetPointsErrors();
              pic.AddDrawable(obj);
              if(iv==0) leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
              fileOutText << step << "\t" << ws << "\t" << fc << "\t" << "pt" << j+1 << "\t" << comp << "\t";
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
          pic.CustomizeLegend(leg1);
          pic.Draw();

          if(is_write_rootfile) {
            fileOut->cd();
            pic.GetCanvas()->Write();
  //           pic.GetCanvas()->SaveAs("cc.C");
          }

          if(is_first_canvas)
            pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
          else
            pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
          is_first_canvas = false;
        } // fitcoeffs
      } // weightstatus
    } // steps
    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");
    if(is_write_rootfile) fileOut->Close();
    fileOutText.close();
  } // particles
}
