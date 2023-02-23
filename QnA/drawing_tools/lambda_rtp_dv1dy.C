#include "lambda.h"

void lambda_rtp_dv1dy() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  enum DrawOption {
    kPlain,
    kDifference,
    kChi2,
    kRatio
  };
//   DrawOption drawOption = kPlain;
//   DrawOption drawOption = kDifference;
  DrawOption drawOption = kChi2;
//   DrawOption drawOption = kRatio;

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  std::string step = "PLAIN";
//   std::string step = "RESCALED";

//   std::string weightstatus = "wei";
  std::string weightstatus = "now";

  std::string rebinnedstatus = ".rebinned";
//   std::string rebinnedstatus = "";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/rectrackspsi/" + evegen + "/dv1dy" + rebinnedstatus + ".rtp." + evegen + "." + weightstatus + ".root";

  std::vector<std::string> particles{"lambda", "kshort"};

  bool average_comp{false};
  float y_lo, y_hi;

  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title;

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::vector<std::string> fitcoeffs{"slope", "intercept"};
  if(drawOption == kRatio) fitcoeffs.pop_back();
  if(drawOption == kDifference) fitcoeffs.erase(fitcoeffs.begin());

  axes.at(kSlice).sim_name_ = "SimParticles_pT";
  axes.at(kProjection).sim_name_ = "RecEventHeader_centrality_tracks";
  axes.at(kSlice).reco_name_ = "ReconstructedParticles_pT";
  axes.at(kProjection).reco_name_ = "RecEventHeader_centrality_tracks";

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  if(!fileIn) throw std::runtime_error("fileIn is absent");
  TFile* fileOut{nullptr};

  for(auto& particle : particles) {

    std::string fileOutName;
    if(drawOption == kPlain) fileOutName = "dv1dy." + particle;
    if(drawOption == kChi2) fileOutName = "chi2.dv1dy." + particle;
    if(drawOption == kDifference) fileOutName = "diff.dv1dy." + particle;
    if(drawOption == kRatio) fileOutName = "ratio.dv1dy." + particle;

    fileOutName += rebinnedstatus + "." + step + "." + weightstatus;

    std::ofstream fileOutText;
    fileOutText.open((fileOutName + ".txt").c_str());
    fileOutText << evegen << "\n";

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

    for(auto fc : fitcoeffs) {

      if(evegen == "dcmqgsm") {
        if(particle == "lambda") {
          if(fc == "slope")     { y_lo = -0.2; y_hi = 0.4; }
          if(fc == "intercept") { y_lo = -0.02; y_hi = 0.02; }
        }
        if(particle == "kshort") {
          if(fc == "slope")     { y_lo = -0.1; y_hi = 0.1; }
          if(fc == "intercept") { y_lo = -0.02; y_hi = 0.02; }
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
      }

      if(fc == "slope") y_axis_title = "dv_{1}/dy";
      if(fc == "intercept") y_axis_title = "v_{1}|_{y=0}";
      if(drawOption == kChi2) y_axis_title = "#chi^{2} of " + y_axis_title;
      if(drawOption == kRatio) y_axis_title = "REC / MC of " + y_axis_title;

      std::vector<DoubleDifferentialCorrelation> v1_rec;
      if(average_comp) {
        v1_rec.resize(1);
        v1_rec.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + step + ".ave").c_str()} );
        v1_rec.at(0).SetMarker(kFullSquare);
      } else {
        v1_rec.resize(2);
        v1_rec.at(0) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + step + ".x1x1").c_str()} );
        v1_rec.at(1) = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1rec_" + fc + "." + step + ".y1y1").c_str()} );
        v1_rec.at(0).SetMarker(kFullSquare);
        v1_rec.at(1).SetMarker(kOpenSquare);
      }

      for(auto& vc : v1_rec) {
        vc.Scale(2.);
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

      v1_rec.at(0).SlightShiftProjectionAxis(1);
      if(!average_comp) {
        v1_rec.at(1).SlightShiftProjectionAxis(1, 0.5);
      }

      auto v1_mc = DoubleDifferentialCorrelation( fileName.c_str(), {(particle + "/v1sim_" + fc + ".psi.ave").c_str()} );
      v1_mc.Scale(2.);
      v1_mc.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_mc.SetMarker(-1);
      v1_mc.SetIsFillLine();
      v1_mc.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                            kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
      v1_mc.SetBiasPalette(false);
      v1_mc.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_mc.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_mc.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_mc.Calculate();
      v1_mc.ShiftProjectionAxis(axes.at(kProjection).shift_);

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
      pic.AddText({0.2, 0.75, weightstatus.c_str()}, 0.025);

      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);
      leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

      TLegendEntry* entry;
      entry = leg1->AddEntry("", "MC, x&y averaged", "L");
      entry->SetMarkerSize(2);
      if(average_comp) {
        entry = leg1->AddEntry("", "Reco, x&y averaged", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kFullSquare);
      } else {
        entry = leg1->AddEntry("", "Reco, X", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kFullSquare);
        entry = leg1->AddEntry("", "Reco, Y", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kOpenSquare);
      }

      if(drawOption == kPlain) {
        for(int i=0; i<v1_rec.size(); i++) {
          for( auto obj : v1_rec.at(i).GetProjections() ){
            pic.AddDrawable( obj );
            if(i==0) leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
          }
        }
        for( auto obj : v1_mc.GetProjections() ){
          pic.AddDrawable( obj );
        }
      } else {
        if(drawOption != kPlain && drawOption != kDifference) {
          pic.DrawZeroLine(false);
          pic.AddHorizontalLine(1);
        }
        if(drawOption == kChi2) pic.AddHorizontalLine(-1);
        for(int i=0; i<v1_rec.size(); i++) { // i: X, Y, ave
          int j{0}; // j: slice axis
          for( auto obj_R_MC : v1_rec.at(i).GetProjections() ) {
            std::string comp;
            if(v1_rec.size() == 2) {
              if(i==0) comp = "X";
              if(i==1) comp = "Y";
            } else {
              comp = "AVE";
            }
            auto obj_psiRP = v1_mc.GetProjections().at(j);
            Graph* obj_R_MC_psiRP{nullptr};
            std::vector<float> vec_values{};
            std::vector<float> vec_errors{};
            if(drawOption == kChi2 || drawOption == kDifference) {
              GraphSubtractor gr_sub;
              gr_sub.SetMinuend(obj_R_MC);
              gr_sub.SetSubtrahend(obj_psiRP);
              if(drawOption == kDifference) gr_sub.DivideOverError(false);
              gr_sub.Calculate();
              obj_R_MC_psiRP = gr_sub.GetResult();
              vec_values = gr_sub.GetPointsValues();
              vec_errors = gr_sub.GetPointsErrors();
            }
            if(drawOption == kRatio) {
              GraphDivider gr_div;
              gr_div.SetNumerator(obj_R_MC);
              gr_div.SetDenominator(obj_psiRP);
              gr_div.Calculate();
              obj_R_MC_psiRP = gr_div.GetResult();
              vec_values = gr_div.GetPointsValues();
              vec_errors = gr_div.GetPointsErrors();
            }
            pic.AddDrawable(obj_R_MC_psiRP);
            fileOutText << fc << "\t" << "pt" << j+1 << "\t" << comp << "\t";
            for(int iv=0; iv<vec_values.size(); iv++) {
              fileOutText << vec_values.at(iv) << "\t" << vec_errors.at(iv) << "\t";
            }
            fileOutText << "\n";
            if(i==0) leg1->AddEntry( obj_R_MC->GetPoints(), obj_R_MC->GetTitle().c_str(), "P" );
            j++;
          } // j
        } // i
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
    }

    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

    if(is_write_rootfile) fileOut->Close();

    fileOutText.close();
  }
}
