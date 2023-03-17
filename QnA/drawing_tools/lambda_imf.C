#include "lambda.h"

void lambda_imf() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

//   std::string particle = "lambda";
  std::string particle = "kshort";

  std::string fileNameMc = "/home/oleksii/cbmdir/working/qna/inv_mass_flow/cl.imf." + evegen + ".12agev.root";
  std::string fileNameImf = "/home/oleksii/cbmdir/working/qna/inv_mass_flow/out.fitter." + evegen + ".12agev." + particle + ".root";

//   bool average_comp = true; std::vector<std::string> components{"ave"};
  bool average_comp = false; std::vector<std::string> components{"x1x1", "y1y1"};

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  float y_lo, y_hi;
  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title = "v_{1}";

  TFile* fileIn = TFile::Open(fileNameMc.c_str(), "open");
  if(!fileIn) throw std::runtime_error("fileIn does not exist");
  auto* dc = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>>("sim/u_lambda_sim_PLAIN.Q_psi_PLAIN.x1x1");
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  TFile* fileOut{nullptr};

  if(evegen == "dcmqgsm") {
    if(particle == "lambda") {
      y_lo = -0.2;
      y_hi = 0.3;
    }
    if(particle == "kshort") {
      y_lo = -0.1;
      y_hi = 0.2;
    }
  }
  if(evegen == "urqmd") {
    if(particle == "lambda") {
      y_lo = -0.2;
      y_hi = 0.2;
    }
    if(particle == "kshort") {
      y_lo = -0.2;
      y_hi = 0.2;
    }
  }

  bool is_first_canvas = true;
  std::string fileOutName;

  fileOutName = "v1.imf." + particle;

  if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

    for(auto& co : components) {
      DoubleDifferentialCorrelation v1_rec_imf;
      if(average_comp) v1_rec_imf = DoubleDifferentialCorrelation(fileNameImf.c_str(), {"Params/signal.x1x1",
                                                                                        "Params/signal.y1y1"});
      else             v1_rec_imf = DoubleDifferentialCorrelation(fileNameImf.c_str(), {("Params/signal." + co).c_str()});
      v1_rec_imf.SetMarker(kFullSquare);
      v1_rec_imf.Scale(2.);
      v1_rec_imf.SetErrorType(error_mode);
      v1_rec_imf.SetMeanType(mean_mode);
      v1_rec_imf.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_rec_imf.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                             kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
      v1_rec_imf.Rebin({{axes.at(kSelect).reco_fit_name_.c_str(),
                        {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
      v1_rec_imf.SetProjectionAxis({axes.at(kProjection).reco_fit_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_rec_imf.SetSliceAxis({axes.at(kSlice).reco_fit_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_rec_imf.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_rec_imf.ShiftProjectionAxis(axes.at(kProjection).shift_);
      v1_rec_imf.SlightShiftProjectionAxis(0.02, 0.01);
      v1_rec_imf.Calculate();

      DoubleDifferentialCorrelation v1_rec_mc;
      if(average_comp) v1_rec_mc = DoubleDifferentialCorrelation(fileNameMc.c_str(), {"rec/u_" + particle + "_rec_mc_RESCALED.Q_psi_PLAIN.x1x1",
                                                                                      "rec/u_" + particle + "_rec_mc_RESCALED.Q_psi_PLAIN.y1y1"});
      else             v1_rec_mc = DoubleDifferentialCorrelation(fileNameMc.c_str(), {("rec/u_" + particle + "_rec_mc_RESCALED.Q_psi_PLAIN." + co).c_str()});
      v1_rec_mc.SetMarker(kOpenSquare);
      v1_rec_mc.Scale(2.);
      v1_rec_mc.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_rec_mc.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                            kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
      v1_rec_mc.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                       {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
      v1_rec_mc.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_rec_mc.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_rec_mc.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_rec_mc.ShiftProjectionAxis(axes.at(kProjection).shift_);
      v1_rec_mc.SlightShiftProjectionAxis(0.02);
      v1_rec_mc.Calculate();

      auto v1_ref = DoubleDifferentialCorrelation( fileNameMc.c_str(), {("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.x1x1").c_str(),
                                                                        ("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.y1y1").c_str()});
      v1_ref.Scale(2.);
      v1_ref.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_ref.SetMarker(-1);
      v1_ref.SetIsFillLine();
      v1_ref.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                         kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
      v1_ref.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                    {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
      v1_ref.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_ref.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_ref.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_ref.ShiftProjectionAxis(axes.at(kProjection).shift_);
      v1_ref.Calculate();

      HeapPicture pic(particle, {1000, 1000});
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
      pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                               " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
      pic.AddText({0.2, 0.75, co.c_str()}, 0.025);

      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);
      leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

      TLegendEntry* entry;
      entry = leg1->AddEntry("", "REF: sim-tracks", "L");
      entry->SetMarkerSize(2);
      entry = leg1->AddEntry("", "EST: reco-tracks, MC-match", "P");
      entry->SetMarkerSize(2);
      entry->SetMarkerStyle(kOpenSquare);
      entry = leg1->AddEntry("", "EST: reco-tracks, invmassfit", "P");
      entry->SetMarkerSize(2);
      entry->SetMarkerStyle(kFullSquare);

      for( auto obj : v1_rec_imf.GetProjections() ){
        pic.AddDrawable( obj );
        leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
      }
      for( auto obj : v1_rec_mc.GetProjections() ){
        pic.AddDrawable( obj );
      }
      for( auto obj : v1_ref.GetProjections() ){
        pic.AddDrawable( obj );
      }

      pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title});
      pic.CustomizeXRange();
      pic.CustomizeYRangeWithLimits(y_lo, y_hi);
      pic.AddLegend(leg1);
      pic.CustomizeLegend(leg1);
      pic.Draw();

      if(is_write_rootfile) {
        fileOut->cd();
        pic.GetCanvas()->Write();
      }

      if(is_first_canvas) pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
      else pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
      is_first_canvas = false;
    }
  }
  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

  if(is_write_rootfile) fileOut->Close();
}
