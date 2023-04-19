#include "lambda.h"

void lambda_error_sqrtN() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string fileInName = "/home/oleksii/cbmdir/working/qna/rectrackspsi/dcmqgsm/cl.rtp.dcmqgsm.root";

  float y_lo, y_hi;

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");

//   std::string y_axis_title = "v_{1}";
  std::string y_axis_title = "#sigma#upoint#sqrt{2N}";
  std::string particle = "lambda";
//   std::string particle = "kshort";

  axes.at(0).sim_name_ = "SimParticles_pT";
  axes.at(1).sim_name_ = "SimParticles_rapidity";
  axes.at(2).sim_name_ = "RecEventHeader_centrality_tracks";
  axes.at(0).reco_name_ = "ReconstructedParticles_pT";
  axes.at(1).reco_name_ = "ReconstructedParticles_rapidity";
  axes.at(2).reco_name_ = "RecEventHeader_centrality_tracks";

  TFile* fileIn = TFile::Open(fileInName.c_str(), "open");
  if(!fileIn) throw std::runtime_error("fileIn does not exist");
  auto* dc = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>>("sim/u_lambda_sim_PLAIN.Q_psi_PLAIN.x1x1");
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  std::string fileOutName = "fileOut";

  bool is_first_canvas = true;

  SetSliceAxisBinEdges({0, 0.4, 0.8, 1.2, 1.6});

  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

    DoubleDifferentialCorrelation v1cos( fileInName, {("rec/u_" + particle + "_rec_wei_RESCALED.Q_psi_PLAIN.x1x1").c_str(),
                                                      ("rec/u_" + particle + "_rec_wei_RESCALED.Q_psi_PLAIN.y1y1").c_str()} );
//     DoubleDifferentialCorrelation v1cos( fileInName, {("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.x1x1").c_str(),
//                                                       ("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.y1y1").c_str()} );
    v1cos.SetMarker(kFullSquare);
    v1cos.Scale(2.*std::sqrt(2.));
    v1cos.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1cos.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                    kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
    v1cos.SetBiasPalette(false);
    v1cos.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                 {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1cos.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1cos.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1cos.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1cos.ShiftProjectionAxis(axes.at(kProjection).shift_);
    v1cos.SlightShiftProjectionAxis(0.025);
    v1cos.SetDrawErrorAsMean(true, true);
    v1cos.Calculate();

    DoubleDifferentialCorrelation v1sin( fileInName, {("rec/u_" + particle + "_rec_wei_RESCALED.Q_psi_PLAIN.x1y1").c_str(),
                                                      ("rec/u_" + particle + "_rec_wei_RESCALED.Q_psi_PLAIN.y1x1").c_str()} );
//     DoubleDifferentialCorrelation v1sin( fileInName, {("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.x1y1").c_str(),
//                                                       ("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.y1x1").c_str()} );
    v1sin.SetMarker(kOpenSquare);
    v1sin.Scale(2.*std::sqrt(2.));
    v1sin.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1sin.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                    kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
    v1sin.SetBiasPalette(false);
    v1sin.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                 {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1sin.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1sin.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1sin.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1sin.ShiftProjectionAxis(axes.at(kProjection).shift_);
    v1sin.SlightShiftProjectionAxis(0.025, 0.025);
    v1sin.SetDrawErrorAsMean(true, true);
    v1sin.Calculate();

    HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});
    pic.AddText({0.2, 0.90, particle.c_str()}, 0.025);
    pic.AddText({0.2, 0.87, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                            " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);
    pic.DrawZeroLine(false);

    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);
    leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());
    TLegendEntry* entry;
    entry = leg1->AddEntry("", "cos", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerStyle(kFullSquare);
    entry = leg1->AddEntry("", "sin", "P");
    entry->SetMarkerSize(2);
    entry->SetMarkerStyle(kOpenSquare);

    for( auto obj : v1cos.GetProjections() ) {
      pic.AddDrawable( obj );
      leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
    }
    for( auto obj : v1sin.GetProjections() ) {
      pic.AddDrawable( obj );
    }

    pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title});
    pic.AddHorizontalLine(1);
    pic.CustomizeXRange();
    pic.CustomizeYRange();
    pic.AddLegend(leg1);
    pic.CustomizeLegend(leg1);
    pic.Draw();

    if(is_first_canvas) pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
    else pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
    is_first_canvas = false;
  }

  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");
}
