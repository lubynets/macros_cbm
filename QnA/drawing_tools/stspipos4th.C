#include "lambda.h"

void stspipos4th() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );
  gStyle->SetNdivisions(206, "xy");

  std::string fileInName = "/home/oleksii/cbmdir/working/qna/aXmass/cl.stspipos4th.root";

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");

  TFile* fileIn = TFile::Open(fileInName.c_str(), "open");
  auto* dc = (QnDcSCo*)fileIn->Get<QnDcSCo>("sim/u_pipos_sim_PLAIN.Q_psi_PLAIN.x1x1");
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  std::string fileOutName = "stspipos4th";

  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  SetSliceAxisBinEdges({0, 0.4, 0.8, 1.4});

  bool is_first_canvas = true;
  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

    DoubleDifferentialCorrelation v1sim(fileInName.c_str(), {"sim/u_pipos_sim_PLAIN.Q_psi_PLAIN.x1x1",
                                                             "sim/u_pipos_sim_PLAIN.Q_psi_PLAIN.y1y1" });

    v1sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1sim.SetMarker(-1);
    v1sim.SetIsFillLine();
    v1sim.SetPalette(Helper::palette1);
    v1sim.SetBiasPalette(false);
    v1sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                 {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1sim.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1sim.ShiftProjectionAxis(axes.at(kProjection).shift_);
    v1sim.Calculate();

    DoubleDifferentialCorrelation v1rec(fileInName.c_str(), {"rec/u_pipos_rec_RESCALED.Q_psi_PLAIN.x1x1",
                                                             "rec/u_pipos_rec_RESCALED.Q_psi_PLAIN.y1y1" });

    v1rec.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1rec.SetMarker(kFullSquare);
    v1rec.SetPalette(Helper::palette1);
    v1rec.SetBiasPalette(false);
    v1rec.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                 {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1rec.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1rec.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1rec.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1rec.ShiftProjectionAxis(axes.at(kProjection).shift_);
    v1rec.Calculate();

    HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});

    const float text_y = 0.94;
    pic.AddText({0.2, text_y     , "#pi^{+}"}, 0.025);
    pic.AddText({0.2, text_y-0.03, "5M Au+Au"}, 0.025);
    pic.AddText({0.2, text_y-0.06, "DCM-QGSM-SMM"}, 0.025);
    pic.AddText({0.2, text_y-0.09, "12A GeV/c"}, 0.025);
    pic.AddText({0.2, text_y-0.12, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                             " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);

    auto leg = new TLegend(0.55, 0.85, 0.9, 0.95);

//     auto* entry = leg->AddEntry("", "MC input", "F");
//     entry->SetFillColorAlpha(kBlack, 0.2);
//     entry->SetLineColor(kWhite);
//     entry->SetFillStyle(1000);
//
//     entry = leg->AddEntry("", "Reconstructed", "P");
//     entry->SetMarkerSize(2);
//     entry->SetMarkerStyle(kFullSquare);

    auto entry = leg->AddEntry("", "p_{T}, GeV/c", "");

    for( auto obj : v1sim.GetProjections() ){
//       pic.AddDrawable( obj );
    }
    for( auto obj : v1rec.GetProjections() ){
      pic.AddDrawable( obj );
      leg->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
    }

    pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), "v_{1}"});
    pic.CustomizeXRange();
    pic.CustomizeYRange();
    pic.AddLegend(leg);

    pic.Draw();

    if(is_first_canvas)
      pic.GetCanvas()->Print((fileOutName + ".pdf(").c_str(), "pdf");
    else
      pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");
    is_first_canvas = false;

  }

  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");
}
