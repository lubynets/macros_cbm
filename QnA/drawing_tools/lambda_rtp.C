#include "lambda.h"

void lambda_rtp() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  std::string fileNameWei = "/home/oleksii/cbmdir/working/qna/rectrackspsi/" + evegen + "/cl.rtp." + evegen + ".wei.root";

  std::vector<std::string> particles{
                                     "lambda",
                                     "kshort",
                                    };

  std::string step;
  bool average_comp = true;
//   bool average_comp = false;

  float y_lo, y_hi;

  SetAxis("centrality", "select");
  SetAxis("rapidity", "projection");
  SetAxis("pT", "slice");

  std::string y_axis_title = "v_{1}";

  axes.at(0).sim_name_ = "SimParticles_pT";
  axes.at(1).sim_name_ = "SimParticles_rapidity";
  axes.at(2).sim_name_ = "RecEventHeader_centrality_tracks";
  axes.at(0).reco_name_ = "ReconstructedParticles_pT";
  axes.at(1).reco_name_ = "ReconstructedParticles_rapidity";
  axes.at(2).reco_name_ = "RecEventHeader_centrality_tracks";

  TFile* fileIn = TFile::Open(fileNameWei.c_str(), "open");
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

  for(auto& particle : particles) {

    if(evegen == "dcmqgsm") {
      if(particle == "lambda") {
        SetSliceAxisBinEdges({0, 0.4, 0.8, 1.2, 1.6});
        y_lo = -0.2;
        y_hi = 0.2;
      }
      if(particle == "kshort") {
        SetSliceAxisBinEdges({0, 0.4, 0.8, 1.6});
        y_lo = -0.04;
        y_hi = 0.04;
      }
    }
    if(evegen == "urqmd") {
      if(particle == "lambda") {
        SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});
        y_lo = -0.2;
        y_hi = 0.2;
      }
      if(particle == "kshort") {
        SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});
        y_lo = -0.2;
        y_hi = 0.2;
      }
    }

    bool is_first_canvas = true;
    std::string fileOutName;

    fileOutName = "v1." + particle;

    if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

    for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {
      std::vector<DoubleDifferentialCorrelation> v1_rec;
      if(average_comp) {
        v1_rec.resize(1);
        v1_rec.at(0) = DoubleDifferentialCorrelation( fileNameWei.c_str(), {("rec/u_" + particle + "_rec_RESCALED.Q_psi_PLAIN.x1x1").c_str(),
                                                                            ("rec/u_" + particle + "_rec_RESCALED.Q_psi_PLAIN.y1y1").c_str()} );
        v1_rec.at(0).SetMarker(kFullSquare);
      } else {
        v1_rec.resize(2);
        v1_rec.at(0) = DoubleDifferentialCorrelation( fileNameWei.c_str(), {("rec/u_" + particle + "_rec_RESCALED.Q_psi_PLAIN.x1x1").c_str()} );
        v1_rec.at(1) = DoubleDifferentialCorrelation( fileNameWei.c_str(), {("rec/u_" + particle + "_rec_RESCALED.Q_psi_PLAIN.y1y1").c_str()} );
        v1_rec.at(0).SetMarker(kFullSquare);
        v1_rec.at(1).SetMarker(kOpenSquare);
      }

      for(auto& vc : v1_rec) {
        vc.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
        vc.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                       kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
        vc.SetBiasPalette(false);
        vc.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                  {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
        vc.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
        vc.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
        vc.ShiftSliceAxis(axes.at(kSlice).shift_);
        vc.Calculate();
        vc.ShiftProjectionAxis(axes.at(kProjection).shift_);
      }

      v1_rec.at(0).SlightShiftProjectionAxis(0.025);
      if(!average_comp) {
        v1_rec.at(1).SlightShiftProjectionAxis(0.025, 0.0125);
      }

      auto v1_sim = DoubleDifferentialCorrelation( fileNameWei.c_str(), {("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.x1x1").c_str(),
                                                                         ("sim/u_" + particle + "_sim_PLAIN.Q_psi_PLAIN.y1y1").c_str()} );
      v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      v1_sim.SetMarker(-1);
      v1_sim.SetIsFillLine();
      v1_sim.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                         kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed});
      v1_sim.SetBiasPalette(false);
      v1_sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                      {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
      v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
      v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
      v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
      v1_sim.Calculate();
      v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);

      HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});
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

      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);
      leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

      TLegendEntry* entry;
      entry = leg1->AddEntry("", "Sim, x&y averaged", "L");
      entry->SetMarkerSize(2);
      if(average_comp) {
        entry = leg1->AddEntry("", "Reco, x&y averaged", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kFullSquare);
      } else {
        entry = leg1->AddEntry("", "Sim, X", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kFullSquare);
        entry = leg1->AddEntry("", "Sim, Y", "P");
        entry->SetMarkerSize(2);
        entry->SetMarkerStyle(kOpenSquare);
      }

      for(int i=0; i<v1_rec.size(); i++) {
        for( auto obj : v1_rec.at(i).GetProjections() ) {
          pic.AddDrawable( obj );
          if(i==0) leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
        }
      }
      for( auto obj : v1_sim.GetProjections() ) {
        pic.AddDrawable( obj );
      }
      pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), y_axis_title});

      pic.CustomizeXRange();
//       pic.CustomizeYRange();
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
    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((fileOutName + ".pdf]").c_str(), "pdf");

    if(is_write_rootfile) fileOut->Close();
  }
}
