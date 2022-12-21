#include "shape_fit_qa.h"

void flow_fit_qa() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string fileName = "/home/oleksii/cbmdir/working/qna/correlations/aXmass/out.fitter.dcmqgsm.apr20.recpid.lightcuts1.3122.root";
  std::string particle = "#Lambda";

  std::vector<std::pair<std::string, std::string>> values {
    {"Params", "signal.x1x1"},
    {"Params", "signal.y1y1"},
    {"Params", "bckgr_0.x1x1"},
    {"Params", "bckgr_0.y1y1"},
    {"Params", "bckgr_1.x1x1"},
    {"Params", "bckgr_1.y1y1"},
    {"Chi2s", "fit_chi2ndf.x1x1"},
    {"Chi2s", "fit_chi2ndf.y1y1"},
    {"Entries", "entries_sgnl"},
    {"Entries", "entries_bckgr"}
  };

  SetAxis("centrality", "select");
  SetAxis("y", "projection");
  SetAxis("pT", "slice");

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>>("Params/signal.x1x1");
  if(dc==nullptr){
    throw std::runtime_error("Data container is nullptr");
  }
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }
  fileIn->Close();

  for(auto& va : values) {
    bool is_first_canvas = true;

    std::string dcName = va.first + "/" + va.second;

    for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){
      auto param = DoubleDifferentialCorrelation( fileName.c_str(), {dcName.c_str()});
      param.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
      param.SetMarker(kFullSquare);
      param.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed});
      param.SetBiasPalette(false);
      param.Rebin({{axes.at(kSelect).name_.c_str(),
                      {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
      param.SetProjectionAxis({axes.at(kProjection).name_.c_str(), axes.at(kProjection).bin_edges_});
      param.SetSliceAxis({axes.at(kSlice).name_.c_str(), axes.at(kSlice).bin_edges_});
      param.ShiftSliceAxis(axes.at(kSlice).shift_);
      param.Calculate();
      param.ShiftProjectionAxis(axes.at(kProjection).shift_);
      param.SlightShiftProjectionAxis(0.03);

      HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});

      pic.AddText({0.2, 0.90, particle.c_str()}, 0.035);
      pic.AddText({0.2, 0.87, "5M Au+Au"}, 0.025);
      pic.AddText({0.2, 0.84, "DCM-QGSM-SMM"}, 0.025);
      pic.AddText({0.2, 0.81, "12A GeV/c"}, 0.025);

      pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                              " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);

      auto leg1 = new TLegend();
      leg1->SetBorderSize(1);
      leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

      for( auto obj : param.GetProjections() ){
        pic.AddDrawable( obj );
        leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "P" );
      }

      pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), va.second.c_str()});

      pic.CustomizeXRange();
      pic.CustomizeYRange(0.95);
      pic.AddLegend(leg1);
      pic.CustomizeLegend(leg1);
      pic.Draw();

      if(is_first_canvas)
        pic.GetCanvas()->Print((va.second + ".pdf(").c_str(), "pdf");
      else
        pic.GetCanvas()->Print((va.second + ".pdf").c_str(), "pdf");
      is_first_canvas = false;
    }
    TCanvas emptycanvas("", "", 1000, 1000);
    emptycanvas.Print((va.second + ".pdf)").c_str(), "pdf");
  }
}
