#include "lambda.h"

void lambda_ffm(int iSetup=1, int iParticle=1, int iFit=1) {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen, pbeam, setup, particle;

  if(iSetup==1) {evegen = "dcmqgsm"; pbeam = "12";  setup = "d12";}
  if(iSetup==2) {evegen = "dcmqgsm"; pbeam = "3.3"; setup = "d3";  axes.at(1).shift_ = -0.985344;}
  if(iSetup==3) {evegen = "urqmd";   pbeam = "12";  setup = "u12";}

  if(iParticle==1) particle = "lambda";
  if(iParticle==2) particle = "kshort";

  bool draw_fit;
  std::string pol;
  std::string proj;

  if(iFit==1) {draw_fit = false; pol = "";}
  if(iFit==2) {draw_fit = true ; pol = "pol1";}
  if(iFit==3) {draw_fit = true ; pol = "pol3";}

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileName = "/home/oleksii/cbmdir/working/qna/flowfrommodel/cl.ffm." + setup + ".root";
  std::string fileDv1DyName = "/home/oleksii/cbmdir/working/qna/flowfrommodel/cl.dv1dy_" + pol + ".ffm." + setup + ".root";

  SetAxis("centrality", "select");

//   SetAxis("rapidity", "projection"); proj = "y";
//   SetAxis("pT", "slice");

  SetAxis("rapidity", "slice");  proj = "pt";
  SetAxis("pT", "projection");

  axes.at(0).sim_name_ = "SimParticles_pT";
  axes.at(1).sim_name_ = "SimParticles_rapidity";
  axes.at(2).sim_name_ = "SimEventHeader_centrality_impactpar";
  axes.at(0).reco_name_ = "SimParticles_pT";
  axes.at(1).reco_name_ = "SimParticles_rapidity";
  axes.at(2).reco_name_ = "SimEventHeader_centrality_impactpar";

  TFile* fileMc = TFile::Open(fileName.c_str(), "open");
  if(fileMc == nullptr) throw std::runtime_error("fileMc == nullptr");
  auto* dc = (QnDcSCo*)fileMc->Get<QnDcSCo>(("uPsi/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.x1x1").c_str());
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  float y_lim;
  if(axes.at(kProjection).name_ == "rapidity") {
    if(setup == "d12") {
      if(particle == "lambda") {SetSliceAxisBinEdges({0, 0.4, 0.8, 1.2, 1.6}); y_lim = 0.15;}
      if(particle == "kshort") {SetSliceAxisBinEdges({0, 0.4, 0.8, 1.6});      y_lim = 0.05;}
    }
    if(setup == "d3") {
      if(particle == "lambda") {SetSliceAxisBinEdges({0, 0.4, 0.8, 1.6});      y_lim = 0.13;}
      if(particle == "kshort") {SetSliceAxisBinEdges({0, 0.4, 0.8, 1.6});      y_lim = 0.13;}
    }
    if(setup == "u12") {
      if(particle == "lambda") {SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});      y_lim = 0.08;}
      if(particle == "kshort") {SetSliceAxisBinEdges({0, 0.8, 1.2, 1.6});      y_lim = 0.12;}
    }
    y_lim *= 2;
  } else {
    const float midrap = -axes.at(kSlice).shift_ ;
    SetSliceAxisBinEdges({0+midrap, 0.4+midrap, 0.8+midrap, 1.2+midrap, 2.0+midrap});
  }

  TFile* fileDv1Dy = TFile::Open(fileDv1DyName.c_str(), "open");
  if(draw_fit && fileDv1Dy == nullptr) throw std::runtime_error("fileDv1Dy == nullptr");

  bool is_first_canvas = true;
  if(pol!="") pol+=".";
  std::string fileOutName = "v1_" + proj + ".ffm." + particle + "." + setup + "." + pol;
  if(pol.size() == 5) pol.pop_back();

  auto* sim_slope = draw_fit ? (QnDcSD*)fileDv1Dy->Get<QnDcSD>(("uPsi/slope/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.ave").c_str()) : nullptr;
  auto* sim_offset = draw_fit ? (QnDcSD*)fileDv1Dy->Get<QnDcSD>(("uPsi/intercept/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.ave").c_str()) : nullptr;
  auto* sim_third = draw_fit ? (QnDcSD*)fileDv1Dy->Get<QnDcSD>(("uPsi/third/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.ave").c_str()) : nullptr;

  if(draw_fit && (sim_slope==nullptr || sim_offset==nullptr)) throw std::runtime_error("sim_slope==nullptr || sim_offset==nullptr");
  if(draw_fit && pol == "pol3" && sim_third == nullptr) throw std::runtime_error("sim_third == nullptr");

  bool select_then_slice;
  if(draw_fit) {
    if(sim_slope->GetAxes().at(0).Name() == axes.at(kSelect).sim_name_.c_str() &&
      sim_slope->GetAxes().at(1).Name() == axes.at(kSlice).sim_name_.c_str()) {
      select_then_slice = true;
    } else if (sim_slope->GetAxes().at(1).Name() == axes.at(kSelect).sim_name_.c_str() &&
              sim_slope->GetAxes().at(0).Name() == axes.at(kSlice).sim_name_.c_str()) {
      select_then_slice  = false;
    } else throw std::runtime_error("Something wrong with Select and Slice axes in dv1dy");
  }

  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++) {

    auto v1_sim = DoubleDifferentialCorrelation( fileName.c_str(),
                                                {("uPsi/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.x1x1").c_str(),
                                                 ("uPsi/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.y1y1").c_str() } );
    v1_sim.SetErrorType(error_mode);
    v1_sim.SetMeanType(mean_mode);
    v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1_sim.SetMarker(-1);
    v1_sim.SetIsFillLine(0.4);
    v1_sim.SetPalette(Helper::palette1);
    v1_sim.SetBiasPalette(false);
    v1_sim.Rebin({{axes.at(kSelect).sim_name_.c_str(),
                  {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);
    v1_sim.Calculate();

    HeapPicture pic( (axes.at(kSelect).name_ + "_" + std::to_string(iEdge)).c_str(), {1000, 1000});

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
    pic.AddText({0.2, 0.78, (axes.at(kSelect).title_ + ": " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) +
                            " - " + to_string_with_precision(axes.at(kSelect).bin_edges_.at(iEdge+1) + axes.at(kSelect).shift_, axes.at(kSelect).precision_) + axes.at(kSelect).unit_).c_str()}, 0.025);

    auto leg1 = new TLegend(/*0.4, 0.8, 0.6, 0.95*/);

    leg1->SetHeader((axes.at(kSlice).title_+axes.at(kSlice).unit_).c_str());

    for( auto obj : v1_sim.GetProjections() ){
      pic.AddDrawable( obj );
      leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "L" );
    }

    if(draw_fit) {
      for(unsigned long iSlice=0; iSlice<v1_sim.GetProjections().size(); iSlice++) {
        float k, b, p3;
        const unsigned long i1 = select_then_slice ? iEdge : iSlice;
        const unsigned long i2 = select_then_slice ? iSlice : iEdge;
        k = sim_slope->At({i1, i2}).Mean();
        b = sim_offset->At({i1, i2}).Mean();
        p3 = pol=="pol3" ? sim_third->At({i1, i2}).Mean() : 0;
        TF1* fsim = new TF1("fsim", "[0]+[1]*x+[2]*x*x*x", -10, 10);
        fsim->SetParameters(b, k, p3);
        fsim->SetLineStyle(1);
        fsim->SetLineColor(Helper::palette1.at(iSlice));
        pic.AddFunction(fsim);
      }
    }

    pic.SetAxisTitles({(axes.at(kProjection).title_ + axes.at(kProjection).unit_).c_str(), "v_{1}"});

//     pic.SetXRange({-1.05, 1.05});
    pic.CustomizeXRange();
    pic.CustomizeYRange();
//     pic.CustomizeYRangeWithLimits(-y_lim, y_lim);
    pic.AddLegend(leg1);
    pic.SetIsCustomizeLegend();
//         pic.SetGridXY();
    pic.Draw();

    if(is_first_canvas)
      pic.GetCanvas()->Print((fileOutName + "pdf(").c_str(), "pdf");
    else
      pic.GetCanvas()->Print((fileOutName + "pdf").c_str(), "pdf");
    is_first_canvas = false;

  }

  TCanvas emptycanvas("", "", 1000, 1000);
  emptycanvas.Print((fileOutName + "pdf]").c_str(), "pdf");

}
