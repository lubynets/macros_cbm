#include "lambda.h"

void lambda_ffm_dv1dy(int iSetup=1, int iPdg=1, int iFit=2, int iPm=2) {
  bool verbose{false};
//   bool verbose{true};

  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );

  std::string evegen, pbeam, pdg, particle, setup, pol, pm;

  if(iSetup==1) {evegen = "dcmqgsm"; pbeam = "12";  setup = "d12";}
  if(iSetup==2) {evegen = "dcmqgsm"; pbeam = "3.3"; setup = "d3";  axes.at(1).shift_ = -0.985344;}
  if(iSetup==3) {evegen = "urqmd";   pbeam = "12";  setup = "u12";}

  if(iPdg==1) {pdg = "3122"; particle = "lambda";}
  if(iPdg==2) {pdg = "310";  particle = "kshort";}

  if(iFit==1) pol = "pol1";
  if(iFit==2) pol = "pol3";

  if(iPm==1) pm = "pm1";
  if(iPm==2) pm = "pm2";

  bool is_write_rootfile = false;
//   bool is_write_rootfile = true;

  Qn::Stat::ErrorType mean_mode{Qn::Stat::ErrorType::PROPAGATION};
  Qn::Stat::ErrorType error_mode{Qn::Stat::ErrorType::BOOTSTRAP};

  std::string fileDv1DyName = "/home/oleksii/cbmdir/working/qna/flowfrommodel/cl.dv1dy_" + pol + "." + pm + ".ffm." + setup + ".root";

  SetAxis("centrality", "projection");
  SetAxis("pT", "slice");
//   SetAxis("centrality", "slice");
//   SetAxis("pT", "projection");

  axes.at(0).sim_name_ = "SimParticles_pT";
  axes.at(1).sim_name_ = "SimParticles_rapidity";
  axes.at(2).sim_name_ = "SimEventHeader_centrality_impactpar";
  axes.at(0).reco_name_ = "SimParticles_pT";
  axes.at(1).reco_name_ = "SimParticles_rapidity";
  axes.at(2).reco_name_ = "SimEventHeader_centrality_impactpar";

  std::string y_axis_title;

  std::vector<std::string> components{"x1x1", "y1y1"};
//   std::vector<std::string> components{"ave"};

  std::vector<std::string> fitcoeffs{"slope"/*, "intercept"*/};

  TFile* fileDv1Dy = TFile::Open(fileDv1DyName.c_str(), "open");
  auto* dc = (QnDcSD*)fileDv1Dy->Get<QnDcSD>(("uPsi/slope/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.ave").c_str());
  if(dc == nullptr) throw std::runtime_error("dc is nullptr");

  if(verbose) std::cout << "Start axes mark up\n";
  for(auto& ax : axes) {
    if(ax.name_ == "rapidity") continue;
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }
  if(verbose) std::cout << "Finish axes mark up\n";

  float slightprojshift;
  if(axes.at(kSlice).name_ == "centrality") {
    SetSliceAxisBinEdges({0, 5, 20, 40, 70});
    slightprojshift = 0.015;
  } else {
    slightprojshift = 1;
  }

  TFile* fileOut;

  for(auto& fc : fitcoeffs) {
    if(verbose) std::cout << "fc = " << fc << "\n";

    std::string fileOutName = fc + "." + pol + "." + pm + "." + pdg + "." + setup;
    if(is_write_rootfile) fileOut = TFile::Open((fileOutName + ".root").c_str(), "recreate");

    if(fc == "slope") y_axis_title = "dv_{1}/dy";
    if(fc == "intercept") y_axis_title = "v_{1}|_{y=0}";

    auto v1_sim = DoubleDifferentialCorrelation( fileDv1DyName.c_str(), {("uPsi/" + fc + "/u_sim_" + particle + "_PLAIN.Q_psi_PLAIN.ave").c_str()} );
    if(verbose) std::cout << "v1_sim created\n";
    v1_sim.SetErrorType(error_mode);
    v1_sim.SetMeanType(mean_mode);
    v1_sim.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    v1_sim.SetMarker(-1);
    v1_sim.SetIsFillLine(0.4);
    v1_sim.SetPalette(Helper::palette1);
    v1_sim.SetBiasPalette(false);
    v1_sim.SetProjectionAxis({axes.at(kProjection).sim_name_.c_str(), axes.at(kProjection).bin_edges_});
    v1_sim.SetSliceAxis({axes.at(kSlice).sim_name_.c_str(), axes.at(kSlice).bin_edges_});
    v1_sim.ShiftSliceAxis(axes.at(kSlice).shift_);
    v1_sim.ShiftProjectionAxis(axes.at(kProjection).shift_);
    v1_sim.SetSliceAxisPrecision(axes.at(kSlice).precision_);
    if(verbose) std::cout << "v1_sim start Calculate\n";
    v1_sim.Calculate();
    if(verbose) std::cout << "v1_sim calculated\n";

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
    }

    auto leg1 = new TLegend();
    leg1->SetBorderSize(1);

    leg1->AddEntry("", (axes.at(kSlice).title_+axes.at(kSlice).unit_ + ":").c_str(), "");

    for( auto obj : v1_sim.GetProjections() ){
      pic.AddDrawable( obj );
      obj->GetPoints()->SetLineWidth(0);
      leg1->AddEntry( obj->GetPoints(), obj->GetTitle().c_str(), "F" );
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
      pic.Write("heap_picture");
      pic.GetCanvas()->Write("canvas");
    }

    pic.GetCanvas()->Print((fileOutName + ".pdf").c_str(), "pdf");

    if(is_write_rootfile) fileOut->Close();
  }// fitcoeffs
}
