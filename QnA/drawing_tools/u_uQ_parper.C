#include "lambda.h"

bool float_equality(float a, float b);
float C_via_b(float b);

void u_uQ_parper() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style.cc" );
  bool is_first_canvas = true;

//   const bool is_draw_difference = false;
  const bool is_draw_difference = true;

  std::string particle = "lambda";
//   std::string particle = "pipos";
//   std::string particle = "pineg";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/covariances_scol." + particle + ".root";

  axes.at(2).reco_name_ = "impactpar";

  TFile* fileIn = TFile::Open(fileName.c_str(), "open");
  auto* dc = (Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>*)fileIn->Get<Qn::DataContainer<Qn::StatCalculate,Qn::Axis<double>>>("upar.all");
  assert(dc!=nullptr);
  for(auto& ax : axes) {
    Qn::Axis<double> qnaxis = dc->GetAxis(ax.sim_name_);
    for(int i=0; i<=qnaxis.size(); i++) {
      ax.bin_edges_.push_back(qnaxis.GetLowerBinEdge(i));
    }
  }

  TF1 oneline("oneline", "1", -100, 100);
  TF1 minusoneline("minusoneline", "-1", -100, 100);

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");

  for(int iEdge=0; iEdge<axes.at(kSelect).bin_edges_.size()-1; iEdge++){

    auto value = DoubleDifferentialCorrelation( fileName.c_str(), {"uperp.all"} );
//     auto value = DoubleDifferentialCorrelation( fileName.c_str(), {"uQperp.all"} );
    value.SetSliceVariable(axes.at(kSlice).title_.c_str(), axes.at(kSlice).unit_.c_str());
    value.Scale(2);
    value.SetMarker(kFullSquare);
    value.SetPalette({kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed,
                      kOrange+1, kBlue, kGreen+2, kAzure-4, kSpring, kViolet, kRed});
    value.Rebin({{axes.at(kSelect).reco_name_.c_str(),
                   {axes.at(kSelect).bin_edges_.at(iEdge), axes.at(kSelect).bin_edges_.at(iEdge+1)}}});
    value.SetProjectionAxis({axes.at(kProjection).reco_name_.c_str(), axes.at(kProjection).bin_edges_});
    value.SetSliceAxis({axes.at(kSlice).reco_name_.c_str(), axes.at(kSlice).bin_edges_});
    value.ShiftSliceAxis(axes.at(kSlice).shift_);
    value.Calculate();
    value.ShiftProjectionAxis(axes.at(kProjection).shift_);
  }



}

bool float_equality(float a, float b) {
  if (std::fabs(a-b)<1e-5)
    return true;
  else
    return false;
}

float C_via_b(float b) {
  if (float_equality(b, 0.0))  return 0.;
  if (float_equality(b, 4.2))  return 10.;
  if (float_equality(b, 6.0))  return 20.;
  if (float_equality(b, 7.3))  return 30.;
  if (float_equality(b, 8.4))  return 40.;
  if (float_equality(b, 9.4))  return 50.;
  if (float_equality(b, 10.3)) return 60.;
  if (float_equality(b, 11.2)) return 70.;
  return -1;
}
