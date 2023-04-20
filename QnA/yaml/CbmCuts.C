#include "CutsRegistry.hpp"

int CbmCuts() {
  using AnalysisTree::RegisterCuts;
  using namespace AnalysisTree::Version1;

  {
    std::string branch = "RecEventHeader";

    std::vector<SimpleCut> cuts;
    SimpleCut vtx_x_cut = RangeCut((branch + ".vtx_x").c_str(), -0.5, 0.5);
    SimpleCut vtx_y_cut = RangeCut((branch + ".vtx_y").c_str(), -0.5, 0.5);
    SimpleCut vtx_z_cut = RangeCut((branch + ".vtx_z").c_str(), -0.03, 0.03);
    SimpleCut vtx_chi2_cut = RangeCut((branch + ".vtx_chi2").c_str(), 0.8, 1.7);

    const char *cuts_name = "goodevents";
    RegisterCuts(cuts_name, Cuts("CbmGoodEvent", {
        vtx_x_cut,
        vtx_y_cut,
        vtx_z_cut,
        vtx_chi2_cut}));
  }

  return 0;
}
