#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 10000;
const float HugeValue = 1e9;
const float TinyValue = 2e-5;
const float PI = TMath::Pi();

using namespace AnalysisTree;

void chi2_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("chi2_qa.root");
  
//   SimpleCut sgnlcut = RangeCut("Candidates.generation",  0.9, 12.1);
//   SimpleCut bckgrcut = RangeCut("Candidates.generation", -0.1, 0.1);
  
//   Cuts* SgnlCut = new Cuts("signal", {sgnlcut});
  
//   task->AddH1({"#chi^{2}_{prim_first}", Variable::FromString("Candidates.chi2_prim_first"), {nbins, -400, 400}});
//   task->AddH1({"#chi^{2}_{prim_second}", Variable::FromString("Candidates.chi2_prim_second"), {nbins, -400, 400}});
//   task->AddH1({"#chi^{2}_{geo}", Variable::FromString("Candidates.chi2_geo"), {nbins, -400, 400}});
//   
//   task->AddH1({"#chi^{2}_{prim_first}", Variable::FromString("Candidates.chi2_prim_first"), {nbins, 0, 4e7}});
//   task->AddH1({"#chi^{2}_{prim_second}", Variable::FromString("Candidates.chi2_prim_second"), {nbins, 0, 5e7}});
//   task->AddH1({"#chi^{2}_{geo}", Variable::FromString("Candidates.chi2_geo"), {nbins, 0, 1e7}});
//   
//   task->AddH1({"#chi^{2}_{prim_first}", Variable::FromString("Candidates.chi2_prim_first"), {nbins, 0, 6e6}}, SgnlCut);
//   task->AddH1({"#chi^{2}_{prim_second}", Variable::FromString("Candidates.chi2_prim_second"), {nbins, 0, 4e5}}, SgnlCut);
//   task->AddH1({"#chi^{2}_{geo}", Variable::FromString("Candidates.chi2_geo"), {nbins, 0, 400}}, SgnlCut);
  
  Variable det_prim_first("dets_prim_first", {{"Candidates", "chi2_prim_first_det"}, {"Candidates", "chi2_prim_first_detinv"}}, []( std::vector<double>& var ) { return var.at(0)*var.at(1); });
  Variable det_prim_second("dets_prim_second", {{"Candidates", "chi2_prim_second_det"}, {"Candidates", "chi2_prim_second_detinv"}}, []( std::vector<double>& var ) { return var.at(0)*var.at(1); });
  Variable det_geo("dets_geo", {{"Candidates", "chi2_geo_det"}, {"Candidates", "chi2_geo_detinv"}}, []( std::vector<double>& var ) { return var.at(0)*var.at(1); });
  
//   task->AddH1({"detC*detC^{-1}_{prim_first}", det_prim_first, {nbins, -0.5, 1.5}});
//   task->AddH1({"detC*detC^{-1}_{prim_second}", det_prim_second, {nbins, -0.5, 1.5}});
//   task->AddH1({"detC*detC^{-1}_{geo}", det_geo, {nbins, -0.5, 1.5}});
   
  SimpleCut det_prim_first_zerocut = RangeCut(det_prim_first, -TinyValue, TinyValue);
  SimpleCut det_prim_second_zerocut = RangeCut(det_prim_second, -TinyValue, TinyValue);
  SimpleCut det_geo_zerocut = RangeCut(det_geo, -TinyValue, TinyValue);
  
  Cuts* DetPrimFirstZerocut = new Cuts("DetPrimFirstZerocut", {det_prim_first_zerocut});
  Cuts* DetPrimSecondZerocut = new Cuts("DetPrimSecondZerocut", {det_prim_second_zerocut});
  Cuts* DetGeoZerocut = new Cuts("DetGeoZerocut", {det_geo_zerocut});
  
//   task->AddH1({"#chi^{2}_{prim_first_cov_xx}", Variable::FromString("Candidates.chi2_prim_first_cov_xx"), {nbins, -0.01, 0.01}}, DetPrimFirstZerocut);
//   task->AddH1({"#chi^{2}_{prim_first_cov_xy}", Variable::FromString("Candidates.chi2_prim_first_cov_xy"), {nbins, -0.01, 0.01}}, DetPrimFirstZerocut);
//   task->AddH1({"#chi^{2}_{prim_first_cov_yy}", Variable::FromString("Candidates.chi2_prim_first_cov_yy"), {nbins, -0.01, 0.01}}, DetPrimFirstZerocut);
//   task->AddH1({"#chi^{2}_{prim_first_cov_xz}", Variable::FromString("Candidates.chi2_prim_first_cov_xz"), {nbins, -0.01, 0.01}}, DetPrimFirstZerocut);
//   task->AddH1({"#chi^{2}_{prim_first_cov_yz}", Variable::FromString("Candidates.chi2_prim_first_cov_yz"), {nbins, -0.01, 0.01}}, DetPrimFirstZerocut);
//   task->AddH1({"#chi^{2}_{prim_first_cov_zz}", Variable::FromString("Candidates.chi2_prim_first_cov_zz"), {nbins, -0.01, 0.01}}, DetPrimFirstZerocut);
//   
//   task->AddH1({"#chi^{2}_{prim_second_cov_xx}", Variable::FromString("Candidates.chi2_prim_second_cov_xx"), {nbins, -0.01, 0.01}}, DetPrimSecondZerocut);
//   task->AddH1({"#chi^{2}_{prim_second_cov_xy}", Variable::FromString("Candidates.chi2_prim_second_cov_xy"), {nbins, -0.01, 0.01}}, DetPrimSecondZerocut);
//   task->AddH1({"#chi^{2}_{prim_second_cov_yy}", Variable::FromString("Candidates.chi2_prim_second_cov_yy"), {nbins, -0.01, 0.01}}, DetPrimSecondZerocut);
//   task->AddH1({"#chi^{2}_{prim_second_cov_xz}", Variable::FromString("Candidates.chi2_prim_second_cov_xz"), {nbins, -0.01, 0.01}}, DetPrimSecondZerocut);
//   task->AddH1({"#chi^{2}_{prim_second_cov_yz}", Variable::FromString("Candidates.chi2_prim_second_cov_yz"), {nbins, -0.01, 0.01}}, DetPrimSecondZerocut);
//   task->AddH1({"#chi^{2}_{prim_second_cov_zz}", Variable::FromString("Candidates.chi2_prim_second_cov_zz"), {nbins, -0.01, 0.01}}, DetPrimSecondZerocut);
//   
//   task->AddH1({"#chi^{2}_{geo_cov_xx}", Variable::FromString("Candidates.chi2_prim_geo_cov_xx"), {nbins, -0.01, 0.01}}, DetGeoZerocut);
//   task->AddH1({"#chi^{2}_{geo_cov_xy}", Variable::FromString("Candidates.chi2_prim_geo_cov_xy"), {nbins, -0.01, 0.01}}, DetGeoZerocut);
//   task->AddH1({"#chi^{2}_{geo_cov_yy}", Variable::FromString("Candidates.chi2_prim_geo_cov_yy"), {nbins, -0.01, 0.01}}, DetGeoZerocut);
//   task->AddH1({"#chi^{2}_{geo_cov_xz}", Variable::FromString("Candidates.chi2_prim_geo_cov_xz"), {nbins, -0.01, 0.01}}, DetGeoZerocut);
//   task->AddH1({"#chi^{2}_{geo_cov_yz}", Variable::FromString("Candidates.chi2_prim_geo_cov_yz"), {nbins, -0.01, 0.01}}, DetGeoZerocut);
//   task->AddH1({"#chi^{2}_{geo_cov_zz}", Variable::FromString("Candidates.chi2_prim_geo_cov_zz"), {nbins, -0.01, 0.01}}, DetGeoZerocut);
  
  
//   task->AddH1({"detC_{first}", Variable::FromString("Candidates.chi2_prim_first_det"), {nbins, -TinyValue, TinyValue}});
//   task->AddH1({"detC_{second}", Variable::FromString("Candidates.chi2_prim_second_det"), {nbins, -TinyValue, TinyValue}});
//   task->AddH1({"detC_{geo}", Variable::FromString("Candidates.chi2_geo_det"), {nbins, -TinyValue, TinyValue}});
//   
//   task->AddH1({"detC^{-1}_{first}", Variable::FromString("Candidates.chi2_prim_first_detinv"), {nbins, -HugeValue, HugeValue}});
//   task->AddH1({"detC^{-1}_{second}", Variable::FromString("Candidates.chi2_prim_second_detinv"), {nbins, -HugeValue, HugeValue}});
//   task->AddH1({"detC^{-1}_{geo}", Variable::FromString("Candidates.chi2_geo_detinv"), {nbins, -HugeValue, HugeValue}});
 
  task->AddH1({"detC_{first}", Variable::FromString("Candidates.chi2_prim_first_det"), {nbins, -1e-12, 1e-12}}, DetPrimFirstZerocut);
  task->AddH1({"detC_{second}", Variable::FromString("Candidates.chi2_prim_second_det"), {nbins, -1e-12, 1e-12}}, DetPrimSecondZerocut);
  task->AddH1({"detC_{geo}", Variable::FromString("Candidates.chi2_geo_det"), {nbins, -1e-12, 1e-12}}, DetGeoZerocut);
  
  task->AddH1({"detC^{-1}_{first}", Variable::FromString("Candidates.chi2_prim_first_detinv"), {nbins, -1e-12, 1e-12}}, DetPrimFirstZerocut);
  task->AddH1({"detC^{-1}_{second}", Variable::FromString("Candidates.chi2_prim_second_detinv"), {nbins, -1e-12, 1e-12}}, DetPrimSecondZerocut);
  task->AddH1({"detC^{-1}_{geo}", Variable::FromString("Candidates.chi2_geo_detinv"), {nbins, -1e-12, 1e-12}}, DetGeoZerocut);
  
  man->AddTask(task);
  
  man->SetVerbosityPeriod(8);

  man->Init({filelist}, {"pTree"});
  man->Run(-1);
  man->Finish();

}

int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./chi2_qa filelist" << std::endl;
    return -1;
  }  
  
  chi2_qa(argv[1]);
  return 0;
}