#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 1000;

using namespace AnalysisTree;

void topo_corr(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("topo_corr.root");

  struct Generation {
    int generation_id_;
    std::string name_;
  };

  std::vector<Generation> generations {
    {0, "bckgr"},
    {1, "sgnl"}
  };

  for(auto& ge : generations) {
    auto cu = EqualsCut("Candidates.generation", ge.generation_id_);
    Cuts* generation_cuts = new Cuts(ge.name_, {cu});
    task->AddH2({"#chi^{2}_{prim #pi^-}", {"Candidates", "chi2_prim_first"}, {nbins, 0, 50}}, {"#chi^{2}_{geo}", {"Candidates", "chi2_geo"}, {nbins, 0, 50}}, generation_cuts);
    task->AddH2({"#chi^{2}_{prim p}", {"Candidates", "chi2_prim_second"}, {nbins, 0, 50}}, {"#chi^{2}_{geo}", {"Candidates", "chi2_geo"}, {nbins, 0, 50}}, generation_cuts);
    task->AddH2({"L/#DeltaL", {"Candidates", "l_over_dl"}, {nbins, 0, 50}}, {"#chi^{2}_{geo}", {"Candidates", "chi2_geo"}, {nbins, 0, 50}}, generation_cuts);
    task->AddH2({"#chi^{2}_{prim #pi^-}", {"Candidates", "chi2_prim_first"}, {nbins, 0, 50}}, {"DCA, cm", {"Candidates", "distance"}, {nbins, 0, 20}}, generation_cuts);
    task->AddH2({"#chi^{2}_{prim p}", {"Candidates", "chi2_prim_second"}, {nbins, 0, 50}}, {"DCA, cm", {"Candidates", "distance"}, {nbins, 0, 20}}, generation_cuts);
    task->AddH2({"L/#DeltaL", {"Candidates", "l_over_dl"}, {nbins, 0, 50}}, {"DCA, cm", {"Candidates", "distance"}, {nbins, 0, 20}}, generation_cuts);
  }

  man->AddTask(task);
  man->SetVerbosityPeriod(100);
  man->Init({filelist}, {"pTree"});
  man->Run(-1);
  man->Finish();
}


int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./topo_corr filelist" << std::endl;
    return -1;
  }

  topo_corr(argv[1]);
  return 0;
}
