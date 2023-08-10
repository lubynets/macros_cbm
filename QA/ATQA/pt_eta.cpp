
#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void pt_eta(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("pt_eta.root");

  const double midrapidity = 0.9853;

  SimpleCut pdg_cut = EqualsCut("RecParticles.pid", 211);
  SimpleCut y_separator = SimpleCut({"RecParticles.rapidity"}, [midrapidity] (std::vector<double>& var) {return std::fabs(var.at(0)-midrapidity - 0. ) > 0.02 &&
                                                                                                                std::fabs(var.at(0)-midrapidity - 0.4) > 0.02 &&
                                                                                                                std::fabs(var.at(0)-midrapidity - 0.8) > 0.02 &&
                                                                                                                std::fabs(var.at(0)-midrapidity - 1.2) > 0.02 &&
                                                                                                                std::fabs(var.at(0)-midrapidity - 1.6) > 0.02 &&
                                                                                                                std::fabs(var.at(0)-midrapidity - 2.0) > 0.02; });
  Cuts* all_cuts = new Cuts("all_cuts", {pdg_cut, y_separator});

  task->AddH2({"#eta", Variable::FromString("RecParticles.eta"), {500, 0, 5}}, {"p_{T}, GeV/c", Variable::FromString("RecParticles.pT"), {400, 0, 2}}, all_cuts);

  man->AddTask(task);

  man->Init({filelist}, {"aTree"});
  man->SetVerbosityPeriod(100);
  man->Run(-1);
  man->Finish();
}

int main(int argc, char* argv[]){
    if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./pt_eta filename" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filename = argv[1];
  pt_eta(filename);

  return 0;
}
