#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 1000;

using namespace AnalysisTree;

void threebranches_qa(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("threebranches_qa.root");

  const float scaleTOF = 29.9792458;
  Variable m2_calc("TofHits_m2", {{"TofHits", "t"}, {"RecEventHeader", "T0"}, {"TofHits", "l"}, {"VtxTracks", "p"}}, [scaleTOF](std::vector<double>& var){
        const float time = var[0]-var[1];
        const float beta = var[2] / (time * scaleTOF);
        return var[3]*var[3] * ( 1. / (beta*beta) -1. ); });

  task->AddH2({"m_{calc}^{2}, GeV^{2}/c^{4}", m2_calc, {500, -1, 2}},
         {"m_{ATConv}^{2}, GeV^{2}/c^{4}", {"TofHits", "mass2"}, {500, -1, 2}});

  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->Run(-1);
  man->Finish();

}

int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./threebranches_qa filelist" << std::endl;
    return -1;
  }

  threebranches_qa(argv[1]);
  return 0;
}
