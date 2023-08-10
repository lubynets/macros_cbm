#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

#include "TMath.h"

const int nbins = 1000;

using namespace AnalysisTree;

void generation_debug(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("generation_debug.root");

  task->AddH1({("m_{" + qp.axisname_ + "}, GeV/c^{2}").c_str(),                                Variable::FromString("Candidates.mass"),             {nbins, qp.leftrange_, qp.rightrange_}}, selection_cuts);




  man->AddTask(task);

  man->Init({filelist}, {"pTree"});
  man->Run(-1);
  man->Finish();

}

int main(int argc, char* argv[]){
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./generation_debug filelist" << std::endl;
    return -1;
  }

  generation_debug(argv[1]);
  return 0;
}
