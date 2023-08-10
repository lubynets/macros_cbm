#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "BasicQA.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void example(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("g4_proc_id.root");

  struct Particle {
    std::string name_;
    SimpleCut cut_;
  };

  std::vector<Particle> particles {
    {"lambda", EqualsCut("SimParticles.pid", 3122)},
    {"kshort", EqualsCut("SimParticles.pid", 310)},
    {"xi", EqualsCut("SimParticles.pid", 3312)},
    {"omega", EqualsCut("SimParticles.pid", 3334)},
  };

  task->AddH1({"b, fm", Variable::FromString("SimEventHeader.b"), {100, 0, 20}});

  for(auto& pa : particles) {
    Cuts* pid_cut = new Cuts(pa.name_.c_str(), {pa.cut_});
    task->AddH1({"G4 process ID", Variable::FromString("SimParticles.geant_process_id"), {50, -0.5, 49.5}}, pid_cut);
  }

  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->SetVerbosityPeriod(1000);
  man->Run(-1);
  man->Finish();
}


int main(int argc, char* argv[]){
    if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./g4_proc_id filename" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filename = argv[1];
  example(filename);

  return 0;
}
