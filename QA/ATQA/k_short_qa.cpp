#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/TaskManager.hpp>

#include "EntryConfig.hpp"
#include "Task.hpp"
#include "BasicQA.hpp"

using namespace AnalysisTree;

const std::string sim_particles = "SimParticles";
const std::string sim_event_header = "SimEventHeader";

int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cout << "Not enough arguments! Please use:" << std::endl;
    std::cout << "   ./k_short_qa filelist" << std::endl;
    return -1;
  }
  int n_events = -1;
  
  if (argc == 3) {
    n_events = atoi(argv[2]);
  }

  const std::string filelist = argv[1];

  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("k_short_qa.root");
  
  auto k_short_cut = EqualsCut("SimParticles.pid", 310);
  auto prim_cut = EqualsCut("SimParticles.mother_id", -1); //-1 for primaries and any other value for secondaries. remove this cut for all lambdas
  auto z_cut = RangeCut("SimParticles.z",-10,80);
  //auto multplcty_cut = RangeCut("RecEventHeader.M", 200, 400);
  //Cuts* mc_lambda_primiry = new Cuts("McLambda" + std::to_string(200) + "_" + std::to_string(400), {lambda_cut, z_cut, multplcty_cut});
  //AddParticleQA(task, sim_particles, mc_lambda_primiry);
  //task->AddH1({"Rec M", {"RecEventHeader", "M"}, {100, 0, 1000}});

for (int i=200;i<=210;i+=50){

    //auto prim_cut = RangeCut("SimParticles.mother_id", -1, 1000000); //Cuts* centr = new Cuts("SemiCentral", {centr_cut});//     Cuts* multp = new Cuts("Multp", {multplcty_cut});    //task->SetEventCuts(centr);//  task->SetEventCuts(multp);
    //auto multplcty_cut = RangeCut("RecEventHeader.M", 200, 400);
    Cuts* mc_k_short_primiry = new Cuts("McK_short" , {k_short_cut,prim_cut , z_cut});
    
    AddParticleQA(task, sim_particles, mc_k_short_primiry);
  }
  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->Run(n_events);
  man->Finish();

  return 0;
}


