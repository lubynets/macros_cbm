
#include <string>

#include "AnalysisTree/TaskManager.hpp"
#include "AnalysisTree/Variable.hpp"
#include "Task.hpp"

using namespace AnalysisTree;

void example(const std::string& filelist){
  auto* man = TaskManager::GetInstance();

  auto* task = new QA::Task;
  task->SetOutputFileName("cbm_qa.root");

  // 1D histo
  task->AddH1({"p_{T}, GeV/c", Variable::FromString("VtxTracks.pT"), {100, 0, 3}});

  // 1D histo
  task->AddH1({"#Psi_{RP}", Variable::FromString("SimEventHeader.psi_RP"), {100, 0, 6.28}});

  // 1D histo with cut
  Cuts* pT_cut =  new Cuts("pT_cut", {RangeCut("VtxTracks.pT", 1, 1.5)});
  task->AddH1({"p_{T}, GeV/c", Variable::FromString("VtxTracks.pT"), {100, 0, 3}}, pT_cut);

  // AnalysisTree::Variable in case of more complicated plot
  Variable chi2_over_ndf("#chi^{2}/NDF", {{"VtxTracks", "chi2"}, {"VtxTracks", "ndf"}}, []( std::vector<double>& var ) { return var.at(0)/var.at(1); });
  task->AddH1({"#chi^{2}/NDF", chi2_over_ndf, {100, 0, 10}});

  // 2D histo
  task->AddH2({"#eta", Variable::FromString("VtxTracks.eta"), {100, -1, 4}}, {"p_{T}, GeV/c", Variable::FromString("VtxTracks.pT"), {100, 0, 3}});

  Variable qp_sts("qp_reco", {{"VtxTracks", "q"}, {"VtxTracks", "p"}}, [](std::vector<double>& qp) { return qp.at(0) * qp.at(1); });
  task->AddH2({"sign(q)*p, GeV/c", qp_sts, {500, -10, 10}},{"m^{2}, GeV^{2}/c^{2}", {"TofHits", "mass2"}, {500, -1, 2}});

  // Histo with additional cuts:
  Cuts* mc_protons =  new Cuts("McProtons", {EqualsCut("SimParticles.pid", 2212)});
  Cuts* mc_pions =  new Cuts("McPions", {EqualsCut("SimParticles.pid", 211)});
  task->AddH1({"MC-protons #chi^{2}/NDF", chi2_over_ndf, {100, 0, 10}}, mc_protons);
  task->AddH1({"MC-pions #chi^{2}/NDF", chi2_over_ndf, {100, 0, 10}}, mc_pions);

  // TProfiles
  const Field psi_RP = Field("SimEventHeader", "psi_RP");
  const Field mc_phi = Field("SimParticles", "phi");
  Variable v1("v1", {mc_phi, psi_RP}, [](std::vector<double> phi) { return cos(phi[0] - phi[1]); });
  Variable v2("v2", {mc_phi, psi_RP}, [](std::vector<double> phi) { return cos(2 * (phi[0] - phi[1])); });
  task->AddProfile({"#it{y}", Variable::FromString("SimParticles.rapidity"), {20, 0.5, 2.5}}, {"MC-protons v_{1}", v1, {}}, mc_protons);
  task->AddProfile({"#it{y}", Variable::FromString("SimParticles.rapidity"), {20, 0.5, 2.5}}, {"MC-protons v_{1}", v1, {}}, mc_pions);

  man->AddTask(task);

  man->Init({filelist}, {"rTree"});
  man->SetVerbosityPeriod(100);
  man->Run(-1);
  man->Finish();
}

int main(int argc, char* argv[]){
    if (argc < 2) {
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./example filename" << std::endl;
    exit(EXIT_FAILURE);
  }

  const std::string filename = argv[1];
  example(filename);

  return 0;
}
