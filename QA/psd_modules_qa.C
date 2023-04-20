#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

void psd_modules_qa(std::string filelist) {
  std::cout << "Macro started\n";

  const int nbinsx = 8;
  const int nbinsy = 6;
  const float length = 20;

  TH2F* h = new TH2F("h", "", nbinsx, -length*nbinsx/2, length*nbinsx/2, nbinsy, -length*nbinsy/2, length*nbinsy/2);

  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist}), std::vector<std::string>({"rTree"}));
  auto* config = treeIn->GetConfiguration();
  auto* dataheader = treeIn->GetDataHeader();

  auto* psd_modules = new AnalysisTree::ModuleDetector();

  std::string dot = ".";  // brex
//   std::string dot = "";  // nobrex

  treeIn -> SetBranchAddress(("PsdModules" + dot).c_str(), &psd_modules);

  const int Nevents = treeIn->GetEntries();

  std::vector<float> X, Y;

  for(auto& ch : dataheader->GetModulePositions(0)) {
    X.push_back(ch.GetPosition().X());
    Y.push_back(ch.GetPosition().Y());
  }

  for(int iEvent=0; iEvent<Nevents; iEvent++) {
    treeIn -> GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << "\n";

    int ich{0};
    for(auto& psd_module : *psd_modules) {
      const float module_signal = psd_module.GetSignal();
      h->Fill(X.at(ich), Y.at(ich), module_signal);
      ich++;
    }
  }

  TFile fileOut("psd_modules_qa.root", "recreate");
  h->Write();
  fileOut.Close();

  std::cout << "Macro finished successfully\n";
}
