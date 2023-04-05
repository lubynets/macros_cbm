#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

void tof_purity_sgnl_bckgr(std::string filelist_cbm, std::string filelist_pfs, int pdg_decay=3122)
{
  std::cout << "Macro started\n";
  
  const TString x_name = "purity";
  const int x_nbins = 2100;
  const float x_low = -1.05;
  const float x_up = 1.05;
  const TString y_name = "Entries";
  
  std::string name_second_daughter;
  if(pdg_decay == 3122)
    name_second_daughter = "p";
  if(pdg_decay == 310)
    name_second_daughter = "pi";  
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_cbm, filelist_pfs}), std::vector<std::string>({"aTree", "pTree"}));
  auto* config = treeIn->GetConfiguration();
  
  auto* candidates = new AnalysisTree::Particles();
  auto* vtx_tracks = new AnalysisTree::Particles();

  treeIn -> SetBranchAddress("Candidates", &candidates);
  treeIn -> SetBranchAddress("RecParticles.", &vtx_tracks);
  
  const int generation_id = config->GetBranchConfig("Candidates").GetFieldId("generation");
  const int daughter1_id_id = config->GetBranchConfig("Candidates").GetFieldId("daughter1_id");
  const int daughter2_id_id = config->GetBranchConfig("Candidates").GetFieldId("daughter2_id");
  const std::vector<int> prob_id = {config->GetBranchConfig("RecParticles").GetFieldId("prob_pi"),
                                    config->GetBranchConfig("RecParticles").GetFieldId(("prob_" + name_second_daughter).c_str())};
                                    
  
  const int Nevents = treeIn->GetEntries();
  
  std::vector<TH1F*> hSgnl;
  std::vector<TH1F*> hBckgr;
  hSgnl.resize(2);
  hBckgr.resize(2);
  
  for(int i=0; i<2; i++) {
    hSgnl.at(i) = new TH1F(("hSgnl_" + std::to_string(i+1)).c_str(), ("HSgnl_" + std::to_string(i+1)).c_str(), x_nbins, x_low, x_up);
    SetAxesNames(hSgnl.at(i), x_name, y_name);
    hBckgr.at(i) = new TH1F(("hBckgr_" + std::to_string(i+1)).c_str(), ("HBckgr_" + std::to_string(i+1)).c_str(), x_nbins, x_low, x_up);
    SetAxesNames(hBckgr.at(i), x_name, y_name);
  }
  
  for(int iEvent=0; iEvent<Nevents; iEvent++) {  
    treeIn->GetEntry(iEvent);
    
    if(iEvent%100==0)
      std::cout << iEvent << "\n";

    for(const auto& cand : *candidates) {
      
      if(cand.GetPid() != pdg_decay) continue;
      
      const int is_signal = cand.GetField<int>(generation_id);
      
      int i{0};
      for(auto& daughter_id_id : {daughter1_id_id, daughter2_id_id}) {
        
        const int daughter_id = cand.GetField<int>(daughter_id_id);
        
        auto& vtx_track = vtx_tracks->GetChannel(daughter_id);
        
        const float purity = vtx_track.GetField<float>(prob_id.at(i));
        
        if(is_signal == 0)
          hBckgr.at(i)->Fill(purity);
        else if(is_signal>0)
          hSgnl.at(i)->Fill(purity);
        
        i++;
      }
    }
  }
  
  TFile fileOut("tof_purity_sgnl_bckgr.root", "recreate");
  fileOut.mkdir("Signal");
  fileOut.cd("Signal");
  hSgnl.at(0)->Write();
  hSgnl.at(1)->Write();
  fileOut.mkdir("Bckgr");
  fileOut.cd("Bckgr");
  hBckgr.at(0)->Write();
  hBckgr.at(1)->Write();  
  fileOut.Close();  
  
  std::cout << "Macro finished successfully\n";  
  
}