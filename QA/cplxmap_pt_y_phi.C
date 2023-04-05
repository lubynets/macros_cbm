#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

void cplxmap_pt_y_phi(std::string filelist_sim, std::string filelist_rec, float pbeam=12.) {
  
  std::cout << "Macro started\n";
  
  const TString y_name = "y_{LAB}";
  const TString pT_name = "p_{T}, GeV/c";
  const TString phi_name = "#varphi, rad";
  const int y_nbins = 22;
  const int pT_nbins = 26;
  const int phi_nbins = 30;
//   const int y_nbins = 11;
//   const int pT_nbins = 13;
//   const int phi_nbins = 15;
  const float midrapidity = MidRapidityByPbeam(pbeam);
  const float y_low = midrapidity - 0.8;
  const float y_up = midrapidity + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.6;
  const float phi_low = -TMath::Pi();
  const float phi_up = TMath::Pi();
  
  std::cout << midrapidity << std::endl;
  
  std::vector<int> pdgs{3122, 310}; std::string recotree = "pTree";
//   std::vector<int> pdgs{3312, 3334}; std::string recotree = "aTree";
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_sim, filelist_rec}), std::vector<std::string>({"rTree", recotree.c_str()}));
  
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* reco_tracks = new AnalysisTree::Particles();
  
  treeIn -> SetBranchAddress("SimParticles", &sim_tracks);
  treeIn -> SetBranchAddress("Candidates", &reco_tracks);
  
  std::vector<TH2F*> hsim_y_pt, hsim_phi_pt, hsim_phi_y;
  hsim_y_pt.resize(pdgs.size());
  hsim_phi_pt.resize(pdgs.size());
  hsim_phi_y.resize(pdgs.size());
  
  std::vector<TH2F*> hrec_y_pt, hrec_phi_pt, hrec_phi_y;
  hrec_y_pt.resize(pdgs.size());
  hrec_phi_pt.resize(pdgs.size());
  hrec_phi_y.resize(pdgs.size());

  for(int i=0; i<pdgs.size(); i++) {
    hsim_y_pt.at(i) = new TH2F("hsim_y_pt", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    hsim_phi_pt.at(i) = new TH2F("hsim_phi_pt", std::to_string(pdgs.at(i)).c_str(), phi_nbins, phi_low, phi_up, pT_nbins, pT_low, pT_up);
    hsim_phi_y.at(i) = new TH2F("hsim_phi_y", std::to_string(pdgs.at(i)).c_str(), phi_nbins, phi_low, phi_up, y_nbins, y_low, y_up);
    
    SetAxesNames(hsim_y_pt.at(i), y_name, pT_name);
    SetAxesNames(hsim_phi_pt.at(i), phi_name, pT_name);
    SetAxesNames(hsim_phi_y.at(i), phi_name, y_name);
    
    hrec_y_pt.at(i) = new TH2F("hrec_y_pt", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    hrec_phi_pt.at(i) = new TH2F("hrec_phi_pt", std::to_string(pdgs.at(i)).c_str(), phi_nbins, phi_low, phi_up, pT_nbins, pT_low, pT_up);
    hrec_phi_y.at(i) = new TH2F("hrec_phi_y", std::to_string(pdgs.at(i)).c_str(), phi_nbins, phi_low, phi_up, y_nbins, y_low, y_up);
    
    SetAxesNames(hrec_y_pt.at(i), y_name, pT_name);
    SetAxesNames(hrec_phi_pt.at(i), phi_name, pT_name);
    SetAxesNames(hrec_phi_y.at(i), phi_name, y_name);
  }
  
  const int mother_id_id = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("mother_id");
  const int z_id = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("z");
  const int g4_process_id_id = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("geant_process_id");
  const int generation_id = treeIn->GetConfiguration()->GetBranchConfig("Candidates").GetFieldId("generation");
  
  const int n_entries = treeIn->GetEntries();
//   const int n_entries = 1000;
  
  for(int iEvent=0; iEvent<n_entries; iEvent++) {
    treeIn->GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << std::endl;
    for(const auto& simtrack : *(sim_tracks->GetChannels())) {
      
      const int mother_id = simtrack.GetField<int>(mother_id_id);
      if(mother_id != -1) continue;
      
//       const int g4_process_id = simtrack.GetField<int>(g4_process_id_id);
//       if (g4_process_id>10) continue;
      
      const int pid = simtrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;
      
      hsim_y_pt.at(index)->Fill(simtrack.GetRapidity(), simtrack.GetPt());
      hsim_phi_pt.at(index)->Fill(simtrack.GetPhi(), simtrack.GetPt());
      hsim_phi_y.at(index)->Fill(simtrack.GetPhi(), simtrack.GetRapidity());
    }
    
    for(const auto& recotrack : *(reco_tracks->GetChannels())) {
      const int generation = recotrack.GetField<int>(generation_id);
      if(generation != 1) continue;
      const int pid = recotrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;
      
      hrec_y_pt.at(index)->Fill(recotrack.GetRapidity(), recotrack.GetPt());
      hrec_phi_pt.at(index)->Fill(recotrack.GetPhi(), recotrack.GetPt());
      hrec_phi_y.at(index)->Fill(recotrack.GetPhi(), recotrack.GetRapidity());
    }      
  }
  
  TFile fileOut("cplxmap_pt_y_phi.root", "recreate");
  
  for(int i=0; i<pdgs.size(); i++) {
    fileOut.mkdir(std::to_string(pdgs.at(i)).c_str());
    fileOut.cd(std::to_string(pdgs.at(i)).c_str());
    
    hsim_y_pt.at(i)->Write();
    hsim_phi_pt.at(i)->Write();
    hsim_phi_y.at(i)->Write();
    
    hrec_y_pt.at(i)->Write();
    hrec_phi_pt.at(i)->Write();
    hrec_phi_y.at(i)->Write();
  }
  
  fileOut.Close();
  
  std::cout << "Macro finished successfully\n";  
}