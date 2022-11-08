#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

int FindBin(std::vector<float>& v, float value);

void cplxmap_pt_y_C(std::string filelist_sim, std::string filelist_rec, float pbeam=12.) {
  
  std::cout << "Macro started\n";
  
  const TString y_name = "y_{LAB}";
  const TString pT_name = "p_{T}, GeV/c";
  const int y_nbins = 22;
  const int pT_nbins = 26;
//   const int y_nbins = 11;
//   const int pT_nbins = 13;
  const float midrapidity = MidRapidityByPbeam(pbeam);
  const float y_low = midrapidity - 0.8;
  const float y_up = midrapidity + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.6;
  
  std::cout << midrapidity << std::endl;
  
  std::vector<int> pdgs{3122, 310}; std::string recotree = "pTree";
//   std::vector<int> pdgs{3312, 3334}; std::string recotree = "aTree";
  
  std::vector<float> C_edges{0, 10, 20, 40, 70, 100};
  const int C_nbins = C_edges.size() - 1;
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_sim, filelist_rec}), std::vector<std::string>({"aTree", recotree.c_str()}));
  
  auto* eve_header=  new AnalysisTree::EventHeader();
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* reco_tracks = new AnalysisTree::Particles();
  
  treeIn -> SetBranchAddress("AnaEventHeader", &eve_header);
  treeIn -> SetBranchAddress("SimParticles", &sim_tracks);
  treeIn -> SetBranchAddress("Candidates", &reco_tracks);
  
  std::vector<std::vector<TH2F*>> hsim_y_pt;
  hsim_y_pt.resize(pdgs.size());
  for(auto& h : hsim_y_pt)
    h.resize(C_nbins);
  
  std::vector<std::vector<TH2F*>> hrec_y_pt;
  hrec_y_pt.resize(pdgs.size());
  for(auto& h : hrec_y_pt)
    h.resize(C_nbins);
  
  for(int i=0; i<pdgs.size(); i++) {
    for(int j=0; j< C_nbins; j++) {
      std::string histoname = std::to_string(pdgs.at(i)) + "_" + std::to_string(j);
      hsim_y_pt.at(i).at(j) = new TH2F("hsim_y_pt", histoname.c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
      SetAxesNames(hsim_y_pt.at(i).at(j), y_name, pT_name);
      
      hrec_y_pt.at(i).at(j) = new TH2F("hrec_y_pt", histoname.c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
      SetAxesNames(hrec_y_pt.at(i).at(j), y_name, pT_name);
    }
  }
  
  const int centrality_id = treeIn->GetConfiguration()->GetBranchConfig("AnaEventHeader").GetFieldId("centrality_tracks");
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
    
    const float centrality = eve_header->GetField<float>(centrality_id);
    const int  C_bin = FindBin(C_edges, centrality);
    
    for(const auto& simtrack : *(sim_tracks->GetChannels())) {
      
      const int mother_id = simtrack.GetField<int>(mother_id_id);
      if(mother_id != -1) continue;
      
//       const int g4_process_id = simtrack.GetField<int>(g4_process_id_id);
//       if (g4_process_id>10) continue;
      
      const int pid = simtrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;
      
      hsim_y_pt.at(index).at(C_bin)->Fill(simtrack.GetRapidity(), simtrack.GetPt());
    }
    
    for(const auto& recotrack : *(reco_tracks->GetChannels())) {
      const int generation = recotrack.GetField<int>(generation_id);
      if(generation != 1) continue;
      const int pid = recotrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;
      
      hrec_y_pt.at(index).at(C_bin)->Fill(recotrack.GetRapidity(), recotrack.GetPt());
    }      
  }
  
  TFile fileOut("cplxmap_pt_y_C.root", "recreate");
  
  for(int i=0; i<pdgs.size(); i++) {
    for(int j=0; j< C_nbins; j++) {
      fileOut.mkdir((std::to_string(pdgs.at(i)) + "/" + std::to_string(j)).c_str());
      fileOut.cd((std::to_string(pdgs.at(i)) + "/" + std::to_string(j)).c_str());
      
      hsim_y_pt.at(i).at(j)->Write();
      hrec_y_pt.at(i).at(j)->Write();
    }
  }
  
  fileOut.Close();
  
  std::cout << "Macro finished successfully\n";  
}

int FindBin(std::vector<float>& v, float value) {
  if(!std::is_sorted(v.begin(), v.end()))
    throw std::runtime_error("Vector of centrality is not ordered!");
  
  int bin = -1;
  for(auto& ele : v) {
    if(value>ele)
      bin++;
    else
      break;
  }
  
  return bin;
}