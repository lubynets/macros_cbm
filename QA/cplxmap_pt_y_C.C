#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

void cplxmap_pt_y_C(std::string filelist_sim, std::string filelist_rec, float pbeam=12.) {
  
  std::cout << "Macro started\n";
  
  const TString y_name = "y_{LAB}";
  const TString pT_name = "p_{T}, GeV/c";
  const TString C_name = "Centrality, %";
  const int y_nbins = 22;
  const int pT_nbins = 26;
//   const int y_nbins = 11;
//   const int pT_nbins = 13;
  const float midrapidity = MidRapidityByPbeam(pbeam);
  const float y_low = midrapidity - 0.8;
  const float y_up = midrapidity + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.6;
  std::vector<float> C_edges{0, 5, 10, 20, 30, 40, 70, 100};
  const int C_nbins = C_edges.size() - 1;

  std::cout << midrapidity << std::endl;
  
  std::vector<int> pdgs{3122, 310}; std::string recotree = "pTree";
//   std::vector<int> pdgs{3312, 3334}; std::string recotree = "aTree";
  
  float* y_bins = new float[y_nbins+1];
  for(int i=0; i<y_nbins+1; i++)
    y_bins[i] = y_low + i*(y_up - y_low) / y_nbins;

  float* pT_bins = new float[pT_nbins+1];
  for(int i=0; i<pT_nbins+1; i++)
    pT_bins[i] = pT_low + i*(pT_up - pT_low) / pT_nbins;

  float* C_bins = new float[C_nbins+1];
  for(int i=0; i<C_nbins+1; i++)
    C_bins[i] = C_edges.at(i);
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_sim, filelist_rec}), std::vector<std::string>({"aTree", recotree.c_str()}));
  auto* config = treeIn->GetConfiguration();

  auto* eve_header=  new AnalysisTree::EventHeader();
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* reco_tracks = new AnalysisTree::Particles();
  auto* simulated_copied = new AnalysisTree::Particles();
  auto* reco_sim_matching = new AnalysisTree::Matching();
  
  treeIn -> SetBranchAddress("RecEventHeader.", &eve_header);
  treeIn -> SetBranchAddress("SimParticles.", &sim_tracks);
  treeIn -> SetBranchAddress("Candidates.", &reco_tracks);
  treeIn -> SetBranchAddress("Simulated.", &simulated_copied);
  treeIn -> SetBranchAddress((config->GetMatchName("Candidates", "Simulated") + ".").c_str(), &reco_sim_matching);

  std::vector<TH3F*> hsim_y_pt_C;
  hsim_y_pt_C.resize(pdgs.size());
  
  std::vector<TH3F*> hrec_y_pt_C;
  hrec_y_pt_C.resize(pdgs.size());
  
  for(int i=0; i<pdgs.size(); i++) {
    std::string histoname = std::to_string(pdgs.at(i));
    hsim_y_pt_C.at(i) = new TH3F("hsim_y_pt_C", histoname.c_str(), y_nbins, y_bins, pT_nbins, pT_bins, C_nbins, C_bins);
    SetAxesNames(hsim_y_pt_C.at(i), y_name, pT_name, C_name);

    hrec_y_pt_C.at(i) = new TH3F("hrec_y_pt_C", histoname.c_str(), y_nbins, y_bins, pT_nbins, pT_bins, C_nbins, C_bins);
    SetAxesNames(hrec_y_pt_C.at(i), y_name, pT_name, C_name);
  }
  
  const int centrality_id = treeIn->GetConfiguration()->GetBranchConfig("RecEventHeader").GetFieldId("centrality_tracks");
  const int mother_id_id = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("mother_id");
  const int g4_process_id_id = treeIn->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("geant_process_id");
  const int generation_id = treeIn->GetConfiguration()->GetBranchConfig("Candidates").GetFieldId("generation");
  
  const int n_entries = treeIn->GetEntries();
//   const int n_entries = 1000;
  
  for(int iEvent=0; iEvent<n_entries; iEvent++) {
    treeIn->GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << std::endl;
    
    const float centrality = eve_header->GetField<float>(centrality_id);
//     const int  C_bin = FindBin(C_edges, centrality);
    
    for(const auto& simtrack : *(sim_tracks->GetChannels())) {
      
      const int mother_id = simtrack.GetField<int>(mother_id_id);
      if(mother_id != -1) continue;
      
//       const int g4_process_id = simtrack.GetField<int>(g4_process_id_id);
//       if (g4_process_id>10) continue;
      
      const int pid = simtrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;
      
      hsim_y_pt_C.at(index)->Fill(simtrack.GetRapidity(), simtrack.GetPt(), centrality);
    }
    
    for(const auto& recotrack : *(reco_tracks->GetChannels())) {
      const int generation = recotrack.GetField<int>(generation_id);
      if(generation != 1) continue;
      const int pid = recotrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;

      auto matched_sim_track_id = reco_sim_matching->GetMatch(recotrack.GetId());
      auto& matched_sim_track = simulated_copied->GetChannel(matched_sim_track_id);
      
      hrec_y_pt_C.at(index)->Fill(matched_sim_track.GetRapidity(), matched_sim_track.GetPt(), centrality);
    }      
  }
  
  TFile fileOut("cplxmap_pt_y_C.root", "recreate");
  
  for(int i=0; i<pdgs.size(); i++) {
    fileOut.mkdir(std::to_string(pdgs.at(i)).c_str());
    fileOut.cd(std::to_string(pdgs.at(i)).c_str());

    hsim_y_pt_C.at(i)->Write();
    hrec_y_pt_C.at(i)->Write();
  }
  
  fileOut.Close();
  
  std::cout << "Macro finished successfully\n";  
}
