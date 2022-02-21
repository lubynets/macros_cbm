void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname);

void tof_purity_compare(std::string filelist_pfs, std::string filelist_std, std::string filelist_dec, int pdg_decay=3122)
{
  std::cout << "Macro started\n";
  
  const TString x_name = "purity standard";
  const int x_nbins = 1000;
  const float x_low = 0;
  const float x_up = 1;
  const TString y_name = "purity decay";
  const int y_nbins = 1000;
  const float y_low = 0;
  const float y_up = 1;
  
  std::string name_second_daughter;
  if(pdg_decay == 3122)
    name_second_daughter = "p";
  if(pdg_decay == 310)
    name_second_daughter = "pi";
    
  AnalysisTree::Chain* tree_pfs = new AnalysisTree::Chain(std::vector<std::string>({filelist_pfs}), std::vector<std::string>({"pTree"}));
  AnalysisTree::Chain* tree_std = new AnalysisTree::Chain(std::vector<std::string>({filelist_std}), std::vector<std::string>({"aTree"}));
  AnalysisTree::Chain* tree_dec = new AnalysisTree::Chain(std::vector<std::string>({filelist_dec}), std::vector<std::string>({"aTree"}));
  
  auto* config_pfs = tree_pfs->GetConfiguration();
  auto* config_std = tree_std->GetConfiguration();
  
  auto* candidates = new AnalysisTree::Particles();
  auto* vtx_tracks_std = new AnalysisTree::Particles();
  auto* vtx_tracks_dec = new AnalysisTree::Particles();
  
  tree_pfs -> SetBranchAddress("Candidates", &candidates);
  tree_std -> SetBranchAddress("RecParticles.", &vtx_tracks_std);
  tree_dec -> SetBranchAddress("RecParticles.", &vtx_tracks_dec);
  
  const int daughter1_id_id = config_pfs->GetBranchConfig("Candidates").GetFieldId("daughter1_id");
  const int daughter2_id_id = config_pfs->GetBranchConfig("Candidates").GetFieldId("daughter2_id");
  const std::vector<int> prob_id = {config_std->GetBranchConfig("RecParticles").GetFieldId("prob_pi"),
                                    config_std->GetBranchConfig("RecParticles").GetFieldId(("prob_" + name_second_daughter).c_str())};
  
  const int Nevents = tree_pfs->GetEntries();
  
  std::vector<TH2F*> h_pur;
  h_pur.resize(2);
  for(int i=0; i<2; i++) {
    h_pur.at(i) = new TH2F(("h"+std::to_string(i+1)).c_str(), ("H"+std::to_string(i+1)).c_str(), x_nbins, x_low, x_up, y_nbins, y_low, y_up);
    SetAxesNames(h_pur.at(i), x_name, y_name);
  }
  
  for(int iEvent=0; iEvent<Nevents; iEvent++) {
    tree_pfs -> GetEntry(iEvent);
    tree_std -> GetEntry(iEvent);
    tree_dec -> GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << "\n";
    
    for(const auto& cand : *candidates) {
      
      if(cand.GetPid() != pdg_decay) continue;
      
      int i{0};
      for(auto& daughter_id_id : {daughter1_id_id, daughter2_id_id}) {
        
        const int daughter_id = cand.GetField<int>(daughter_id_id);
        
        auto& vtx_track_std = vtx_tracks_std->GetChannel(daughter_id);
        auto& vtx_track_dec = vtx_tracks_dec->GetChannel(daughter_id);
        
        const float pur_std = vtx_track_std.GetField<float>(prob_id.at(i));
        const float pur_dec = vtx_track_dec.GetField<float>(prob_id.at(i));
        
        h_pur.at(i)->Fill(pur_std, pur_dec);
        i++;
      }
    }
  }    
  
  TFile fileOut("tof_purity_compare.root", "recreate");
  h_pur.at(0)->Write();
  h_pur.at(1)->Write();  
  fileOut.Close();  
  
  std::cout << "Macro finished successfully\n";
}

void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname)
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
}