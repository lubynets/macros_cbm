void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname);

void m2_pq(std::string filelist)
{
  std::cout << "Macro started\n";
  
  const TString p_name = "p/q, GeV/c";
  const int p_nbins = 2000;
  const float p_low = -25;
  const float p_up = 25;
  const TString m2_name = "m^{2}/q^{2}, (GeV/c^{2})^{2}";
  const int m2_nbins = 2000;
  const float m2_low = -6;
  const float m2_up = 6;
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist}), std::vector<std::string>({"rTree"}));
  auto* config = treeIn->GetConfiguration();
  
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* tof_hits = new AnalysisTree::HitDetector();
  auto* tof_sim_matching = new AnalysisTree::Matching();
  
  treeIn -> SetBranchAddress("SimParticles", &sim_tracks);
  treeIn -> SetBranchAddress("TofHits", &tof_hits);
  treeIn -> SetBranchAddress(config->GetMatchName("TofHits", "SimParticles").c_str(), &tof_sim_matching);
  
  const int mass2_id = config->GetBranchConfig("TofHits").GetFieldId("mass2");
  const int qp_id = config->GetBranchConfig("TofHits").GetFieldId("qp_tof");
  
  const int Nevents = treeIn->GetEntries();
  
  struct Particle
  {
    std::string name_;
    int pdg_;
  };
  
  std::vector<Particle> particles
  {
    {"h2TofM2_kminus", -321},
    {"h2TofM2_piminus", -211},
    {"h2TofM2_kplus", 321},
    {"h2TofM2_piplus", 211},
    {"h2TofM2_p", 2212},
    {"h2TofM2_d", 1000010020},
    {"h2TofM2_He3", 1000020030},
  };
  
  TH2F* h_all = new TH2F("h2TofM2", "all", p_nbins, p_low, p_up, m2_nbins, m2_low, m2_up);
  SetAxesNames(h_all, p_name, m2_name);
  std::vector<TH2F*> histo;
  histo.resize(particles.size());
  
  int i=0;
  for(auto& pa : particles)
  {
    histo.at(i) = new TH2F(pa.name_.c_str(), std::to_string(pa.pdg_).c_str(), p_nbins, p_low, p_up, m2_nbins, m2_low, m2_up);
    SetAxesNames(histo.at(i), p_name, m2_name);
    i++;
  }
  
  for(int iEvent=0; iEvent<Nevents; iEvent++)
  {
    treeIn -> GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << "\n";
    
    for(const auto& tof_hit : *tof_hits)
    {
      const float mass2 = tof_hit.GetField<float>(mass2_id);
      const float pq = tof_hit.GetField<float>(qp_id);
      
      h_all -> Fill(pq, mass2);
      
      const int sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      if(sim_id<0) continue;
      
      auto& sim_track = sim_tracks->GetChannel(sim_id);
      const int sim_pdg = sim_track.GetPid();
      
      int j=0;
      for(auto& pa : particles)
      {
        if(sim_pdg == pa.pdg_)
        {
          histo.at(j) -> Fill(pq, mass2);
          break;
        }
        j++;
      }
    }    
  }  
  
  TFile fileOut("m2_pq.root", "recreate");
  fileOut.mkdir("reco_info");
  fileOut.cd("reco_info");
  h_all -> Write();
  fileOut.mkdir("reco_vs_sim_info");
  fileOut.cd("reco_vs_sim_info");
  for(auto& h : histo)
    h -> Write();
  
  fileOut.Close();
    
  std::cout << "Macro finished successfully\n";
}

void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname)
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
}