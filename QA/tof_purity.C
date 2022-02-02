void tof_purity(std::string filelist, std::string filemapmapname)
{
  
  std::cout << "Macro started\n";
  
  TFile* fileMap = TFile::Open(filemapmapname.c_str());
  
  const TString xaxisname = "purity";
  const int nbins = 100;
  const float low = 0;
  const float up = 1;
  const TString yaxisname = "Entries";
  
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
    {"pi_plus", 211},
    {"pi_minus", -211},
    {"K_plus", 321},
    {"K_minus", -321},
    {"proton_plus", 2212},
    {"deuteron", 1000010020},
  };
  
  std::vector<TH2F*> histo_purmap;
  std::vector<TH2F*> histo_errmap;
  std::vector<TH1F*> histo_pur;
  
  histo_purmap.resize(particles.size());
  histo_errmap.resize(particles.size());
  histo_pur.resize(particles.size());
  
  int i = 0;
  for(auto& pa : particles)
  {
    histo_purmap.at(i) = fileMap->Get<TH2F>((pa.name_ + "/pur").c_str());
    histo_errmap.at(i) = fileMap->Get<TH2F>((pa.name_ + "/err").c_str());
    histo_pur.at(i) = new TH1F(("h_" + pa.name_).c_str(), (std::to_string(pa.pdg_)).c_str(), nbins, low, up);
    histo_pur.at(i) -> GetXaxis() -> SetTitle(xaxisname);
    histo_pur.at(i) -> GetYaxis() -> SetTitle(yaxisname);
    i++;
  }
  
  for(int iEvent=0; iEvent<Nevents; iEvent++)
  {
    treeIn -> GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << "\n";
    
    for(const auto& tof_hit : *tof_hits)
    {
      const int sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      if(sim_id<0) continue;
      
      auto& sim_track = sim_tracks->GetChannel(sim_id);
      const int sim_pdg = sim_track.GetPid();
      
      int j=0;
      for(auto& pa : particles)
      {
        if(sim_pdg == pa.pdg_)
        {
          const float mass2 = tof_hit.GetField<float>(mass2_id);
          const float pq = tof_hit.GetField<float>(qp_id);
          const float pur = histo_purmap.at(j)->GetBinContent(histo_purmap.at(j)->FindBin(pq, mass2));
          const float err = histo_errmap.at(j)->GetBinContent(histo_errmap.at(j)->FindBin(pq, mass2));
//           if(err < 0.05)
            histo_pur.at(j) -> Fill(pur);
          break;
        }
        j++;
      }
    }    
  }    
  
  TFile fileOut("tof_purity.root", "recreate");
  for(auto& h : histo_pur)
    h -> Write();
  
  fileOut.Close();
    
  std::cout << "Macro finished successfully\n";  
}