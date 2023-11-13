void tofhits_check(std::string filelist) {

  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist}), std::vector<std::string>({"aTree"}));
  auto* config = treeIn->GetConfiguration();

  auto* vtx_tracks = new AnalysisTree::Particles();
  auto* sim_particles = new AnalysisTree::Particles();
  auto* tof_hits = new AnalysisTree::HitDetector();
  auto* vtx_tof_matching = new AnalysisTree::Matching();
  auto* tof_sim_matching = new AnalysisTree::Matching();
  auto* vtx_sim_matching = new AnalysisTree::Matching();

  const int Nevents = treeIn->GetEntries();

  std::string dot = ".";  // brex
//   std::string dot = "";  // nobrex

  treeIn -> SetBranchAddress(("RecParticles" + dot).c_str(), &vtx_tracks);
  treeIn -> SetBranchAddress(("SimParticles" + dot).c_str(), &sim_particles);
  treeIn -> SetBranchAddress(("TofHits" + dot).c_str(), &tof_hits);
  treeIn -> SetBranchAddress((config->GetMatchName("RecParticles", "TofHits") + dot).c_str(), &vtx_tof_matching);
  treeIn -> SetBranchAddress((config->GetMatchName("TofHits", "SimParticles") + dot).c_str(), &tof_sim_matching);
  treeIn -> SetBranchAddress((config->GetMatchName("RecParticles", "SimParticles") + dot).c_str(), &vtx_sim_matching);

  int Nall{0};
  int Ngood{0};

  for(int iEvent=0; iEvent<Nevents; iEvent++) {
    treeIn -> GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << "\n";

    for(const auto& tof_hit : *tof_hits) {
      Nall++;

      const int direct_sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());


      const int vtx_track_id = vtx_tof_matching->GetMatchInverted(tof_hit.GetId());
      if(vtx_track_id < 0) continue;  // not legal situation, just for oversure

      const int indirect_sim_id = vtx_sim_matching->GetMatch(vtx_track_id);

      if(direct_sim_id == indirect_sim_id) Ngood++;


    }

  }

  std::cout << Nall << "\n";
  std::cout << Ngood << "\n";

}
