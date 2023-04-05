#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

void m2_pq_vtx(std::string filelist) {
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
  
  auto* rec_event_header = new AnalysisTree::EventHeader();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* sim_particles = new AnalysisTree::Particles();
  auto* tof_hits = new AnalysisTree::HitDetector();
  auto* vtx_tof_matching = new AnalysisTree::Matching();
  auto* tof_sim_matching = new AnalysisTree::Matching();
  auto* vtx_sim_matching = new AnalysisTree::Matching();

  std::string dot = ".";  // brex
//   std::string dot = "";  // nobrex

  treeIn -> SetBranchAddress(("RecEventHeader" + dot).c_str(), &rec_event_header);
  treeIn -> SetBranchAddress(("VtxTracks" + dot).c_str(), &vtx_tracks);
  treeIn -> SetBranchAddress(("SimParticles" + dot).c_str(), &sim_particles);
  treeIn -> SetBranchAddress(("TofHits" + dot).c_str(), &tof_hits);
  treeIn -> SetBranchAddress((config->GetMatchName("VtxTracks", "TofHits") + dot).c_str(), &vtx_tof_matching);
  treeIn -> SetBranchAddress((config->GetMatchName("TofHits", "SimParticles") + dot).c_str(), &tof_sim_matching);
  treeIn -> SetBranchAddress((config->GetMatchName("VtxTracks", "SimParticles") + dot).c_str(), &vtx_sim_matching);

  const int vtx_chi2_id = config->GetBranchConfig("RecEventHeader").GetFieldId("vtx_chi2");
  const int mc_pdg_id = config->GetBranchConfig("VtxTracks").GetFieldId("mc_pdg");
  const int vtx_track_chi2_id = config->GetBranchConfig("VtxTracks").GetFieldId("vtx_chi2");
  const int mass2_id = config->GetBranchConfig("TofHits").GetFieldId("mass2");
  const int qp_id = config->GetBranchConfig("TofHits").GetFieldId("qp_tof");
  const int dx_id = config->GetBranchConfig("TofHits").GetFieldId("dx");
  const int dy_id = config->GetBranchConfig("TofHits").GetFieldId("dy");
    
  const int Nevents = treeIn->GetEntries();
//   const int Nevents = 1000;
  
  struct Particle {
    std::string name_;
    int pdg_;
  };
  
  std::vector<Particle> particles {
    {"h2TofM2_eminus", 11},
    {"h2TofM2_eplus", -11},
    {"h2TofM2_muminus", 13},
    {"h2TofM2_muplus", -13},
    {"h2TofM2_kminus", -321},
    {"h2TofM2_piminus", -211},
    {"h2TofM2_kplus", 321},
    {"h2TofM2_piplus", 211},
    {"h2TofM2_p", 2212},
    {"h2TofM2_d", 1000010020},
    {"h2TofM2_He3", 1000020030}
  };
  
  TH2F* h_all = new TH2F("h2TofM2", "all", p_nbins, p_low, p_up, m2_nbins, m2_low, m2_up);
  SetAxesNames(h_all, p_name, m2_name);
  std::vector<TH2F*> histo_from_tof;
  histo_from_tof.resize(particles.size());
  std::vector<TH2F*> histo_via_vtx;
  histo_via_vtx.resize(particles.size());

  int i=0;
  for(auto& pa : particles) {
    histo_from_tof.at(i) = new TH2F(pa.name_.c_str(), std::to_string(pa.pdg_).c_str(), p_nbins, p_low, p_up, m2_nbins, m2_low, m2_up);
    SetAxesNames(histo_from_tof.at(i), p_name, m2_name);
    histo_via_vtx.at(i) = new TH2F(pa.name_.c_str(), std::to_string(pa.pdg_).c_str(), p_nbins, p_low, p_up, m2_nbins, m2_low, m2_up);
    SetAxesNames(histo_via_vtx.at(i), p_name, m2_name);
    i++;
  }
  
  for(int iEvent=0; iEvent<Nevents; iEvent++) {
    treeIn -> GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << "\n";
    
    if( rec_event_header->GetField<float>(vtx_chi2_id) >=3 ) continue;

    for(const auto& tof_hit : *tof_hits) {
      const float dx = tof_hit.GetField<float>(dx_id);
      const float dy = tof_hit.GetField<float>(dy_id);
      if(dx*dx + dy*dy > 2.25) continue;

      const float mass2 = tof_hit.GetField<float>(mass2_id);
      const float p_over_q = tof_hit.GetField<float>(qp_id);

      const int vtx_track_id = vtx_tof_matching->GetMatchInverted(tof_hit.GetId());
      if(vtx_track_id < 0) continue;  // not legal situation, just for oversure

      const auto& vtx_track = vtx_tracks->GetChannel(vtx_track_id);

      if( vtx_track.GetField<float>(vtx_track_chi2_id) >=18 ) continue;

      h_all -> Fill(p_over_q, mass2);

      const int sim_track_from_tof_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      int sim_track_from_tof_pdg = -999;

      if(sim_track_from_tof_id >= 0) {
        const auto& sim_track_from_tof = sim_particles->GetChannel(sim_track_from_tof_id);
        sim_track_from_tof_pdg = sim_track_from_tof.GetPid();
      }

      const int sim_track_via_vtx_id = vtx_sim_matching->GetMatch(vtx_track.GetId());
      int sim_track_via_vtx_pdg = -999;

      if(sim_track_via_vtx_id >= 0) {
        const auto& sim_track_via_vtx = sim_particles->GetChannel(sim_track_via_vtx_id);
        sim_track_via_vtx_pdg = sim_track_via_vtx.GetPid();
      }

      int j=0;
      for(auto& pa : particles) {
        if(sim_track_from_tof_pdg == pa.pdg_) {
          histo_from_tof.at(j) -> Fill(p_over_q, mass2);
        }
        if(sim_track_via_vtx_pdg == pa.pdg_) {
          histo_via_vtx.at(j) -> Fill(p_over_q, mass2);
        }
        j++;
      }
    }
  }
  
  TFile fileOut("m2_pq_vtx.root", "recreate");
  fileOut.mkdir("reco_info");
  fileOut.cd("reco_info");
  h_all -> Write();

  fileOut.mkdir("reco_vs_sim_info_from_tof");
  fileOut.cd("reco_vs_sim_info_from_tof");
  for(auto& h : histo_from_tof) {
    h -> Write();
  }

  fileOut.mkdir("reco_vs_sim_info_via_vtx");
  fileOut.cd("reco_vs_sim_info_via_vtx");
  for(auto& h : histo_via_vtx) {
    h -> Write();
  }

  fileOut.Close();
    
  std::cout << "Macro finished successfully\n";
}
