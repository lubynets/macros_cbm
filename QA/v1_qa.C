void v1_qa(const TString filename) {
  
  TFile* fileIn = TFile::Open(filename, "read");
  TTree* treeIn = fileIn->Get<TTree>("aTree");
  
  auto* sim_tracks = new AnalysisTree::Particles();
  auto* sim_event_header = new AnalysisTree::EventHeader();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) fileIn->Get("Configuration");
  
  treeIn->SetBranchAddress("SimParticles", &sim_tracks);
  treeIn->SetBranchAddress("SimEventHeader", &sim_event_header);
  
  const int psi_RP_id_ = config->GetBranchConfig("SimEventHeader").GetFieldId("psi_RP");
  const int mother_id_id_ = config->GetBranchConfig("SimParticles").GetFieldId("mother_id");
  
  const float y_min = 1.02179;
  const float y_max = 2.62179;
  const int y_nbins = 4;
  
  const float flow_min = -1;
  const float flow_max = 1;
  const int flow_nbins = 2000;
  
  const float pT_min = 0.2;
  const float pT_max = 1.4;
  
  TH2D* histoall = new TH2D("histoall", "histoall", y_nbins, y_min, y_max, flow_nbins, flow_min, flow_max);
  histoall->GetXaxis()->SetTitle("y_{CM}");
  histoall->GetYaxis()->SetTitle("cos(#varphi - #Psi_{RP})");  
  
  TH2D* histoprim = new TH2D("histoprim", "histoprim", y_nbins, y_min, y_max, flow_nbins, flow_min, flow_max);
  histoprim->GetXaxis()->SetTitle("y_{CM}");
  histoprim->GetYaxis()->SetTitle("cos(#varphi - #Psi_{RP})");  
  
  const int n_entries = treeIn->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event) {
    treeIn->GetEntry(i_event);
    
    const float psi_RP = sim_event_header->GetField<float>(psi_RP_id_);
    
    for(const auto& mc_track : *(sim_tracks->GetChannels()) ) {
      
      const int pid = mc_track.GetPid();
      if(pid != 3122) continue;
      
      const float pT = mc_track.GetPt();
      if(pT < pT_min || pT > pT_max) continue;
      
      const float phi = mc_track.GetPhi();
      const float y = mc_track.GetRapidity();
      const float flow = TMath::Cos(phi-psi_RP);

      histoall->Fill(y, flow);
      
      const int mother_id = mc_track.GetField<int>(mother_id_id_);
      if(mother_id != -1) continue;
      
      histoprim->Fill(y, flow);      
    }
  }
  
  TFile* fileOut = TFile::Open("v1_qa.root", "recreate");
  histoall->Write();
  histoprim->Write();
  auto* profileall = histoall->ProfileX();
  auto* profileprim = histoprim->ProfileX();
  profileall->GetYaxis()->SetTitle("#LT cos(#varphi - #Psi_{RP}) #GT");
  profileall->Write();
  profileprim->GetYaxis()->SetTitle("#LT cos(#varphi - #Psi_{RP}) #GT");
  profileprim->Write();
  fileOut->Close();    
  
}