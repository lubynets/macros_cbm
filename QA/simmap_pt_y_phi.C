float MidRapidityByPbeam(float pbeam);

void simmap_pt_y_phi(std::string filelist, float pbeam=12.)
{
  std::cout << "Macro started\n";
  
  const TString y_name = "y_{LAB}";
  const TString pT_name = "pT, GeV/c";
  const TString phi_name = "#varphi, rad";
  const int y_nbins = 22;
  const int pT_nbins = 26;
  const int phi_nbins = 30;
  const float midrapidity = MidRapidityByPbeam(pbeam);
  const float y_low = midrapidity - 0.8;
  const float y_up = midrapidity + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.6;
  const float phi_low = -TMath::Pi();
  const float phi_up = TMath::Pi();
  
  std::cout << midrapidity << std::endl;
  
  std::vector<int> pdgs{3122, 310};
  
  AnalysisTree::Chain* sim_tree = new AnalysisTree::Chain(std::vector<std::string>({filelist}), std::vector<std::string>({"aTree"}));
  
  auto* sim_tracks = new AnalysisTree::Particles();
  
  sim_tree -> SetBranchAddress("SimParticles", &sim_tracks);
  
  std::vector<TH2F*> hsim_y_pt, hsim_phi_pt, hsim_phi_y;
  hsim_y_pt.resize(pdgs.size());
  hsim_phi_pt.resize(pdgs.size());
  hsim_phi_y.resize(pdgs.size());
  
  for(int i=0; i<pdgs.size(); i++)
  {
    hsim_y_pt.at(i) = new TH2F("hsim_y_pt", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    hsim_phi_pt.at(i) = new TH2F("hsim_phi_pt", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    hsim_phi_y.at(i) = new TH2F("hsim_phi_y", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    
    hsim_y_pt.at(i)->GetXaxis()->SetTitle(y_name);
    hsim_phi_pt.at(i)->GetXaxis()->SetTitle(phi_name);
    hsim_phi_y.at(i)->GetXaxis()->SetTitle(phi_name);
    
    hsim_y_pt.at(i)->GetYaxis()->SetTitle(pT_name);
    hsim_phi_pt.at(i)->GetYaxis()->SetTitle(pT_name);
    hsim_phi_y.at(i)->GetYaxis()->SetTitle(y_name);    
  }
  
  const int mother_id_id = sim_tree->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("mother_id");
  
  const int sim_entries = sim_tree->GetEntries();
  for(int iEvent=0; iEvent<sim_entries; iEvent++)
  {
    sim_tree->GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << std::endl;
    for(const auto& simtrack : *(sim_tracks->GetChannels()))
    {
      const int mother_id = simtrack.GetField<int>(mother_id_id);
      if(mother_id != -1) continue;
      
      const int pid = simtrack.GetPid();
      const int index = std::find(pdgs.begin(), pdgs.end(), pid)-pdgs.begin();
      if (index>=pdgs.size()) continue;
      
      hsim_y_pt.at(index)->Fill(simtrack.GetRapidity(), simtrack.GetPt());
      hsim_phi_pt.at(index)->Fill(simtrack.GetPhi(), simtrack.GetPt());
      hsim_phi_y.at(index)->Fill(simtrack.GetPhi(), simtrack.GetRapidity());
    }    
  }
  
  TFile fileOut("simmap_pt_y_phi.root", "recreate");
  
  for(int i=0; i<pdgs.size(); i++)
  {
    fileOut.mkdir(std::to_string(pdgs.at(i)).c_str());
    fileOut.cd(std::to_string(pdgs.at(i)).c_str());
    
    hsim_y_pt.at(i)->Write();
    hsim_phi_pt.at(i)->Write();
    hsim_phi_y.at(i)->Write();    
  }
  
  fileOut.Close();
  
  std::cout << "Macro finished successfully\n";
}

float MidRapidityByPbeam(float pbeam)
{
  const float m = 0.938;
  const float E = std::sqrt(m*m + pbeam*pbeam);
  TLorentzVector v(0, 0, pbeam, E);
  
  return v.Rapidity()/2;
}
  