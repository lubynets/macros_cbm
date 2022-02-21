float MidRapidityByPbeam(float pbeam);

void recmap_pt_y_phi(std::string filelist, float pbeam=12.)
{
  std::cout << "Macro started\n";
  
  const TString y_name = "y_{LAB}";
  const TString pT_name = "pT, GeV/c";
  const TString phi_name = "#varphi, rad";
//   const int y_nbins = 22;
//   const int pT_nbins = 26;
//   const int phi_nbins = 30;
  const int y_nbins = 11;
  const int pT_nbins = 13;
  const int phi_nbins = 15;
  const float midrapidity = MidRapidityByPbeam(pbeam);
  const float y_low = midrapidity - 0.8;
  const float y_up = midrapidity + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.6;
  const float phi_low = -TMath::Pi();
  const float phi_up = TMath::Pi();
  
//   std::vector<int> pdgs{3312};
  std::vector<int> pdgs{3334};
//   std::vector<int> pdgs{3122, 310};

  AnalysisTree::Chain* reco_tree = new AnalysisTree::Chain(std::vector<std::string>({filelist}), std::vector<std::string>({"aTree"}));
  
  auto* reco_tracks = new AnalysisTree::Particles();
  
  reco_tree -> SetBranchAddress("Candidates", &reco_tracks);
  
  std::vector<TH2F*> hrec_y_pt, hrec_phi_pt, hrec_phi_y;
  hrec_y_pt.resize(pdgs.size());
  hrec_phi_pt.resize(pdgs.size());
  hrec_phi_y.resize(pdgs.size());
  
  for(int i=0; i<pdgs.size(); i++)
  {
    hrec_y_pt.at(i) = new TH2F("hrec_y_pt", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    hrec_phi_pt.at(i) = new TH2F("hrec_phi_pt", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    hrec_phi_y.at(i) = new TH2F("hrec_phi_y", std::to_string(pdgs.at(i)).c_str(), y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
    
    hrec_y_pt.at(i)->GetXaxis()->SetTitle(y_name);
    hrec_phi_pt.at(i)->GetXaxis()->SetTitle(phi_name);
    hrec_phi_y.at(i)->GetXaxis()->SetTitle(phi_name);
    
    hrec_y_pt.at(i)->GetYaxis()->SetTitle(pT_name);
    hrec_phi_pt.at(i)->GetYaxis()->SetTitle(pT_name);
    hrec_phi_y.at(i)->GetYaxis()->SetTitle(y_name);    
  }
  
  const int generation_id = reco_tree->GetConfiguration()->GetBranchConfig("Candidates").GetFieldId("generation");
  
  const int rec_entries = reco_tree->GetEntries();
  for(int iEvent=0; iEvent<rec_entries; iEvent++)
  {
    reco_tree->GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << std::endl;
    for(const auto& recotrack : *(reco_tracks->GetChannels()))
    {
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
  
  TFile fileOut("recmap_pt_y_phi.root", "recreate");
  
  for(int i=0; i<pdgs.size(); i++)
  {
    fileOut.mkdir(std::to_string(pdgs.at(i)).c_str());
    fileOut.cd(std::to_string(pdgs.at(i)).c_str());
    
    hrec_y_pt.at(i)->Write();
    hrec_phi_pt.at(i)->Write();
    hrec_phi_y.at(i)->Write();    
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
  