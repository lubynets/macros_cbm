void recmap_pt_y_phi(const TString recofilename)
{
  std::cout << "Macro started\n";
  
  const TString y_name = "y_{LAB}";
  const TString pT_name = "pT, GeV/c";
  const TString phi_name = "#varphi, rad";
  const int y_nbins = 11;
  const int pT_nbins = 13;
  const int phi_nbins = 15;
  const float y_beam = 1.62179;
  const float y_low = y_beam - 0.8;
  const float y_up = y_beam + 1.4;
  const float pT_low = 0.;
  const float pT_up = 2.6;
  const float phi_low = -TMath::Pi();
  const float phi_up = TMath::Pi();

  TFile* reco_file = TFile::Open(recofilename);
  TTree* reco_tree = reco_file->Get<TTree>("aTree");
  
  auto* reco_tracks = new AnalysisTree::Particles();
  
  reco_tree -> SetBranchAddress("Candidates", &reco_tracks);
  
  TH2F hrec_y_pt("hrec_y_pt", "", y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);
  TH2F hrec_phi_pt("hrec_phi_pt", "", phi_nbins, phi_low, phi_up, pT_nbins, pT_low, pT_up);
  TH2F hrec_phi_y("hrec_phi_y", "", phi_nbins, phi_low, phi_up, y_nbins, y_low, y_up);
  
  hrec_y_pt.GetXaxis()->SetTitle(y_name);
  hrec_phi_pt.GetXaxis()->SetTitle(phi_name);
  hrec_phi_y.GetXaxis()->SetTitle(phi_name);
  
  hrec_y_pt.GetYaxis()->SetTitle(pT_name);
  hrec_phi_pt.GetYaxis()->SetTitle(pT_name);
  hrec_phi_y.GetYaxis()->SetTitle(y_name);
  
  const int rec_entries = reco_tree->GetEntries();
  for(int iEvent=0; iEvent<rec_entries; iEvent++)
  {
    reco_tree->GetEntry(iEvent);
    for(const auto& recotrack : *(reco_tracks->GetChannels()))
    {
      if(recotrack.GetPid() != 3312) continue;
      if(recotrack.GetField<int>(3) == 0) continue;
      
      hrec_y_pt.Fill(recotrack.GetRapidity(), recotrack.GetPt());
      hrec_phi_pt.Fill(recotrack.GetPhi(), recotrack.GetPt());
      hrec_phi_y.Fill(recotrack.GetPhi(), recotrack.GetRapidity());
    }
  }
  
  TFile fileOut("recmap_pt_y_phi.root", "recreate");
  fileOut.cd();
  
  hrec_y_pt.Write();
  hrec_phi_pt.Write();
  hrec_phi_y.Write();
  
  fileOut.Close();
  
  std::cout << "Macro finished successfully\n";
}
  