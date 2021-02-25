std::string StringBinNumber(int number);

void mass3D(const TString& infile)
{
  TFile* file = TFile::Open(infile);
  TTree* tree = file->Get<TTree>("cTree");
  
  auto* reco_tracks = new AnalysisTree::Particles();
  auto* event_header = new AnalysisTree::EventHeader();
  
  AnalysisTree::Configuration* config = (AnalysisTree::Configuration*) file->Get("Configuration");
  
  tree->SetBranchAddress("RecParticlesMcPid", &reco_tracks);
  tree->SetBranchAddress("AnaEventHeader", &event_header);
  
  int centrality_id_ = config->GetBranchConfig("AnaEventHeader").GetFieldId("tracks_centrality");
  int is_signal_id_ = config->GetBranchConfig("RecParticlesMcPid").GetFieldId("is_signal");
  
//   std::cout << centrality_id_ << "\t" << is_signal_id_ << "\n";
  
  const float y_beam = 1.62179;
  
  const int   y_nbins = 5;
  const float y_low = y_beam-0.45;
  const float y_up = y_beam+1.05;
  const float y_bin_width = (y_up-y_low)/y_nbins;

  const int   pT_nbins = 5;
  const float pT_low = 0.3;
  const float pT_up = 1.3;
  const float pT_bin_width = (pT_up-pT_low)/pT_nbins;
  
  std::vector<double> C_binranges {0, 20, 40, 100};
  const int C_nbins = C_binranges.size()-1;
  
  double* C_edges = &C_binranges[0]; 
  double* y_edges = new double[y_nbins+1];
  double* pT_edges = new double[pT_nbins+1];
  
  for(int i=0; i<=y_nbins; i++)
    y_edges[i] = y_low + i*y_bin_width;
  for(int i=0; i<=pT_nbins; i++)
    pT_edges[i] = pT_low + i*pT_bin_width;
  
//   for(int i=0; i<=C_nbins; i++)
//     std::cout << C_edges[i] << "\t";
//   std::cout << "\n";
//   for(int i=0; i<=y_nbins; i++)
//     std::cout << y_edges[i] << "\t";
//   std::cout << "\n";
//   for(int i=0; i<=pT_nbins; i++)
//     std::cout << pT_edges[i] << "\t";
//   std::cout << "\n";
  
  TH3C hframe("hframe", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
  
  std::vector<TH1F*> hmass;
  hmass.resize(hframe.GetNcells());
  std::vector<bool> histo_is_inside;
  histo_is_inside.resize(hframe.GetNcells());
  for(auto is : histo_is_inside)
    is = false;  
  
  for(int iC = 1; iC<=C_nbins; iC++)
    for(int iy = 1; iy<=y_nbins; iy++)
      for(int ipT = 1; ipT<=pT_nbins; ipT++)
      {
        std::string histoname = "C" + StringBinNumber(iC) + "_y" + StringBinNumber(iy) + "_pT" + StringBinNumber(ipT);
        int binnumber = hframe.GetBin(iC, iy, ipT);
        std::string histotitle = std::to_string(binnumber);
        hmass.at(binnumber) = new TH1F(histoname.c_str(), histotitle.c_str(), 600, 1.08, 1.15);
        hmass.at(binnumber) -> GetXaxis() -> SetTitle("m_{inv}, GeV");
        histo_is_inside.at(binnumber) = true;
      }
  
  TFile* fileOut = TFile::Open("out.mass3D.root", "recreate");
  
  const int n_entries = tree->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    tree->GetEntry(i_event);
    const float centrality = event_header->GetField<float>(centrality_id_);
    for(const auto& reco_track : *(reco_tracks->GetChannels()) )
    {
      if(!(reco_track.GetField<int>(is_signal_id_)==1 || reco_track.GetField<int>(is_signal_id_)==2)) continue;
      const float pT_reco = reco_track.GetPt();
      const float y_reco = reco_track.GetRapidity();
      
      const int binnumber = hframe.FindBin(centrality, y_reco, pT_reco);
      if(histo_is_inside.at(binnumber) == true)
        hmass.at(binnumber) -> Fill(reco_track.GetMass());
    }    
  }  
  
  for(int iC = 1; iC<=C_nbins; iC++)
    for(int iy = 1; iy<=y_nbins; iy++)
      for(int ipT = 1; ipT<=pT_nbins; ipT++)
      {
        int binnumber = hframe.GetBin(iC, iy, ipT);
        hmass.at(binnumber) -> Write();
      }  
  fileOut -> Close();
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}