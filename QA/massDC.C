void massDC(std::string filelist_eve, std::string filelist_pfs, int pdg=3122)
{
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_eve, filelist_pfs}), std::vector<std::string>({"aTree", "pTree"}));
  auto* config = treeIn->GetConfiguration();
  
  auto* reco_tracks = new AnalysisTree::Particles();
  auto* event_header = new AnalysisTree::EventHeader();
   
  treeIn->SetBranchAddress("Candidates", &reco_tracks);
  treeIn->SetBranchAddress("AnaEventHeader", &event_header);
  
  int centrality_id_ = config->GetBranchConfig("AnaEventHeader").GetFieldId("centrality_tracks");
  int is_signal_id_ = config->GetBranchConfig("Candidates").GetFieldId("generation");
  
  const float y_beam = 1.62179;
  
  const int   y_nbins = 4;
  const float y_low = y_beam-0.6;
  const float y_up = y_beam+1.0;

  const int   pT_nbins = 4;
  const float pT_low = 0.2;
  const float pT_up = 1.4;
  
  std::vector<double> C_binranges {0, 10, 20, 30, 40, 70, 100};
  const int C_nbins = C_binranges.size()-1;
  
  Qn::Axis<double> C_axis("centrality", C_binranges);
  Qn::Axis<double> pT_axis("pT", pT_nbins, pT_low, pT_up);
  Qn::Axis<double> y_axis("y", y_nbins, y_low, y_up);
  
  Qn::DataContainer<TH1F, Qn::Axis<double>> dcmass;
  dcmass.AddAxes({C_axis, pT_axis, y_axis});
  
  for(int i=0; i<dcmass.size(); i++) {
    dcmass[i].SetName("hmass");
    dcmass[i].SetTitle("hMass");
    dcmass[i].SetBins(1200, 1.045, 1.185);
    dcmass[i].GetXaxis()->SetTitle("m_{inv}, GeV");
  }
  
  TFile* fileOut = TFile::Open("massDC.root", "recreate");
  
  const int n_entries = treeIn->GetEntries();
  for(int i_event=0; i_event<n_entries; ++i_event)
  {
    if(i_event%100 == 0)
      std::cout << i_event << "\n";
    treeIn->GetEntry(i_event);
    const float centrality = event_header->GetField<float>(centrality_id_);
    for(const auto& reco_track : *(reco_tracks->GetChannels()) )
    {     
      const int pid = reco_track.GetPid();
      if(pid != pdg) continue;
      
//       if(!(reco_track.GetField<int>(is_signal_id_)==1 || reco_track.GetField<int>(is_signal_id_)==2)) continue; // ONLY signal
//       if(!(reco_track.GetField<int>(is_signal_id_)==0)) continue;                                               // ONLY background
      const float pT_reco = reco_track.GetPt();
      const float y_reco = reco_track.GetRapidity();
      
      const int linearindex = dcmass.FindBin(std::vector<double>({centrality, pT_reco, y_reco}));
      if(linearindex<0) continue;
      dcmass[linearindex].Fill(reco_track.GetMass());
    }    
  }  
      
  dcmass.Write("dcmass");
  
  fileOut -> Close();
}
