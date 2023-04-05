void multiplicity_qa(std::vector<std::string> filelist)
{
  std::cout << "Macro started" << std::endl;
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(filelist, std::vector<std::string>({"rTree"}));
  
  auto* conf = treeIn->GetConfiguration();
  auto* sim_event_header = new AnalysisTree::EventHeader();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  treeIn->SetBranchAddress("SimEventHeader.", &sim_event_header);
  treeIn->SetBranchAddress("VtxTracks.", &vtx_tracks);
  
  const int b_id = conf->GetBranchConfig("SimEventHeader").GetFieldId("b");
  
//------------------- CbmGoodVertexTrack ------------------------------------------------------------  
  AnalysisTree::SimpleCut vtx_chi2_track_cut = AnalysisTree::RangeCut("VtxTracks.vtx_chi2", 0, 3);
  AnalysisTree::SimpleCut nhits_cut          = AnalysisTree::RangeCut("VtxTracks.nhits", 4, 100);
  AnalysisTree::SimpleCut chi2_cut({"VtxTracks.chi2", "VtxTracks.ndf"},
                                   [](std::vector<double> par) { return par[0] / par[1] < 3; });
  AnalysisTree::SimpleCut eta_cut = AnalysisTree::RangeCut("VtxTracks.eta", 0.2, 6);

  auto* vertex_tracks_cuts =
    new AnalysisTree::Cuts("VtxTracks", {vtx_chi2_track_cut, nhits_cut, chi2_cut, eta_cut});
//----------------------------------------------------------------------------------------------------
  
//   AnalysisTree::SimpleCut vtx_chi2_track_cut = AnalysisTree::RangeCut("VtxTracks.vtx_chi2", 0, 18);
//   
//   auto* vertex_tracks_cuts =
//     new AnalysisTree::Cuts("VtxTracks", {vtx_chi2_track_cut});
    
  vertex_tracks_cuts->Init(*conf);
    
  const int Nevents = treeIn->GetEntries();
  TFile* fileOut = TFile::Open("multiplicity_qa.root", "recreate");
  TH1F hMult("hMult", "", 1001, -0.5, 1000.5);
  TH1F hB("hB", "", 600, 0, 20);
  TH2F hCorr("hCorr", "", 601, -0.5, 600.5, 600, 0, 20);
  hMult.GetXaxis()->SetTitle("Multiplicity");
  hMult.GetYaxis()->SetTitle("Entries");
  hCorr.GetXaxis()->SetTitle("Multiplicity");
  hCorr.GetYaxis()->SetTitle("b, fm");
  
  
  for(int i=0; i<Nevents; i++)
  {
    treeIn -> GetEntry(i);
    int m = 0;
    const float b = sim_event_header->GetField<float>(b_id);
    
    for(const auto& vtx_track : *(vtx_tracks->GetChannels()) )
    {
      if(vertex_tracks_cuts->Apply(vtx_track))
        m++;      
    }
    
    hMult.Fill(m);
    hCorr.Fill(m, b);
    hB.Fill(b);
  }
  
  hMult.Write();
  hB.Write();
  hCorr.Write();
  fileOut->Close();  
  
  std::cout << "Macro finished successfully" << std::endl;
}
