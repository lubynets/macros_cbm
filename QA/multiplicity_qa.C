void multiplicity_qa(std::vector<std::string> filelist)
{
  std::cout << "Macro started" << std::endl;
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(filelist, std::vector<std::string>({"rTree"}));
  
  auto* conf = treeIn->GetConfiguration();
  auto* sim_event_header = new AnalysisTree::EventHeader();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  auto* rec_event_header = new AnalysisTree::EventHeader();
  treeIn->SetBranchAddress("SimEventHeader.", &sim_event_header);
  treeIn->SetBranchAddress("RecEventHeader.", &rec_event_header);
  treeIn->SetBranchAddress("VtxTracks.", &vtx_tracks);
  
  const int b_id = conf->GetBranchConfig("SimEventHeader").GetFieldId("b");
  const int epsd_id = conf->GetBranchConfig("RecEventHeader").GetFieldId("Epsd");
  
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
  TH1F hE("hE", "", 600, 0, 60);
  TH1F hB("hB", "", 600, 0, 20);
  TH2F hB_Mult("hB_Mult", "", 600, 0, 20, 601, -0.5, 600.5);
  TH2F hB_E("hB_E", "", 600, 0, 20, 600, 0, 60);
  hB.GetXaxis()->SetTitle("b, fm");
  hB.GetYaxis()->SetTitle("Entries");
  hE.GetXaxis()->SetTitle("E_{PSD}, GeV");
  hE.GetYaxis()->SetTitle("Entries");
  hMult.GetXaxis()->SetTitle("Multiplicity");
  hMult.GetYaxis()->SetTitle("Entries");
  hB_Mult.GetXaxis()->SetTitle("b, fm");
  hB_Mult.GetYaxis()->SetTitle("Multiplicity");
  hB_E.GetXaxis()->SetTitle("b, fm");
  hB_E.GetYaxis()->SetTitle("E_{PSD}, GeV");

  
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
    hB_Mult.Fill(b, m);
    hB.Fill(b);

    const float e = rec_event_header->GetField<float>(epsd_id);

    hE.Fill(e);
    hB_E.Fill(b, e);
  }
  
  hMult.Write();
  hB.Write();
  hB_Mult.Write();
  hE.Write();
  hB_E.Write();
  fileOut->Close();  
  
  std::cout << "Macro finished successfully" << std::endl;
}
