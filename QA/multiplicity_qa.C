void multiplicity_qa(std::vector<std::string> filelist)
{
  std::cout << "Macro started" << std::endl;
  
  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(filelist, std::vector<std::string>({"rTree"}));
  
  auto* conf = treeIn->GetConfiguration();
  auto* vtx_tracks = new AnalysisTree::TrackDetector();
  treeIn->SetBranchAddress("VtxTracks", &vtx_tracks);
  
  AnalysisTree::SimpleCut vtx_chi2_track_cut = AnalysisTree::RangeCut("VtxTracks.vtx_chi2", 0, 3);
  AnalysisTree::SimpleCut nhits_cut          = AnalysisTree::RangeCut("VtxTracks.nhits", 4, 100);
  AnalysisTree::SimpleCut chi2_cut({"VtxTracks.chi2", "VtxTracks.ndf"},
                                   [](std::vector<double> par) { return par[0] / par[1] < 3; });
  AnalysisTree::SimpleCut eta_cut = AnalysisTree::RangeCut("VtxTracks.eta", 0.2, 6);

  auto* vertex_tracks_cuts =
    new AnalysisTree::Cuts("VtxTracks", {vtx_chi2_track_cut, nhits_cut, chi2_cut, eta_cut});
    
  vertex_tracks_cuts->Init(*conf);
    
  const int Nevents = treeIn->GetEntries();
  TFile* fileOut = TFile::Open("out.mult.root", "recreate");
  TH1F hMult("hMult", "", 1000, 0, 1000);  
  hMult.GetXaxis()->SetTitle("Multiplicity");
  hMult.GetYaxis()->SetTitle("Entries");
  
  for(int i=0; i<Nevents; i++)
  {
    treeIn -> GetEntry(i);
    int m = 0;
    
    for(const auto& vtx_track : *(vtx_tracks->GetChannels()) )
    {
      if(vertex_tracks_cuts->Apply(vtx_track))
        m++;      
    }
    
    hMult.Fill(m);
  }
  
  hMult.Write();
  fileOut->Close();  
  
  std::cout << "Macro finished successfully" << std::endl;
}