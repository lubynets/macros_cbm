void v1_builder(const TString& fileName)
{
//   std::string component = "tracks_centrality";
  std::string component = "pT";
//   std::string component = "rapidity";
  
  TFile* v1file = TFile::Open(fileName);
  auto* mc_lambda_psi_xx = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)v1file -> Get("sim/u_sim_PLAIN.Q_psi_PLAIN.x1x1");
  auto* lambda_psi_xx = (Qn::DataContainer<Qn::StatCollect,Qn::Axis<double>>*)v1file -> Get("rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN.x1x1");
  auto mc_projection = mc_lambda_psi_xx->Projection({"SimParticles_"+component});
  auto projection = lambda_psi_xx->Projection({"RecParticlesMcPid_"+component});
//   auto mc_projection = mc_lambda_psi_xx->Projection({"AnaEventHeader_tracks_centrality"});
//   auto projection = lambda_psi_xx->Projection({"AnaEventHeader_tracks_centrality"});
  auto* mc_graph = Qn::DataContainerHelper::ToTGraph(mc_projection);
  auto* graph = Qn::DataContainerHelper::ToTGraph(projection);
  
  if(component=="tracks_centrality")
    mc_graph -> GetXaxis() -> SetTitle("Centrality, %");
  if(component=="pT")
    mc_graph -> GetXaxis() -> SetTitle("p_{T}, GeV");
  if(component=="rapidity")
    mc_graph -> GetXaxis() -> SetTitle("y_{LAB}");
  graph -> GetXaxis() -> SetTitle(mc_graph -> GetXaxis() -> GetTitle());
  
  mc_graph -> SetName("mc");
  graph -> SetName("rec");
  
  const int N = graph->GetN();
  std::vector<double> x;
  std::vector<double> mc_y;
  std::vector<double> mc_y_err;
  std::vector<double> y;
  std::vector<double> y_err;
  x.resize(N);
  mc_y.resize(N);
  mc_y_err.resize(N);
  y.resize(N);
  y_err.resize(N);
    
  for(int i=0; i<N; i++)
  {
    mc_graph->GetPoint(i, x.at(i), mc_y.at(i));
    graph->GetPoint(i, x.at(i), y.at(i));
    mc_y.at(i) = 2*mc_y.at(i);
    y.at(i) = 2*y.at(i);
    mc_graph->SetPoint(i, x.at(i), mc_y.at(i));
    graph->SetPoint(i, x.at(i), y.at(i));
    mc_graph->SetPointError(i, 0, 2*mc_graph->GetErrorY(i));
    graph->SetPointError(i, 0, 2*graph->GetErrorY(i));
  }
  
  mc_graph -> SetMarkerColor(kRed);
  mc_graph -> SetLineColor(kRed);
  mc_graph -> SetLineWidth(2);
  graph -> SetMarkerColor(kBlack);
  graph -> SetLineColor(kBlack);
  graph -> SetLineWidth(2);
  
  mc_graph->Draw("APL");
  graph->Draw("PL");
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  mc_graph->Write();
  graph->Write();
  fileOut->Close();
  
}