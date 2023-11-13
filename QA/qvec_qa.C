template<typename T>
int FindBin(std::vector<T>& v, float value, bool is_spectator=true);
std::pair<float, float> PhiW(Qn::QVector v);

void qvec_qa(const std::string&& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  TTree* treeIn = fileIn->Get<TTree>("tree");

  std::vector<double> C_edges{0, 5, 10, 15, 20, 30, 40, 70};
  const int C_nbins = C_edges.size() - 1;

  Qn::DataContainer<Qn::QVector>* dc_usim{nullptr};
  Qn::DataContainer<Qn::QVector>* dc_urec_plain{nullptr};
  Qn::DataContainer<Qn::QVector>* dc_urec_recentered{nullptr};
  Qn::DataContainer<Qn::QVector>* dc_urec_twist{nullptr};
  Qn::DataContainer<Qn::QVector>* dc_urec_rescaled{nullptr};
  double centrality;

//   treeIn->SetBranchAddress("u_sim_PLAIN", &dc_usim);
//   treeIn->SetBranchAddress("u_rec_sgnl_PLAIN", &dc_urec_plain);
//   treeIn->SetBranchAddress("u_rec_sgnl_RECENTERED", &dc_urec_recentered);
//   treeIn->SetBranchAddress("u_rec_sgnl_TWIST", &dc_urec_twist);
//   treeIn->SetBranchAddress("u_rec_sgnl_RESCALED", &dc_urec_rescaled);
  treeIn->SetBranchAddress("u_sim_prim_PLAIN", &dc_usim);
  treeIn->SetBranchAddress("u_rec_sgnl_prim_PLAIN", &dc_urec_plain);
  treeIn->SetBranchAddress("u_rec_sgnl_prim_RECENTERED", &dc_urec_recentered);
  treeIn->SetBranchAddress("u_rec_sgnl_prim_TWIST", &dc_urec_twist);
  treeIn->SetBranchAddress("u_rec_sgnl_prim_RESCALED", &dc_urec_rescaled);
  treeIn->SetBranchAddress("RecEventHeader_centrality_tracks", &centrality);
  treeIn->GetEntry(0);

  Qn::DataContainer<TH1F, Qn::AxisD>* th1_usim = new Qn::DataContainer<TH1F, Qn::AxisD>();
  th1_usim->AddAxes({{"centrality", C_edges}});
  th1_usim->AddAxes(dc_usim->GetAxes());

  Qn::DataContainer<TH1F, Qn::AxisD>* th1_rec_sgnl_PLAIN = new Qn::DataContainer<TH1F, Qn::AxisD>();
  th1_rec_sgnl_PLAIN->AddAxes({{"centrality", C_edges}});
  th1_rec_sgnl_PLAIN->AddAxes(dc_usim->GetAxes());

  Qn::DataContainer<TH1F, Qn::AxisD>* th1_rec_sgnl_RECENTERED = new Qn::DataContainer<TH1F, Qn::AxisD>();
  th1_rec_sgnl_RECENTERED->AddAxes({{"centrality", C_edges}});
  th1_rec_sgnl_RECENTERED->AddAxes(dc_usim->GetAxes());

  Qn::DataContainer<TH1F, Qn::AxisD>* th1_rec_sgnl_TWIST = new Qn::DataContainer<TH1F, Qn::AxisD>();
  th1_rec_sgnl_TWIST->AddAxes({{"centrality", C_edges}});
  th1_rec_sgnl_TWIST->AddAxes(dc_usim->GetAxes());

  Qn::DataContainer<TH1F, Qn::AxisD>* th1_rec_sgnl_RESCALED = new Qn::DataContainer<TH1F, Qn::AxisD>();
  th1_rec_sgnl_RESCALED->AddAxes({{"centrality", C_edges}});
  th1_rec_sgnl_RESCALED->AddAxes(dc_usim->GetAxes());


  for(int i=0; i<th1_usim->size(); i++) {
    th1_usim->At(i).SetBins(1200, -TMath::Pi(), TMath::Pi());
    th1_rec_sgnl_PLAIN->At(i).SetBins(1200, -TMath::Pi(), TMath::Pi());
    th1_rec_sgnl_RECENTERED->At(i).SetBins(1200, -TMath::Pi(), TMath::Pi());
    th1_rec_sgnl_TWIST->At(i).SetBins(1200, -TMath::Pi(), TMath::Pi());
    th1_rec_sgnl_RESCALED->At(i).SetBins(1200, -TMath::Pi(), TMath::Pi());

  }

  const int nEvents = treeIn->GetEntries();

  for(int iEvent=0; iEvent<nEvents; iEvent++) {
    if(iEvent%1000==0)
      std::cout << iEvent << "\n";
    treeIn->GetEntry(iEvent);
    const size_t C_bin = FindBin(C_edges, centrality, false);
    if(C_bin<0 || C_bin>C_nbins-1) continue;

    for(int iBin=0; iBin<dc_usim->size(); iBin++) {
      std::vector<size_t> indices = dc_usim->GetIndex(iBin);

      Qn::QVector v_usim = dc_usim->At(iBin);
      auto phiw_usim = PhiW(v_usim);
      th1_usim->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw_usim.first, phiw_usim.second);

      Qn::QVector v_urec_plain = dc_urec_plain->At(iBin);
      auto phiw_urec_plain = PhiW(v_urec_plain);
      th1_rec_sgnl_PLAIN->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw_urec_plain.first, phiw_urec_plain.second);

      Qn::QVector v_urec_recentered = dc_urec_recentered->At(iBin);
      auto phiw_urec_recentered = PhiW(v_urec_recentered);
      th1_rec_sgnl_RECENTERED->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw_urec_recentered.first, phiw_urec_recentered.second);

      Qn::QVector v_urec_twist = dc_urec_twist->At(iBin);
      auto phiw_urec_twist = PhiW(v_urec_twist);
      th1_rec_sgnl_TWIST->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw_urec_twist.first, phiw_urec_twist.second);

      Qn::QVector v_urec_rescaled = dc_urec_rescaled->At(iBin);
      auto phiw_urec_rescaled = PhiW(v_urec_rescaled);
      th1_rec_sgnl_RESCALED->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw_urec_rescaled.first, phiw_urec_rescaled.second);
    }
  }

  TFile* fileOut = TFile::Open("qvec_qa.root", "recreate");
  fileOut->cd();
  th1_usim->Write("th1_usim");
  th1_rec_sgnl_PLAIN->Write("th1_rec_sgnl_PLAIN");
  th1_rec_sgnl_RECENTERED->Write("th1_rec_sgnl_RECENTERED");
  th1_rec_sgnl_TWIST->Write("th1_rec_sgnl_TWIST");
  th1_rec_sgnl_RESCALED->Write("th1_rec_sgnl_RESCALED");
  fileOut->Close();
}


template<typename T>
int FindBin(std::vector<T>& v, float value, bool is_spectator) {
  if(!std::is_sorted(v.begin(), v.end())) {
    throw std::runtime_error("Vector of centrality is not ordered!");
  }

  int bin = -1;
  for(auto& ele : v) {
    if(value>ele)
      bin++;
    else
      break;
  }

  return is_spectator ? bin : v.size()-2-bin;
}

std::pair<float, float> PhiW(Qn::QVector v) {
  const float x = v.x(1);
  const float y = v.y(1);
  const float w = v.sumweights();

  const float phi = TMath::ATan2(y, x);

  return std::make_pair(phi, w);
}
