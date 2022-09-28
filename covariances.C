template<typename T>
int FindBin(std::vector<T>& v, float value);

struct MyQVec{

  MyQVec() = default;
  MyQVec(Qn::QVector v) {
    x_ = v.x(1);
    y_ = v.y(1);
    sumw_ = v.sumweights();
  };

  float x_;
  float y_;
  float sumw_;
};

MyQVec RotateRight(MyQVec mqv);
float ScalarProduct(MyQVec v1, MyQVec v2);

// void covariances(const std::string&& fileName) {
void covariances(const std::string&& particle) {
//   TFile* fileIn = TFile::Open(fileName.c_str(), "read");
//   TTree* treeIn = fileIn->Get<TTree>("tree");

  TChain* treeIn = new TChain("tree");
  std::string inputPath = "/lustre/cbm/users/lubynets/qna/outputs/sim_tracks_flow/";
  for(int iFile=1; iFile<=100; iFile++) {
    std::string fileNum = std::to_string(iFile);
    std::string fileName = inputPath + fileNum + "/correction_out_1.root";
    treeIn->Add(fileName.c_str());
  }

  std::vector<double> C_edges{0, 4.2, 6.0, 7.3, 8.4, 9.4, 10.3, 11.2};
  const int C_nbins = C_edges.size() - 1;

//   std::string particle = "lambda";
//   std::string particle = "pipos";
//   std::string particle = "pineg";

  Qn::DataContainer<Qn::QVector>* dc_usim{nullptr};
  Qn::DataContainer<Qn::QVector>* dc_qspec{nullptr};
  Qn::DataContainer<Qn::QVector>* dc_qpsi{nullptr};
  double impactpar;

  treeIn->SetBranchAddress(("u_sim_" + particle + "_PLAIN").c_str(), &dc_usim);
  treeIn->SetBranchAddress("Q_spec_PLAIN", &dc_qspec);
  treeIn->SetBranchAddress("Q_psi_PLAIN", &dc_qpsi);
  treeIn->SetBranchAddress("SimEventHeader_b", &impactpar);
  treeIn->GetEntry(0);

  const int   par_nbins_u = 200;
  const int   par_nbins_Q = 200;
  const float par_low_u = -1.0;
  const float par_up_u = 1.0;
  const float par_low_Q = -1.0;
  const float par_up_Q = 1.0;

  const int   perp_nbins_u = 200;
  const int   perp_nbins_Q = 200;
  const float perp_low_u = -1.0;
  const float perp_up_u = 1.0;
  const float perp_low_Q = -1.0;
  const float perp_up_Q = 1.0;

  Qn::DataContainer<TH2F, Qn::AxisD>* dc_par = new Qn::DataContainer<TH2F, Qn::AxisD>();
  dc_par->AddAxes({{"impactpar", C_edges}});
  dc_par->AddAxes(dc_usim->GetAxes());

  Qn::DataContainer<TH2F, Qn::AxisD>* dc_perp = new Qn::DataContainer<TH2F, Qn::AxisD>();
  dc_perp->AddAxes({{"impactpar", C_edges}});
  dc_perp->AddAxes(dc_usim->GetAxes());

  for(int i=0; i<dc_par->size(); i++) {
    dc_par->At(i).SetName("hpar");
    dc_par->At(i).SetTitle("hPar");
    dc_par->At(i).SetBins(par_nbins_u, par_low_u, par_up_u, par_nbins_Q, par_low_Q, par_up_Q);
    dc_par->At(i).GetXaxis()->SetTitle("u_{#parallel}");
    dc_par->At(i).GetYaxis()->SetTitle("Q_{#parallel}");

    dc_perp->At(i).SetName("hperp");
    dc_perp->At(i).SetTitle("hPerp");
    dc_perp->At(i).SetBins(perp_nbins_u, perp_low_u, perp_up_u, perp_nbins_Q, perp_low_Q, perp_up_Q);
    dc_perp->At(i).GetXaxis()->SetTitle("u_{#perp}");
    dc_perp->At(i).GetYaxis()->SetTitle("Q_{#perp}");
  }

  const int nEvents = treeIn->GetEntries();
//   const int nEvents = 10;

  for(int iEvent=0; iEvent<nEvents; iEvent++) {
    if(iEvent%1000==0)
      std::cout << iEvent << "\n";
    treeIn->GetEntry(iEvent);
    const size_t C_bin = FindBin(C_edges, impactpar);
    if(C_bin<0 || C_bin>C_nbins-1) continue;
    MyQVec qspec = MyQVec(dc_qspec->At(0));
    MyQVec qpsi = MyQVec(dc_qpsi->At(0));
    MyQVec qpsiR = RotateRight(qpsi);
    const float Qpar = ScalarProduct(qspec, qpsi);
    const float Qperp = ScalarProduct(qspec, qpsiR);

    for(int iBin=0; iBin<dc_usim->size(); iBin++) {
      std::vector<size_t> indices = dc_usim->GetIndex(iBin);
      MyQVec qsim = MyQVec(dc_usim->At(iBin));

      const float upar = ScalarProduct(qsim, qpsi);
      const float uperp = ScalarProduct(qsim, qpsiR);

      dc_par->At({C_bin, indices.at(0), indices.at(1)}).Fill(upar, Qpar, qsim.sumw_);
      dc_perp->At({C_bin, indices.at(0), indices.at(1)}).Fill(uperp, Qperp, qsim.sumw_);
    }
  }

  TFile* fileOut = TFile::Open(("covariances." + particle + ".root").c_str(), "recreate");
  fileOut->cd();
  dc_par->Write("par");
  dc_perp->Write("perp");
  fileOut->Close();
}

template<typename T>
int FindBin(std::vector<T>& v, float value) {
  if(!std::is_sorted(v.begin(), v.end()))
    throw std::runtime_error("Vector of centrality is not ordered!");

  int bin = -1;
  for(auto& ele : v) {
    if(value>ele)
      bin++;
    else
      break;
  }

  return bin;
}

MyQVec RotateRight(MyQVec mqv) {
  MyQVec result;
  result.x_ = mqv.y_;
  result.y_ = -mqv.x_;
  result.sumw_ = mqv.sumw_;

  return result;
}

float ScalarProduct(MyQVec v1, MyQVec v2) {
  return v1.x_*v2.x_ + v1.y_*v2.y_;
}


