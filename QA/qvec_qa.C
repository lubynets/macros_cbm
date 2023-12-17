template<typename T>
int FindBin(std::vector<T>& v, float value, bool is_spectator=true);
std::pair<float, float> PhiW(Qn::QVector v);

void qvec_qa(const std::string&& fileName) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  TTree* treeIn = fileIn->Get<TTree>("tree");

  const int nbins = 1200;

  std::vector<double> C_edges{0, 5, 10, 15, 20, 30, 40, 70};
  const int C_nbins = C_edges.size() - 1;

  enum eStep : short {
    kSim = 0,
    kPlain,
    kRecentered,
    kTwist,
    kRescaled,
    kNumberOfSteps
  };

  struct Step {
    std::string inname_;
    std::string outname_;
  };

//   std::string is_prim = "";
  std::string is_prim = "_prim_";

  std::vector<Step> steps {
    {"u_sim" + is_prim + "PLAIN", "sim"},
    {"u_rec_sgnl" + is_prim + "PLAIN", "plain"},
    {"u_rec_sgnl" + is_prim + "RECENTERED", "recenter"},
    {"u_rec_sgnl" + is_prim + "TWIST", "twist"},
    {"u_rec_sgnl" + is_prim + "RESCALED", "rescale"}
  };

  typedef Qn::DataContainer<Qn::QVector> QnDcV;
  typedef Qn::DataContainer<TH1F, Qn::AxisD> QnDcTH1F;

  std::vector<QnDcV*> dcIn(kNumberOfSteps, nullptr);
  std::vector<QnDcTH1F*> dcTH1(kNumberOfSteps);

  double centrality;

  for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
    treeIn->SetBranchAddress(steps.at(iStep).inname_.c_str(), &dcIn.at(iStep));
  }
  treeIn->SetBranchAddress("RecEventHeader_centrality_tracks", &centrality);
  treeIn->GetEntry(0);

  for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
    dcTH1.at(iStep) = new QnDcTH1F();
    dcTH1.at(iStep)->AddAxes({{"centrality", C_edges}});
    dcTH1.at(iStep)->AddAxes(dcIn.at(0)->GetAxes());
    for(int iCell=0; iCell<dcTH1.at(0)->size(); iCell++) {
      dcTH1.at(iStep)->At(iCell).SetBins(nbins, -TMath::Pi(), TMath::Pi());
    }
  }

  const int nEvents = treeIn->GetEntries();

  for(int iEvent=0; iEvent<nEvents; iEvent++) {
    if(iEvent%1000==0)
      std::cout << iEvent << "\n";
    treeIn->GetEntry(iEvent);
    const size_t C_bin = FindBin(C_edges, centrality, false);
    if(C_bin<0 || C_bin>C_nbins-1) continue;

    for(int iBin=0; iBin<dcIn.at(0)->size(); iBin++) {
      std::vector<size_t> indices = dcIn.at(0)->GetIndex(iBin);

      for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
        Qn::QVector v = dcIn.at(iStep)->At(iBin);
        auto phiw = PhiW(v);
        dcTH1.at(iStep)->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw.first, phiw.second);
      }
    }
  }

  TFile* fileOut = TFile::Open("qvec_qa.root", "recreate");
  fileOut->cd();
  for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
    dcTH1.at(iStep)->Write(("phi_" + steps.at(iStep).outname_).c_str());
  }
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
