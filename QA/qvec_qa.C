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

  std::vector<std::string> subevents{"psd1", "psd2", "psd3"};

  typedef Qn::DataContainer<Qn::QVector> QnDcV;
  typedef Qn::DataContainer<TH1F, Qn::AxisD> QnDcTH1F;

  std::vector<QnDcV*> dcObs(kNumberOfSteps, nullptr);
  std::vector<QnDcTH1F*> dcPhi_o(kNumberOfSteps);
  std::vector<QnDcTH1F*> dcX_o(kNumberOfSteps);
  std::vector<QnDcTH1F*> dcY_o(kNumberOfSteps);

  std::vector<QnDcV*> dcRef(subevents.size(), nullptr);
  std::vector<QnDcTH1F*> dcPhi_r(subevents.size());
  std::vector<QnDcTH1F*> dcMod_r(subevents.size());
  std::vector<QnDcTH1F*> dcX_r(subevents.size());
  std::vector<QnDcTH1F*> dcY_r(subevents.size());

  double centrality;

  for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
    treeIn->SetBranchAddress(steps.at(iStep).inname_.c_str(), &dcObs.at(iStep));
  }
  for(int iSub=0; iSub<subevents.size(); iSub++) {
    treeIn->SetBranchAddress((subevents.at(iSub) + "_PLAIN").c_str(), &dcRef.at(iSub));
  }
  treeIn->SetBranchAddress("RecEventHeader_centrality_tracks", &centrality);
  treeIn->GetEntry(0);

  for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
    int idc{0};
    std::vector<std::vector<QnDcTH1F*>*> dcVec{&dcPhi_o, &dcX_o, &dcY_o};
    for(auto& dc : dcVec) {
      dc->at(iStep) = new QnDcTH1F();
      dc->at(iStep)->AddAxes({{"centrality", C_edges}});
      dc->at(iStep)->AddAxes(dcObs.at(0)->GetAxes());
      const float lim = idc == 0 ? TMath::Pi() : 1.2;
      for(int iCell=0; iCell<dc->at(0)->size(); iCell++) {
        dc->at(iStep)->At(iCell).SetBins(nbins, -lim, lim);
      }
      idc++;
    }
  }
  for(int iSub=0; iSub<subevents.size(); iSub++) {
    int idc{0};
    std::vector<std::vector<QnDcTH1F*>*> dcVec{&dcPhi_r, &dcMod_r, &dcX_r, &dcY_r};
    for(auto& dc : dcVec) {
      dc->at(iSub) = new QnDcTH1F();
      dc->at(iSub)->AddAxes({{"centrality", C_edges}});
      const float lim = idc == 0 ? TMath::Pi() : 1.2;
      for(int iCell=0; iCell<dc->at(0)->size(); iCell++) {
        dc->at(iSub)->At(iCell).SetBins(nbins, -lim, lim);
      }
      idc++;
    }
  }


  const int nEvents = treeIn->GetEntries();

  for(int iEvent=0; iEvent<nEvents; iEvent++) {
    if(iEvent%1000==0)
      std::cout << iEvent << "\n";
    treeIn->GetEntry(iEvent);
    const size_t C_bin = FindBin(C_edges, centrality, false);
    if(C_bin<0 || C_bin>C_nbins-1) continue;

    for(int iBin=0; iBin<dcObs.at(0)->size(); iBin++) {
      std::vector<size_t> indices = dcObs.at(0)->GetIndex(iBin);

      for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
        Qn::QVector v = dcObs.at(iStep)->At(iBin);
        auto phiw = PhiW(v);
        dcPhi_o.at(iStep)->At({C_bin, indices.at(0), indices.at(1)}).Fill(phiw.first, phiw.second);
        dcX_o.at(iStep)->At({C_bin, indices.at(0), indices.at(1)}).Fill(v.x(1), v.sumweights());
        dcY_o.at(iStep)->At({C_bin, indices.at(0), indices.at(1)}).Fill(v.y(1), v.sumweights());
      }
    }
    for(int iSub=0; iSub<subevents.size(); iSub++) {
      Qn::QVector v = dcRef.at(iSub)->At(0);
      auto phiw = PhiW(v);
      dcPhi_r.at(iSub)->At(C_bin).Fill(phiw.first, phiw.second);
      dcMod_r.at(iSub)->At(C_bin).Fill(std::sqrt(v.x(1)*v.x(1) + v.y(1)*v.y(1)), v.sumweights());
      dcX_r.at(iSub)->At(C_bin).Fill(v.x(1), v.sumweights());
      dcY_r.at(iSub)->At(C_bin).Fill(v.y(1), v.sumweights());
    }
  }

  TFile* fileOut = TFile::Open("qvec_qa.root", "recreate");
  fileOut->cd();
//   for(int iStep=0; iStep<kNumberOfSteps; iStep++) {
//     dcPhi_o.at(iStep)->Write(("phi_" + steps.at(iStep).outname_).c_str());
//     dcX_o.at(iStep)->Write(("x_" + steps.at(iStep).outname_).c_str());
//     dcY_o.at(iStep)->Write(("y_" + steps.at(iStep).outname_).c_str());
//   }
  for(int iSub=0; iSub<subevents.size(); iSub++) {
    dcPhi_r.at(iSub)->Write(("phi_" + subevents.at(iSub)).c_str());
    dcMod_r.at(iSub)->Write(("mod_" + subevents.at(iSub)).c_str());
    dcX_r.at(iSub)->Write(("x_" + subevents.at(iSub)).c_str());
    dcY_r.at(iSub)->Write(("y_" + subevents.at(iSub)).c_str());
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
