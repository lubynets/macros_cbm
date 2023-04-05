#include "/lustre/cbm/users/lubynets/soft/QnTools/install/include/QnTools/ReSampleFunctor.hpp"

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


void covariances_scol(const std::string&& fileName, const std::string&& particle) {
  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  TTree* treeIn = fileIn->Get<TTree>("tree");


// void covariances_scol(const std::string&& particle) {
//   TChain* treeIn = new TChain("tree");
//   std::string inputPath = "/lustre/cbm/users/lubynets/qna/outputs/sim_tracks_flow/";
//   for(int iFile=1; iFile<=100; iFile++) {
//     std::string fileNum = std::to_string(iFile);
//     std::string fileName = inputPath + fileNum + "/correction_out_1.root";
//     treeIn->Add(fileName.c_str());
//   }

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

  const int nsam = 50;

  Qn::DataContainer<Qn::StatCollect, Qn::AxisD>* dc_u_par = new Qn::DataContainer<Qn::StatCollect, Qn::AxisD>();
  dc_u_par->AddAxes({{"impactpar", C_edges}});
  dc_u_par->AddAxes(dc_usim->GetAxes());
  for(int iBin=0; iBin<dc_u_par->size(); iBin++)
    dc_u_par->At(iBin).SetNumberOfSamples(nsam);

  Qn::DataContainer<Qn::StatCollect, Qn::AxisD>* dc_Q_par = new Qn::DataContainer<Qn::StatCollect, Qn::AxisD>();
  dc_Q_par->AddAxes({{"impactpar", C_edges}});
  for(int iBin=0; iBin<dc_Q_par->size(); iBin++)
    dc_Q_par->At(iBin).SetNumberOfSamples(nsam);

  Qn::DataContainer<Qn::StatCollect, Qn::AxisD>* dc_uQ_par = new Qn::DataContainer<Qn::StatCollect, Qn::AxisD>();
  dc_uQ_par->AddAxes({{"impactpar", C_edges}});
  dc_uQ_par->AddAxes(dc_usim->GetAxes());
  for(int iBin=0; iBin<dc_uQ_par->size(); iBin++)
    dc_uQ_par->At(iBin).SetNumberOfSamples(nsam);

  Qn::DataContainer<Qn::StatCollect, Qn::AxisD>* dc_u_perp = new Qn::DataContainer<Qn::StatCollect, Qn::AxisD>();
  dc_u_perp->AddAxes({{"impactpar", C_edges}});
  dc_u_perp->AddAxes(dc_usim->GetAxes());
  for(int iBin=0; iBin<dc_u_perp->size(); iBin++)
    dc_u_perp->At(iBin).SetNumberOfSamples(nsam);

  Qn::DataContainer<Qn::StatCollect, Qn::AxisD>* dc_Q_perp = new Qn::DataContainer<Qn::StatCollect, Qn::AxisD>();
  dc_Q_perp->AddAxes({{"impactpar", C_edges}});
  for(int iBin=0; iBin<dc_Q_perp->size(); iBin++)
    dc_Q_perp->At(iBin).SetNumberOfSamples(nsam);

  Qn::DataContainer<Qn::StatCollect, Qn::AxisD>* dc_uQ_perp = new Qn::DataContainer<Qn::StatCollect, Qn::AxisD>();
  dc_uQ_perp->AddAxes({{"impactpar", C_edges}});
  dc_uQ_perp->AddAxes(dc_usim->GetAxes());
  for(int iBin=0; iBin<dc_uQ_perp->size(); iBin++)
    dc_uQ_perp->At(iBin).SetNumberOfSamples(nsam);

  const int nEvents = treeIn->GetEntries();
//   const int nEvents = 10000;

  Qn::Correlation::ReSampleFunctor rsf(nsam);

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

    dc_Q_par->At(C_bin).Fill(Qpar, 1.0, rsf());
    dc_Q_perp->At(C_bin).Fill(Qperp, 1.0, rsf());

    for(int iBin=0; iBin<dc_usim->size(); iBin++) {
      std::vector<size_t> indices = dc_usim->GetIndex(iBin);
      MyQVec qsim = MyQVec(dc_usim->At(iBin));

      const float upar = ScalarProduct(qsim, qpsi);
      const float uperp = ScalarProduct(qsim, qpsiR);

      dc_u_par->At({C_bin, indices.at(0), indices.at(1)}).Fill(upar, qsim.sumw_, rsf());
      dc_u_perp->At({C_bin, indices.at(0), indices.at(1)}).Fill(uperp, qsim.sumw_, rsf());
      dc_uQ_par->At({C_bin, indices.at(0), indices.at(1)}).Fill(upar*Qpar, qsim.sumw_, rsf());
      dc_uQ_perp->At({C_bin, indices.at(0), indices.at(1)}).Fill(uperp*Qperp, qsim.sumw_, rsf());
    }
  }

  TFile* fileOut = TFile::Open(("covariances_scol." + particle + ".root").c_str(), "recreate");
  fileOut->cd();
  dc_u_par->Write("upar.all");
  dc_Q_par->Write("Qpar.all");
  dc_uQ_par->Write("uQpar.all");
  dc_u_perp->Write("uperp.all");
  dc_Q_perp->Write("Qperp.all");
  dc_uQ_perp->Write("uQperp.all");
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
  const float sp = v1.x_*v2.x_ + v1.y_*v2.y_;
  return sp;
}
