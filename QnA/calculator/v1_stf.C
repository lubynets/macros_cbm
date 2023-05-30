void v1_stf() {
//   std::string evegen = "dcmqgsm";
  std::string evegen = "urqmd";

  std::string pbeam = "12";
//   std::string pbeam = "3.3";

//   std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/" + pbeam + "agev/cl.stf." + evegen + "." + pbeam + "agev.root";
  std::string fileName = "/home/oleksii/cbmdir/sandbox/cl.d3.root";

  TFile* fileIn = TFile::Open(fileName.c_str());
  
  std::vector<std::string> particles{/*"lambda", "kshort", "xi", "pipos", "pineg"*/};
  std::vector<std::string> components{"x1x1", "y1y1", "x1y1", "y1x1"};
  std::vector<std::string> components_same{"x1x1", "y1y1"};
  std::vector<std::string> components_cross{"x1y1", "y1x1"};
  std::vector<std::string> subevents{"psd1_RECENTERED",
                                     "psd2_RECENTERED",
                                     "psd3_RECENTERED",
//                                      "etacut_1_charged_PLAIN",
//                                      "etacut_2_charged_PLAIN",
//                                      "etacut_3_charged_PLAIN",
//                                      "etacut_1_all_PLAIN",
//                                      "etacut_2_all_PLAIN",
//                                      "etacut_3_all_PLAIN"
                                     };
  std::vector<std::string> sub4th{"sts_pipos_yS_nocut", "sts_pipos_yS_cut", "sts_pipos_yL_nocut", "sts_pipos_yL_cut"};

  std::vector<std::pair<int, int>> sub_indices{{0, 3}/*, {3, 6}, {6, 9}*/};
  
  TFile* fileOut = TFile::Open(("v1andR1.stf." + evegen + "." + pbeam + "agev.root").c_str(), "recreate");
  fileOut->cd();
  for(auto& pa : particles) {
    fileOut->mkdir(("v1/" + pa + "/uPsi").c_str());
    fileOut->mkdir(("v1/" + pa + "/uQ_R1").c_str());
    fileOut->mkdir(("v1/" + pa + "/uQ").c_str());
  }
  fileOut->mkdir("R1");
  fileOut->mkdir("R1R1");
  fileOut->mkdir("QQ");

  for(auto& pa : particles) {
    Correlation uPsi(fileIn, "uPsi", {("u_sim_" + pa + "_PLAIN").c_str(), "Q_psi_PLAIN"}, components);
    uPsi = uPsi*2.;
    fileOut->cd(("v1/" + pa + "/uPsi").c_str());
    uPsi.Save("v1.uPsi");
  }
  
  fileOut->cd("R1");
  std::vector<Correlation> R1;
  R1.resize(subevents.size() + sub4th.size());
  for(int i=0; i<subevents.size(); i++) {
    R1.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i), "Q_psi_PLAIN"}, components_same) * 2.;
    R1.at(i).Save("res.mc." + subevents.at(i));
  }
  for(int i=0; i<sub4th.size(); i++) {
    R1.at(subevents.size() + i) = Correlation(fileIn, "QPsi", {sub4th.at(i) + "_RESCALED", "Q_psi_PLAIN"}, components_same) * 2.;
    R1.at(subevents.size() + i).Save("res.mc." + sub4th.at(i) + "_RESCALED");
  }

  fileOut->cd("R1R1");
  for(auto& si : sub_indices) {
    for(int i=si.first; i<si.second; i++) {
      for(int j=si.first; j<si.second; j++) {
        auto R1R1 = R1.at(i) * R1.at(j);
        R1R1.Save("r1r1." + subevents.at(i) + "_" + subevents.at(j));
      }
    }
  }

  for(int j=0; j<sub4th.size(); j++){
    for(int i=0; i<3; i++){
      auto R1R1 = R1.at(i) * R1.at(subevents.size() + j);
      fileOut->cd("R1R1");
      R1R1.Save("r1r1." + subevents.at(i) + "_" + sub4th.at(j));

      fileOut->cd("QQ");
      Correlation QQ = Correlation(fileIn, "QQ", {subevents.at(i), sub4th.at(j) + "_RESCALED"}, components_same) * 2.;
      QQ.Save("qq." + subevents.at(i) + "_" + sub4th.at(j));
    }
  }

  for(auto& si : sub_indices) {
    std::vector<std::vector<Correlation>> QQ;
    QQ.resize(si.second - si.first);
    for(auto& qq : QQ) {
      qq.resize(si.second - si.first);
    }
    fileOut->cd("QQ");
    for(int i=si.first; i<si.second; i++) {
      for(int j=si.first; j<si.second; j++) {
        if(i>j) continue;
        int I = i - si.first;
        int J = j - si.first;
        QQ.at(I).at(J) = Correlation(fileIn, "QQ", {subevents.at(i), subevents.at(j)}, components_same) * 2.;
        QQ.at(I).at(J).Save("qq." + subevents.at(i) + "_" + subevents.at(j));
      }
    }
    fileOut->cd("R1");
    std::vector<Correlation> R1_3sub = Functions::VectorResolutions3S(QQ.at(1).at(2), QQ.at(0).at(2), QQ.at(0).at(1));
    for(int i=0; i<3; i++) {
      R1_3sub.at(i).Save("res.sub3." + subevents.at(si.first + i));
    }
  }

  for(auto& s4 : sub4th){
    Correlation Q13 = Correlation(fileIn, "QQ", {subevents.at(0), subevents.at(2)}, components_same) * 2.;
    Correlation Q14 = Correlation(fileIn, "QQ", {subevents.at(0), (s4 + "_RESCALED").c_str()}, components_same) * 2.;
    Correlation Q24 = Correlation(fileIn, "QQ", {subevents.at(1), (s4 + "_RESCALED").c_str()}, components_same) * 2.;
    Correlation Q34 = Correlation(fileIn, "QQ", {subevents.at(2), (s4 + "_RESCALED").c_str()}, components_same) * 2.;

    std::vector<Correlation> R1_4sub = Functions::VectorResolutions4S(Q14, Q24, Q34, Q13);
    for(int i=0; i<3; i++) {
      R1_4sub.at(i).Save("res.sub4." + s4 + "." + subevents.at(i));
    }
  }

  std::vector<Correlation> R1_cross;
  R1_cross.resize(subevents.size());
  fileOut->cd("R1");
  for(int i=0; i<subevents.size(); i++) {
    R1_cross.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i), "Q_psi_PLAIN"}, components_cross) * 2.;
    R1_cross.at(i).Save("res_cross." + subevents.at(i));
  }
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      R1_cross.at(i) = Correlation(fileIn, "QQ", {subevents.at(i), subevents.at(j)}, components_cross) * 2.;
      R1_cross.at(i).Save("res_cross." + subevents.at(i) + "." + subevents.at(j));
    }
  }

  for(auto& pa : particles) {
    for(int i=0; i<subevents.size(); i++) {
      Correlation uQ(fileIn, "uQ", {("u_sim_" + pa + "_PLAIN").c_str(), subevents.at(i)}, components_same);
      auto uQ_R1 = uQ / R1.at(i) * 2.;
      fileOut->cd(("v1/" + pa + "/uQ_R1").c_str());
      uQ_R1.Save("v1.uQ_R1." + subevents.at(i));
      auto uQ_R1_mirrored = Mirror(uQ_R1, "SimParticles_rapidity");
      auto uQ_R1_even = (uQ_R1 - uQ_R1_mirrored) / 2.;
      uQ_R1_even.Save("v1.uQ_R1_even." + subevents.at(i));
      auto uQ_R1_odd = (uQ_R1 + uQ_R1_mirrored) / 2.;
      uQ_R1_odd.Save("v1.uQ_R1_odd." + subevents.at(i));
    }
  }

  for(auto& pa : particles) {
    for(int i=0; i<subevents.size(); i++) {
      Correlation uQ(fileIn, "uQ", {("u_sim_" + pa + "_PLAIN").c_str(), subevents.at(i)}, components_cross);
      fileOut->cd(("v1/" + pa + "/uQ").c_str());
      uQ.Save("v1.uQ." + subevents.at(i));
    }
  }

  fileOut->Close();
}
