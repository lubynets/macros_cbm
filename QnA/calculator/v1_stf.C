void v1_stf() {
//   std::string evegen = "dcmqgsm";
  std::string evegen = "urqmd";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/cl.stf." + evegen + ".root";
  
  TFile* fileIn = TFile::Open(fileName.c_str());
  
  std::vector<std::string> particles{"lambda", "kshort", "pipos", "pineg"};
  std::vector<std::string> components{"x1x1", "y1y1", "x1y1", "y1x1"};
  std::vector<std::string> components_same{"x1x1", "y1y1"};
  std::vector<std::string> components_cross{"x1y1", "y1x1"};
  std::vector<std::string> subevents{"psd1_RECENTERED",
                                     "psd2_RECENTERED",
                                     "psd3_RECENTERED",
//                                      "psdall_RECENTERED",
                                     "spec1_prim_PLAIN",
                                     "spec2_prim_PLAIN",
                                     "spec3_prim_PLAIN",
//                                      "spec1_all_PLAIN",
//                                      "spec2_all_PLAIN",
//                                      "spec3_all_PLAIN"
                                    };
  
  TFile* fileOut = TFile::Open(("v1andR1.stf." + evegen + ".root").c_str(), "recreate");
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
  
  std::vector<Correlation> R1;
  R1.resize(subevents.size());
  fileOut->cd("R1");
  for(int i=0; i<subevents.size(); i++) {
    R1.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i), "Q_psi_PLAIN"}, components_same) * 2.;
    R1.at(i).Save("res." + subevents.at(i));
  }

  fileOut->cd("R1R1");
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      auto R1R1 = R1.at(i) * R1.at(j);
      R1R1.Save("r1r1." + subevents.at(i) + "_" + subevents.at(j));
    }
  }
  for(int i=3; i<6; i++) {
    for(int j=3; j<6; j++) {
      auto R1R1 = R1.at(i) * R1.at(j);
      R1R1.Save("r1r1." + subevents.at(i) + "_" + subevents.at(j));
    }
  }

  fileOut->cd("QQ");
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      auto QQ = Correlation(fileIn, "QQ", {subevents.at(i), subevents.at(j)}, components_same) * 2.;
      QQ.Save("qq." + subevents.at(i) + "_" + subevents.at(j));
    }
  }
  for(int i=3; i<6; i++) {
    for(int j=3; j<6; j++) {
      auto QQ = Correlation(fileIn, "QQ", {subevents.at(i), subevents.at(j)}, components_same) * 2.;
      QQ.Save("qq." + subevents.at(i) + "_" + subevents.at(j));
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
//       uQ_R1_mirrored.Save("v1.uQ_R1_mirrored." + subevents.at(i));
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
