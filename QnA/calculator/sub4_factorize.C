void sub4_factorize() {
  std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/cl.dcmqgsm.12agev.lc1.3122.root";

  TFile* fileIn = TFile::Open(fileName.c_str());

  std::vector<std::string> components{"x1x1", "y1y1"};
  std::vector<std::string> subevents{"psd1_RECENTERED",
                                     "psd2_RECENTERED",
                                     "psd3_RECENTERED"};
  std::string sub4th{"sts_pipos"};
  std::string sub4th_midrap{"sts_pipos_midrap"};

  TFile* fileOut = TFile::Open(("v1andR1.stf." + evegen + "." + pbeam + "agev.root").c_str(), "recreate");
  fileOut->cd();

  fileOut->mkdir("R1");
  fileOut->mkdir("R1R1");
  fileOut->mkdir("QQ");

  fileOut->cd("R1");
  std::vector<Correlation> R1;
  R1.resize(subevents.size() + 1);
  for(int i=0; i<subevents.size(); i++) {
    R1.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i), "Q_psi_PLAIN"}, components) * 2.;
    R1.at(i).Save("res.mc." + subevents.at(i));
  }
  R1.at(subevents.size()) = Correlation(fileIn, "QPsi", {sub4th + "_RESCALED", "Q_psi_PLAIN"}, components) * 2.;
  R1.at(subevents.size()).Save("res.mc." + sub4th.at(i) + "_RESCALED");

  fileOut->cd("R1R1");
  for(int i=0; i<subevents.size(); i++) {
    for(int j=0; j<subevents.size(); j++) {
      if(i == j) continue;
      auto R1R1 = R1.at(i) * R1.at(j);
      R1R1.Save("r1r1." + subevents.at(i) + "_" + subevents.at(j));
    }
    auto R1R1 = R1.at(i) * R1.at(subevents.size());
    R1R1.Save("r1r1." + subevents.at(i) + "_" + sub4th);
  }

  Correlation Q13 = Correlation(fileIn, "QQ", {subevents.at(0), subevents.at(2)}, components) * 2.;

  Correlation Q14 = Correlation(fileIn, "QQ", {subevents.at(0), (sub4th + "_RESCALED").c_str()}, components) * 2.;
  Correlation Q24 = Correlation(fileIn, "QQ", {subevents.at(1), (sub4th + "_RESCALED").c_str()}, components) * 2.;
  Correlation Q34 = Correlation(fileIn, "QQ", {subevents.at(2), (sub4th + "_RESCALED").c_str()}, components) * 2.;

  Correlation Q14_mr = Correlation(fileIn, "QQ", {subevents.at(0), (sub4th_midrap + "_RESCALED").c_str()}, components) * 2.;
  Correlation Q24_mr = Correlation(fileIn, "QQ", {subevents.at(1), (sub4th_midrap + "_RESCALED").c_str()}, components) * 2.;
  Correlation Q34_mr = Correlation(fileIn, "QQ", {subevents.at(2), (sub4th_midrap + "_RESCALED").c_str()}, components) * 2.;

  Correlation Q14_corr = Q14 - Q14_mr;
  Correlation Q24_corr = Q24 - Q24_mr;
  Correlation Q34_corr = Q34 - Q34_mr;

  fileOut->cd("QQ");
  Q14.Save("qq." + subevents.at(0) + "_" + sub4th);
  Q24.Save("qq." + subevents.at(1) + "_" + sub4th);
  Q34.Save("qq." + subevents.at(2) + "_" + sub4th);

  Q14_mr.Save("qq_mr." + subevents.at(0) + "_" + sub4th);
  Q24_mr.Save("qq_mr." + subevents.at(1) + "_" + sub4th);
  Q34_mr.Save("qq_mr." + subevents.at(2) + "_" + sub4th);

  Q14_corr.Save("qq_corr." + subevents.at(0) + "_" + sub4th);
  Q24_corr.Save("qq_corr." + subevents.at(1) + "_" + sub4th);
  Q34_corr.Save("qq_corr." + subevents.at(2) + "_" + sub4th);

  std::vector<Correlation> R1_4sub = Functions::VectorResolutions4S(Q14, Q24, Q34, Q13);
    for(int i=0; i<3; i++) {
      R1_4sub.at(i).Save("res.sub4." + sub4th + "." + subevents.at(i));
  }

}
