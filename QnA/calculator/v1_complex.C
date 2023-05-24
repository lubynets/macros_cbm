void v1_complex()
{
  std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/cl.dcmqgsm.12agev.lc1.3122.root";

  TFile* fileIn = TFile::Open(fileName.c_str());
  
  std::vector<std::string> components{"x1x1", "y1y1"/*, "x1y1", "y1x1"*/};
  std::vector<std::string> subevents{"psd1", "psd2", "psd3"};
  std::vector<std::string> sub4th{"sts_pipos"};
  
  TFile* fileOut = TFile::Open("v1andR1.root", "recreate");
  fileOut->cd();
  fileOut->mkdir("v1/usimPsi");
  fileOut->mkdir("v1/uPsi");
  fileOut->mkdir("R1/MC");
  fileOut->mkdir("R1/sub3");
  fileOut->mkdir("R1/sub4");
  fileOut->mkdir("v1/uQ_R1_MC");
  fileOut->mkdir("v1/uQ_R1_sub3");
  fileOut->mkdir("v1/uQ_R1_sub4");
  
  Correlation usimPsi(fileIn, "uPsi", {"u_sim_PLAIN", "Q_psi_PLAIN"}, components);
  usimPsi = usimPsi*2.;
  fileOut->cd("v1/usimPsi");
  usimPsi.Save("v1.u_sim.Q_psi");
  
  Correlation uPsi(fileIn, "uPsi", {"u_rec_RESCALED", "Q_psi_PLAIN"}, components);
  uPsi = uPsi*2.;
  fileOut->cd("v1/uPsi");
  uPsi.Save("v1.u_rec.Q_psi");
  
  Correlation usgnlPsi(fileIn, "uPsi", {"u_rec_sgnl_RESCALED", "Q_psi_PLAIN"}, components);
  usgnlPsi = usgnlPsi*2.;
  fileOut->cd("v1/uPsi");
  usgnlPsi.Save("v1.u_rec_sgnl.Q_psi");
  
  std::vector<Correlation> R1_MC;
  R1_MC.resize(subevents.size());
  fileOut->cd("R1/MC");
  for(int i=0; i<subevents.size(); i++) {
    R1_MC.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i) + "_RECENTERED", "Q_psi_PLAIN"}, components) * 2.;
    R1_MC.at(i).Save("res_MC." + subevents.at(i));
  }

  std::vector<Correlation> QQ;
  QQ.resize(subevents.size());
  for(int i=0; i<subevents.size(); i++) {
    QQ.at(i) = Correlation(fileIn, "QQ", {subevents.at((i+1)%3) + "_RECENTERED", subevents.at((i+2)%3) + "_RECENTERED"}, components) * 2.;
  }
  
  std::vector<Correlation> R1_sub3;
  R1_sub3.resize(subevents.size());
  fileOut->cd("R1/sub3");
  R1_sub3 = Functions::VectorResolutions3S(QQ.at(0), QQ.at(1), QQ.at(2));
  for(int i=0; i<subevents.size(); i++) {
    R1_sub3.at(i).Save("res_sub3." + subevents.at(i));
  }

  std::vector<std::vector<Correlation>> R1_sub4;
  R1_sub4.resize(sub4th.size());
  for(auto& r : R1_sub4) {
    r.resize(subevents.size());
  }
  fileOut->cd("R1/sub4");
  for(int j=0; j<sub4th.size(); j++) {
    auto s4 = sub4th.at(j);
    Correlation Q13 = Correlation(fileIn, "QQ", {subevents.at(0) + "_RECENTERED", subevents.at(2) + "_RECENTERED"}, components) * 2.;
    Correlation Q14 = Correlation(fileIn, "QQ", {subevents.at(0) + "_RECENTERED", s4 + "_RESCALED"}, components) * 2.;
    Correlation Q24 = Correlation(fileIn, "QQ", {subevents.at(1) + "_RECENTERED", s4 + "_RESCALED"}, components) * 2.;
    Correlation Q34 = Correlation(fileIn, "QQ", {subevents.at(2) + "_RECENTERED", s4 + "_RESCALED"}, components) * 2.;

    R1_sub4.at(j) = Functions::VectorResolutions4S(Q14, Q24, Q34, Q13);
    for(int i=0; i<subevents.size(); i++) {
      R1_sub4.at(j).at(i).Save("res.sub4." + s4 + "." + subevents.at(i));
    }
  }
  
  for(int i=0; i<subevents.size(); i++) {
    Correlation uQ(fileIn, "uQ", {"u_rec_RESCALED", subevents.at(i) + "_RECENTERED"}, components);

    fileOut->cd("v1/uQ_R1_MC");
    auto uQ_R1_MC = uQ / R1_MC.at(i) * 2.;
    uQ_R1_MC.Save("v1.u_rec." + subevents.at(i) + ".res_MC");

    fileOut->cd("v1/uQ_R1_sub3");
    auto uQ_R1_sub3 = uQ / R1_sub3.at(i) * 2.;
    uQ_R1_sub3.Save("v1.u_rec." + subevents.at(i) + ".res_sub3");

    fileOut->cd("v1/uQ_R1_sub4");
    for(int j=0; j<sub4th.size(); j++) {
      auto uQ_R1_sub4 = uQ / R1_sub4.at(j).at(i) * 2.;
      uQ_R1_sub4.Save("v1.u_rec." + subevents.at(i)+ ".res_sub4_" + sub4th.at(j));
    }
  }

  for(int i=0; i<subevents.size(); i++) {
    Correlation uQ(fileIn, "uQ", {"u_rec_sgnl_RESCALED", subevents.at(i) + "_RECENTERED"}, components);

    fileOut->cd("v1/uQ_R1_MC");
    auto uQ_R1_MC = uQ / R1_MC.at(i) * 2.;
    uQ_R1_MC.Save("v1.u_rec_sgnl." + subevents.at(i) + ".res_MC");

    fileOut->cd("v1/uQ_R1_sub3");
    auto uQ_R1_sub3 = uQ / R1_sub3.at(i) * 2.;
    uQ_R1_sub3.Save("v1.u_rec_sgnl." + subevents.at(i) + ".res_sub3");

    fileOut->cd("v1/uQ_R1_sub4");
    for(int j=0; j<sub4th.size(); j++) {
      auto uQ_R1_sub4 = uQ / R1_sub4.at(j).at(i) * 2.;
      uQ_R1_sub4.Save("v1.u_rec_sgnl." + subevents.at(i)+ ".res_sub4_" + sub4th.at(j));
    }
  }

  fileOut->Close();
}
