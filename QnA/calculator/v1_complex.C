void v1_complex()
{
//   std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/cl.dcmqgsm.apr20.recpid.lightcuts1.3122.set4.sgnl.root";
  std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/cl.dcmqgsm.apr20.recpid.lightcuts1.310.set4.root";

  TFile* fileIn = TFile::Open(fileName.c_str());
  
  std::vector<std::string> components{"x1x1", "y1y1", "x1y1", "y1x1"};
  std::vector<std::string> subevents{"psd1_RECENTERED", "psd2_RECENTERED", "psd3_RECENTERED"};
  
  TFile* fileOut = TFile::Open("v1andR1.root", "recreate");
  fileOut->cd();
  fileOut->mkdir("v1/usimPsi");
  fileOut->mkdir("v1/uPsi");
  fileOut->mkdir("R1/MC");
  fileOut->mkdir("R1/sub3");
  fileOut->mkdir("R1/sub4_sts_p");
  fileOut->mkdir("R1/sub4_sts_pipos");
  fileOut->mkdir("v1/uQ_R1_MC");
  fileOut->mkdir("v1/uQ_R1_sub3");
  fileOut->mkdir("v1/uQ_R1_sub4_sts_p");
  fileOut->mkdir("v1/uQ_R1_sub4_sts_pipos");
  
  Correlation usimPsi(fileIn, "uPsi", {"u_sim_PLAIN", "Q_psi_PLAIN"}, components);
  usimPsi = usimPsi*2.;
  fileOut->cd("v1/usimPsi");
  usimPsi.Save("v1.u_sim_PLAIN.Q_psi_PLAIN");
  
  Correlation uPsi(fileIn, "uPsi", {"u_rec_RESCALED", "Q_psi_PLAIN"}, components);
  uPsi = uPsi*2.;
  fileOut->cd("v1/uPsi");
  uPsi.Save("v1.u_rec_RESCALED.Q_psi_PLAIN");
  
  Correlation usgnlPsi(fileIn, "uPsi", {"u_rec_sgnl_RESCALED", "Q_psi_PLAIN"}, components);
  usgnlPsi = usgnlPsi*2.;
  fileOut->cd("v1/uPsi");
  usgnlPsi.Save("v1.u_rec_sgnl_RESCALED.Q_psi_PLAIN");
  
  std::vector<Correlation> R1_MC;
  R1_MC.resize(subevents.size());  
  fileOut->cd("R1/MC");
  for(int i=0; i<subevents.size(); i++) {
    R1_MC.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i), "Q_psi_PLAIN"}, components) * 2.;
    R1_MC.at(i).Save("res_MC." + subevents.at(i));
  }
  
  std::vector<Correlation> R1_sub3;
  R1_sub3.resize(subevents.size());  
  fileOut->cd("R1/sub3");
  for(int i=0; i<subevents.size(); i++)
  {
    R1_sub3.at(i) = Functions::VectorResolutions3S( fileIn, "QQ", subevents.at(i), subevents, components ).at(0);
    R1_sub3.at(i).Save("res_sub3." + subevents.at(i));
  }

  std::vector<Correlation> R1_sub4_sts_p;
  R1_sub4_sts_p.resize(subevents.size());
  fileOut->cd("R1/sub4_sts_p");
  
  R1_sub4_sts_p.at(0) = Functions::VectorResolutions3S( fileIn, "QQ", "psd1_RECENTERED", {"sts_p_RESCALED", "psd3_RECENTERED"}, components ).at(0);
  R1_sub4_sts_p.at(0).Save("res_sub4_sts_p." + subevents.at(0));
  
  auto res_sts_p = Functions::VectorResolutions3S( fileIn, "QQ", "sts_p_RESCALED", {"psd1_RECENTERED", "psd3_RECENTERED"}, components ).at(0);
  Correlation psd2_sts_p( fileIn, "QQ", {"psd2_RECENTERED", "sts_p_RESCALED"}, components);
  R1_sub4_sts_p.at(1) = psd2_sts_p/res_sts_p*(-2.);
  R1_sub4_sts_p.at(1).Save("res_sub4_sts_p." + subevents.at(1));  
  
  R1_sub4_sts_p.at(2) = Functions::VectorResolutions3S( fileIn, "QQ", "psd3_RECENTERED", {"sts_p_RESCALED", "psd1_RECENTERED"}, components ).at(0);
  R1_sub4_sts_p.at(2).Save("res_sub4_sts_p." + subevents.at(2));
  
  std::vector<Correlation> R1_sub4_sts_pipos;
  R1_sub4_sts_pipos.resize(subevents.size());
  fileOut->cd("R1/sub4_sts_pipos");
  
  R1_sub4_sts_pipos.at(0) = Functions::VectorResolutions3S( fileIn, "QQ", "psd1_RECENTERED", {"sts_pipos_RESCALED", "psd3_RECENTERED"}, components ).at(0);
  R1_sub4_sts_pipos.at(0).Save("res_sub4_sts_pipos." + subevents.at(0));
  
  auto res_sts_pipos = Functions::VectorResolutions3S( fileIn, "QQ", "sts_pipos_RESCALED", {"psd1_RECENTERED", "psd3_RECENTERED"}, components ).at(0);
  Correlation psd2_sts_pipos( fileIn, "QQ", {"psd2_RECENTERED", "sts_pipos_RESCALED"}, components);
  R1_sub4_sts_pipos.at(1) = psd2_sts_pipos/res_sts_pipos*(-2.);
  R1_sub4_sts_pipos.at(1).Save("res_sub4_sts_pipos." + subevents.at(1));  
  
  R1_sub4_sts_pipos.at(2) = Functions::VectorResolutions3S( fileIn, "QQ", "psd3_RECENTERED", {"sts_pipos_RESCALED", "psd1_RECENTERED"}, components ).at(0);
  R1_sub4_sts_pipos.at(2).Save("res_sub4_sts_pipos." + subevents.at(2));
  
  for(int i=0; i<subevents.size(); i++) {
    Correlation uQ(fileIn, "uQ", {"u_rec_RESCALED", subevents.at(i)}, components);
    auto uQ_R1_MC = uQ / R1_MC.at(i) * 2.;
    auto uQ_R1_sub3 = uQ / R1_sub3.at(i) * 2.;
    auto uQ_R1_sub4_sts_p = uQ / R1_sub4_sts_p.at(i) * 2.;
    auto uQ_R1_sub4_sts_pipos = uQ / R1_sub4_sts_pipos.at(i) * 2.;
        
    fileOut->cd("v1/uQ_R1_MC");
    uQ_R1_MC.Save("v1.u_rec_RESCALED." + subevents.at(i) + ".res_MC");
    fileOut->cd("v1/uQ_R1_sub3");
    uQ_R1_sub3.Save("v1.u_rec_RESCALED." + subevents.at(i) + ".res_sub3");
    fileOut->cd("v1/uQ_R1_sub4_sts_p");
    uQ_R1_sub4_sts_p.Save("v1.u_rec_RESCALED." + subevents.at(i) + ".res_sub4_sts_p");
    fileOut->cd("v1/uQ_R1_sub4_sts_pipos");
    uQ_R1_sub4_sts_pipos.Save("v1.u_rec_RESCALED." + subevents.at(i) + ".res_sub4_sts_pipos");
  }
  
  for(int i=0; i<subevents.size(); i++) {
    Correlation uQ(fileIn, "uQ", {"u_rec_sgnl_RESCALED", subevents.at(i)}, components);
    auto uQ_R1_MC = uQ / R1_MC.at(i) * 2.;
    auto uQ_R1_sub3 = uQ / R1_sub3.at(i) * 2.;
    auto uQ_R1_sub4_sts_p = uQ / R1_sub4_sts_p.at(i) * 2.;
    auto uQ_R1_sub4_sts_pipos = uQ / R1_sub4_sts_pipos.at(i) * 2.;
        
    fileOut->cd("v1/uQ_R1_MC");
    uQ_R1_MC.Save("v1.u_rec_sgnl_RESCALED." + subevents.at(i) + ".res_MC");
    fileOut->cd("v1/uQ_R1_sub3");
    uQ_R1_sub3.Save("v1.u_rec_sgnl_RESCALED." + subevents.at(i) + ".res_sub3");
    fileOut->cd("v1/uQ_R1_sub4_sts_p");
    uQ_R1_sub4_sts_p.Save("v1.u_rec_sgnl_RESCALED." + subevents.at(i) + ".res_sub4_sts_p");
    fileOut->cd("v1/uQ_R1_sub4_sts_pipos");
    uQ_R1_sub4_sts_pipos.Save("v1.u_rec_sgnl_RESCALED." + subevents.at(i) + ".res_sub4_sts_pipos");
  }
  
  fileOut->Close();
}
