void v2_complex()
{
//   std::string fileName = "/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.12agev.defcuts.3122.set2.sgnl_1.root";
//   std::string fileName = "/home/user/cbmdir/working/qna/correlations/aXmass/cl.dcmqgsm.apr20.recpid.defaultcuts.3312.set4.root";
  std::string fileName = "/home/user/cbmdir/working/qna/correlations/aXmass/cl.dcmqgsm.apr20.recpid.lightcuts1.3122.set4_v2.root";
  
  TFile* fileIn = TFile::Open(fileName.c_str());
  
  std::vector<std::string> components_1{"x1x1", "y1y1"};
  std::vector<std::string> components_2{"x2x2", "y2y2"};
  std::vector<std::string> components_211{"x2x1x1", "x2y1y1", "y2x1y1", "y2y1x1"};
  std::vector<std::string> subevents{"psd1_RECENTERED", "psd2_RECENTERED", "psd3_RECENTERED"};
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  fileOut->cd();
  fileOut->mkdir("v2/usimPsi");
  fileOut->mkdir("v2/uPsi");
  fileOut->mkdir("R1/MC");
  fileOut->mkdir("R1/sub3");
  fileOut->mkdir("R1/sub4_sts_p");
  fileOut->mkdir("R1/sub4_sts_pipos");
  fileOut->mkdir("v2/uQQ_R1_MC");
  fileOut->mkdir("v2/uQQ_R1_sub3");
  fileOut->mkdir("v2/uQQ_R1_sub4_sts_p");
  fileOut->mkdir("v2/uQQ_R1_sub4_sts_pipos");
  
  Correlation usimPsi(fileIn, "uPsi", {"u_sim_PLAIN", "Q_psi_PLAIN"}, components_2);
  usimPsi = usimPsi*2.;
  fileOut->cd("v2/usimPsi");
  usimPsi.Save("v2.u_sim_PLAIN.Q_psi_PLAIN");
  
//   Correlation uPsi(fileIn, "uPsi", {"u_rec_RESCALED", "Q_psi_PLAIN"}, components);
//   uPsi = uPsi*2.;
//   fileOut->cd("v1/uPsi");
//   uPsi.Save("v1.u_rec_RESCALED.Q_psi_PLAIN");
//   
//   Correlation usgnlPsi(fileIn, "uPsi", {"u_rec_sgnl_RESCALED", "Q_psi_PLAIN"}, components);
//   usgnlPsi = usgnlPsi*2.;
//   fileOut->cd("v1/uPsi");
//   usgnlPsi.Save("v1.u_rec_sgnl_RESCALED.Q_psi_PLAIN");
  
//   std::vector<Correlation> R1_MC;
//   R1_MC.resize(subevents.size());  
//   fileOut->cd("R1/MC");
//   for(int i=0; i<subevents.size(); i++) {
//     R1_MC.at(i) = Correlation(fileIn, "QPsi", {subevents.at(i), "Q_psi_PLAIN"}, components_1) * 2.;
//     R1_MC.at(i).Save("res_MC." + subevents.at(i));
//   }
  
  std::vector<Correlation> R1_sub3;
  R1_sub3.resize(subevents.size());  
  fileOut->cd("R1/sub3");
  for(int i=0; i<subevents.size(); i++)
  {
    R1_sub3.at(i) = Functions::VectorResolutions3S( fileIn, "QQ", subevents.at(i), subevents, components_1 ).at(0);
    R1_sub3.at(i).Save("res_sub3." + subevents.at(i));
  }

  std::vector<Correlation> R1_sub4_sts_p;
  R1_sub4_sts_p.resize(subevents.size());
  fileOut->cd("R1/sub4_sts_p");
  
  R1_sub4_sts_p.at(0) = Functions::VectorResolutions3S( fileIn, "QQ", "psd1_RECENTERED", {"sts_p_RESCALED", "psd3_RECENTERED"}, components_1 ).at(0);
  R1_sub4_sts_p.at(0).Save("res_sub4_sts_p." + subevents.at(0));
  
  auto res_sts_p = Functions::VectorResolutions3S( fileIn, "QQ", "sts_p_RESCALED", {"psd1_RECENTERED", "psd3_RECENTERED"}, components_1 ).at(0);
  Correlation psd2_sts_p( fileIn, "QQ", {"psd2_RECENTERED", "sts_p_RESCALED"}, components_1);
  R1_sub4_sts_p.at(1) = psd2_sts_p/res_sts_p*(-2.);
  R1_sub4_sts_p.at(1).Save("res_sub4_sts_p." + subevents.at(1));  
  
  R1_sub4_sts_p.at(2) = Functions::VectorResolutions3S( fileIn, "QQ", "psd3_RECENTERED", {"sts_p_RESCALED", "psd1_RECENTERED"}, components_1 ).at(0);
  R1_sub4_sts_p.at(2).Save("res_sub4_sts_p." + subevents.at(2));
  
  std::vector<Correlation> R1_sub4_sts_pipos;
  R1_sub4_sts_pipos.resize(subevents.size());
  fileOut->cd("R1/sub4_sts_pipos");
  
  R1_sub4_sts_pipos.at(0) = Functions::VectorResolutions3S( fileIn, "QQ", "psd1_RECENTERED", {"sts_pipos_RESCALED", "psd3_RECENTERED"}, components_1 ).at(0);
  R1_sub4_sts_pipos.at(0).Save("res_sub4_sts_pipos." + subevents.at(0));
  
  auto res_sts_pipos = Functions::VectorResolutions3S( fileIn, "QQ", "sts_pipos_RESCALED", {"psd1_RECENTERED", "psd3_RECENTERED"}, components_1 ).at(0);
  Correlation psd2_sts_pipos( fileIn, "QQ", {"psd2_RECENTERED", "sts_pipos_RESCALED"}, components_1);
  R1_sub4_sts_pipos.at(1) = psd2_sts_pipos/res_sts_pipos*(-2.);
  R1_sub4_sts_pipos.at(1).Save("res_sub4_sts_pipos." + subevents.at(1));  
  
  R1_sub4_sts_pipos.at(2) = Functions::VectorResolutions3S( fileIn, "QQ", "psd3_RECENTERED", {"sts_pipos_RESCALED", "psd1_RECENTERED"}, components_1 ).at(0);
  R1_sub4_sts_pipos.at(2).Save("res_sub4_sts_pipos." + subevents.at(2));
  
  for(int i=0; i<subevents.size(); i++) {
    for(int j=0; j<subevents.size(); j++) {
      Correlation uQQ(fileIn, "uQ", {"u_rec_RESCALED", subevents.at(i), subevents.at(j)}, components_211);
      
      for(int k=0; k<3; k++) {
        
        const int ki = k%2;
        int kj = 0;
        if(k==1 || k==2)
          kj = 1;
        
        auto uQ_R1_sub3 = uQQ[k] / R1_sub3.at(i)[ki] / R1_sub3.at(j)[kj] * 4.;
        auto uQ_R1_sub4_sts_p = uQQ[k] / R1_sub4_sts_p.at(i)[ki] / R1_sub4_sts_p.at(j)[kj] * 4.;
        auto uQ_R1_sub4_sts_pipos = uQQ[k] / R1_sub4_sts_pipos.at(i)[ki] / R1_sub4_sts_pipos.at(j)[kj] * 4.;
        
        fileOut->cd("v2/uQQ_R1_sub3");
        uQ_R1_sub3.Write(("v2.u_rec_RESCALED." + subevents.at(i) + "." + subevents.at(j) + ".res_sub3" + "." + components_211[k]).c_str());
        fileOut->cd("v2/uQQ_R1_sub4_sts_p");
        uQ_R1_sub4_sts_p.Write(("v2.u_rec_RESCALED." + subevents.at(i) + "." + subevents.at(j) + ".res_sub4_sts_p" + "." + components_211[k]).c_str());
        fileOut->cd("v2/uQQ_R1_sub4_sts_pipos");
        uQ_R1_sub4_sts_pipos.Write(("v2.u_rec_RESCALED." + subevents.at(i) + "." + subevents.at(j) + ".res_sub4_sts_pipos" + "." + components_211[k]).c_str());
      }
    }
  }
  
//   for(int i=0; i<subevents.size(); i++) {
//     for(int j=0; j<subevents.size(); j++) {
//       Correlation uQQ(fileIn, "uQ", {"u_rec_sgnl_RESCALED", subevents.at(i), subevents.at(j)}, components_211);
//       auto uQ_R1_sub3 = uQQ / R1_sub3.at(i) / R1_sub3.at(j) * 2.;
//       auto uQ_R1_sub4_sts_p = uQQ / R1_sub4_sts_p.at(i) / R1_sub4_sts_p.at(j) * 2.;
//       auto uQ_R1_sub4_sts_pipos = uQQ / R1_sub4_sts_pipos.at(i) / R1_sub4_sts_pipos.at(j) * 2.;
//       
//       fileOut->cd("v2/uQQ_R1_sub3");
//       uQ_R1_sub3.Save("v2.u_rec_sgnl_RESCALED." + subevents.at(i) + "." + subevents.at(j) + ".res_sub3");
//       fileOut->cd("v2/uQQ_R1_sub4_sts_p");
//       uQ_R1_sub4_sts_p.Save("v2.u_rec_sgnl_RESCALED." + subevents.at(i) + "." + subevents.at(j) + ".res_sub4_sts_p");
//       fileOut->cd("v2/uQQ_R1_sub4_sts_pipos");
//       uQ_R1_sub4_sts_pipos.Save("v2.u_rec_sgnl_RESCALED." + subevents.at(i) + "." + subevents.at(j) + ".res_sub4_sts_pipos");      
//     }
//   }  

    
  fileOut->Close();
}