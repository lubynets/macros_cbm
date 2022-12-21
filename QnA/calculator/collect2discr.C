void collect2discr() {
  
  std::string fileName = "/home/user/cbmdir/working/qna/correlations/cl.dcmqgsm.12agev.defcuts.3122.set1.sgnl_1.root";
  TFile* fileIn = TFile::Open(fileName.c_str());
  
  std::vector<std::string> components{"x1x1", "y1y1", "x2x2", "y2y2"};
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  fileOut->mkdir("rec/RESCALED");
  fileOut->mkdir("sim");
  
  for(auto& co : components) {
    // rec/RESCALDED
    std::string containerName = "rec/RESCALED/u_rec_RESCALED.Q_psi_PLAIN." + co;
    Qn::DataContainerStatCollect* scol{nullptr};
    fileIn -> GetObject(containerName.c_str(), scol);
//     Qn::DataContainerStatCalculate* scalc = new Qn::DataContainerStatCalculate(*scol);
    Qn::DataContainerStatDiscriminator* sdisc = new Qn::DataContainerStatDiscriminator(*scol);
    fileOut->cd("rec/RESCALED");
    sdisc -> Write(("u_rec_RESCALED.Q_psi_PLAIN." + co).c_str());
    
    // sim
    containerName = "sim/u_sim_PLAIN.Q_psi_PLAIN." + co;
    fileIn -> GetObject(containerName.c_str(), scol);
//     scalc = new Qn::DataContainerStatCalculate(*scol);
    sdisc = new Qn::DataContainerStatDiscriminator(*scol);
    fileOut->cd("sim");
    sdisc -> Write(("u_sim_PLAIN.Q_psi_PLAIN." + co).c_str());
  }

  fileOut -> Close();
}