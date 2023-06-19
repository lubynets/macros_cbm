void flow_simple() {

  std::string evegen = "dcmqgsm"; std::string pbeam = "3.3";
  const float mu = 0.497611; std::string particle = "kshort"; std::string pdg = "310";

  std::string v1filename = "/home/oleksii/cbmdir/working/qna/aXmass/vR." + evegen + "." + pbeam + "agev.oc1." + pdg + ".root";
  std::string shapefilename = "/home/oleksii/cbmdir/working/qna/aXmass/shapes/shapesimple." + evegen + "." + pbeam + "agev.oc1." + pdg + ".root";

  TFile* shapefile = TFile::Open(shapefilename.c_str(), "read");
  Qn::DataContainerStatDiscriminator* shcntr = (Qn::DataContainerStatDiscriminator*) shapefile->Get( "Purities" );

  TFile* v1file = TFile::Open(v1filename.c_str(), "read");

  std::string invmassaxis = "ReconstructedParticles_mass";

  std::vector<std::string> components = {"x1x1", "y1y1"};

  TFile* fileOut = TFile::Open("of.simple.root", "recreate");

  std::vector<std::pair<std::string, std::string>> inputs {
    {"uPsi", ""},
    {"uQ_R1_MC", "_res_MC"},
    {"uQ_R1_sub3", "_res_sub3"},
    {"uQ_R1_sub4", "_res_sub4_sts_pipos"}
  };

  for(auto& ip : inputs) {
      std::vector<std::string> subevents;
      if(ip.second == "") subevents = {"Q_psi"};
      else                subevents = {"psd1", "psd2", "psd3"};

      for(auto& se : subevents) {
        for (auto& co : components) {

          std::string corr_name = "v1/" + ip.first + "/v1.u_rec." + se + ip.second + "." + co;
          std::cout << corr_name << "\n";
          Qn::DataContainerStatCalculate lambda_psi_pre =
          Qn::DataContainerStatCalculate(*(Qn::DataContainerStatCalculate*) v1file->Get(corr_name.c_str()));
          auto lambda_psi = lambda_psi_pre;

          FitterSimple fitter;
          fitter.SetDataContainer(&lambda_psi);
          fitter.SetPurityContainer(shcntr);
          fitter.SetSelectAxis(invmassaxis.c_str());
          fitter.AutoSetEdgesOfLMR();
          fitter.Calculate();

          Qn::DataContainerStatDiscriminator* dc_signal = fitter.GetResult();
          CD(fileOut, ("Pars/" + ip.first + "/" + se + ip.second).c_str());
          dc_signal->Write(("signal." + co).c_str());
          delete dc_signal;
        } // components
      } // subevents
  } // inputs
  fileOut->Close();
}
