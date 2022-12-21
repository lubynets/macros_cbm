void fit_dv1dy() {
  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen + "/v1andR1.stf." + evegen + ".root";

  TFile* fileIn = TFile::Open(fileName.c_str());

  std::vector<std::string> particles{"lambda", "kshort", "pipos", "pineg"};
  std::vector<std::string> components{"x1x1", "y1y1"};
  std::vector<std::string> subevents{"psd1_RECENTERED",
                                     "psd2_RECENTERED",
                                     "psd3_RECENTERED",
                                     "spec1_prim_PLAIN",
                                     "spec2_prim_PLAIN",
                                     "spec3_prim_PLAIN",
                                    };

  std::string axistofit = "SimParticles_rapidity";
  const float midrapidity = 1.6217901;

  TFile* fileOut = TFile::Open(("dv1dy.stf." + evegen + ".root").c_str(), "recreate");
  fileOut->cd();
  for(auto& pa : particles) {
    fileOut->mkdir(pa.c_str());
    fileOut->cd(pa.c_str());
    for(auto& se : subevents) {
      for (auto& co : components) {
        Qn::DataContainerStatCalculate* v1sim = (Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uPsi/v1.uPsi." + co).c_str());
        Qn::DataContainerStatCalculate* v1rec = (Qn::DataContainerStatCalculate*) fileIn->Get(("v1/" + pa + "/uQ_R1/v1.uQ_R1." + se +"." + co).c_str());

        const double fitaxis_lo = v1sim->GetAxis(axistofit.c_str()).GetFirstBinEdge();
        const double fitaxis_hi = v1sim->GetAxis(axistofit.c_str()).GetLastBinEdge();

        Qn::DataContainerStatCalculate v1sim_rebinned = v1sim->Rebin({axistofit, 1, fitaxis_lo, fitaxis_hi});
        Qn::DataContainerStatCalculate v1sim_reduced = v1sim_rebinned.Select({axistofit, 1, fitaxis_lo, fitaxis_hi});
        Qn::DataContainerStatCalculate v1rec_rebinned = v1rec->Rebin({axistofit, 1, fitaxis_lo, fitaxis_hi});
        Qn::DataContainerStatCalculate v1rec_reduced = v1rec_rebinned.Select({axistofit, 1, fitaxis_lo, fitaxis_hi});

        Qn::DataContainerStatDiscriminator v1sim_slope;
        Qn::DataContainerStatDiscriminator v1sim_intercept;
        Qn::DataContainerStatDiscriminator v1rec_slope;
        Qn::DataContainerStatDiscriminator v1rec_intercept;

        v1sim_slope.AddAxes(v1sim_reduced.GetAxes());
        v1sim_intercept.AddAxes(v1sim_reduced.GetAxes());
        v1rec_slope.AddAxes(v1sim_reduced.GetAxes());
        v1rec_intercept.AddAxes(v1sim_reduced.GetAxes());

        GraphExtractor gex_sim;
        gex_sim.SetDataContainer(v1sim);
        gex_sim.SetSelectAxis(axistofit.c_str());
        GraphExtractor gex_rec;
        gex_rec.SetDataContainer(v1rec);
        gex_rec.SetSelectAxis(axistofit.c_str());

        for (int i = 0; i < v1sim_reduced.size(); i++) {
          TGraphErrors* gr_sim = gex_sim.GetGraph(v1sim_reduced.GetIndex(i));
          TGraphErrors* gr_rec = gex_rec.GetGraph(v1rec_reduced.GetIndex(i));

          TF1* fsim = new TF1("fsim", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
          TF1* frec = new TF1("frec", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
          fsim->FixParameter(2, midrapidity);
          frec->FixParameter(2, midrapidity);

          gr_sim->Fit(fsim, "0");
          gr_rec->Fit(frec, "0");

          v1sim_intercept[i].SetVEW(fsim->GetParameter(0), fsim->GetParError(0));
          v1sim_slope[i].SetVEW(fsim->GetParameter(1), fsim->GetParError(1));
          v1rec_intercept[i].SetVEW(frec->GetParameter(0), frec->GetParError(0));
          v1rec_slope[i].SetVEW(frec->GetParameter(1), frec->GetParError(1));
        }
        v1sim_intercept.Write(("v1sim_intercept." + se + "." + co).c_str());
        v1sim_slope.Write(("v1sim_slope." + se + "." + co).c_str());
        v1rec_intercept.Write(("v1rec_intercept." + se + "." + co).c_str());
        v1rec_slope.Write(("v1rec_slope." + se + "." + co).c_str());
      }
    }
  }
  fileOut->Close();
}
