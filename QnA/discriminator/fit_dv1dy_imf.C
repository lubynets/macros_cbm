void fit_dv1dy_imf() {
  std::string evegen = "dcmqgsm";
//   std::string evegen = "urqmd";

 std::string particle = "lambda";
//  std::string particle = "kshort";

//   bool average_comp{true}; std::vector<std::string> components{"ave"};
  bool average_comp{false}; std::vector<std::string> components{"x1x1", "y1y1"};

  std::string fileNameMc = "/home/oleksii/cbmdir/working/qna/inv_mass_flow/" + evegen + "/cl.imf." + evegen + ".12agev.root";
  std::string fileNameImf = "/home/oleksii/cbmdir/working/qna/inv_mass_flow/" + evegen + "/imfits/out.fitter." + evegen + ".12agev." + particle + ".root";

  TFile* fileInMc = TFile::Open(fileNameMc.c_str());
  if(!fileInMc) throw std::runtime_error("fileInMc absent!");
  TFile* fileInImf = TFile::Open(fileNameImf.c_str());
  if(!fileInImf) throw std::runtime_error("fileInImf absent!");

  std::string axistofit_mcmatch = "ReconstructedParticles_rapidity";
  std::string axistofit_imf = "y";
  const float midrapidity = 1.6217901;

  std::string fileOutName= "dv1dy.imf." + evegen + "." + particle + ".root";
  if(average_comp) {
    fileOutName.erase(fileOutName.length()-4, 4);
    fileOutName += "ave.root";
  }

  TFile* fileOut = TFile::Open(fileOutName.c_str(), "recreate");

  for(auto& co : components) {

    Qn::DataContainerStatCalculate v1mcmatch_in;
    if(average_comp) {
      v1mcmatch_in = Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileInMc->Get(("rec/u_" + particle + "_rec_mc_PLAIN.Q_psi_PLAIN.x1x1").c_str()))) +
                 Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileInMc->Get(("rec/u_" + particle + "_rec_mc_PLAIN.Q_psi_PLAIN.y1y1").c_str())));
      v1mcmatch_in = v1mcmatch_in/2.;
    } else {
      v1mcmatch_in = Qn::DataContainerStatCalculate(*((Qn::DataContainerStatCollect*) fileInMc->Get(("rec/u_" + particle + "_rec_mc_PLAIN.Q_psi_PLAIN." + co).c_str())));
    }

    const int Nsamples = v1mcmatch_in.At(0).GetSampleMeans().size();

    Qn::DataContainerStatCalculate v1mcmatch = v1mcmatch_in;

    double fitaxis_lo = v1mcmatch.GetAxis(axistofit_mcmatch.c_str()).GetFirstBinEdge();
    double fitaxis_hi = v1mcmatch.GetAxis(axistofit_mcmatch.c_str()).GetLastBinEdge();

    Qn::DataContainerStatCalculate v1mcmatch_rebinned = v1mcmatch.Rebin({axistofit_mcmatch, 1, fitaxis_lo, fitaxis_hi});
    Qn::DataContainerStatCalculate v1mcmatch_reduced = v1mcmatch_rebinned.Select({axistofit_mcmatch, 1, fitaxis_lo, fitaxis_hi});

    Qn::DataContainerStatDiscriminator v1mcmatch_slope;
    Qn::DataContainerStatDiscriminator v1mcmatch_intercept;

    v1mcmatch_slope.AddAxes(v1mcmatch_reduced.GetAxes());
    v1mcmatch_intercept.AddAxes(v1mcmatch_reduced.GetAxes());

    GraphExtractor gex_sim;
    gex_sim.SetDataContainer(&v1mcmatch);
    gex_sim.SetSelectAxis(axistofit_mcmatch.c_str());

    for (int i = 0; i < v1mcmatch_reduced.size(); i++) {
      gex_sim.ReduceDataContainerToBin(v1mcmatch_reduced.GetIndex(i));
      TGraphErrors* gr_sim = gex_sim.GetGraph();

      TF1* fmcmatch = new TF1("fmcmatch", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
      fmcmatch->FixParameter(2, midrapidity);

      gr_sim->Fit(fmcmatch, "0");

      v1mcmatch_intercept[i].SetVEW(fmcmatch->GetParameter(0), fmcmatch->GetParError(0));
      v1mcmatch_slope[i].SetVEW(fmcmatch->GetParameter(1), fmcmatch->GetParError(1));
      delete fmcmatch;

      std::vector<TGraphErrors*> gr_sims = gex_sim.GetSamplesGraphs();
      std::vector<double> samples_weights = gex_sim.GetSamplesWeights();
      for(int isample = 0; isample<Nsamples; isample++) {
        fmcmatch = new TF1("fmcmatch", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
        fmcmatch->FixParameter(2, midrapidity);

        gr_sims.at(isample)->Fit(fmcmatch, "0");

        v1mcmatch_intercept[i].AddSampleMean(fmcmatch->GetParameter(0));
        v1mcmatch_slope[i].AddSampleMean(fmcmatch->GetParameter(1));
        v1mcmatch_intercept[i].AddSampleWeight(samples_weights.at(isample));
        v1mcmatch_slope[i].AddSampleWeight(samples_weights.at(isample));
        delete fmcmatch;
        delete gr_sims.at(isample);
      }
    }
    v1mcmatch_intercept.Write(("v1mcmatch_intercept." + co).c_str());
    v1mcmatch_slope.Write(("v1mcmatch_slope." + co).c_str());


    Qn::DataContainerStatDiscriminator v1imf_in;
    if(average_comp) {
      v1imf_in = Qn::DataContainerStatDiscriminator(*((Qn::DataContainerStatDiscriminator*) fileInImf->Get("Params/signal.x1x1"))) +
                 Qn::DataContainerStatDiscriminator(*((Qn::DataContainerStatDiscriminator*) fileInImf->Get("Params/signal.y1y1")));
      v1imf_in = v1imf_in/2.;
    } else {
      v1imf_in = Qn::DataContainerStatDiscriminator(*((Qn::DataContainerStatDiscriminator*) fileInImf->Get(("Params/signal." + co).c_str())));
    }

    Qn::DataContainerStatDiscriminator v1imf = v1imf_in;

    fitaxis_lo = v1imf.GetAxis(axistofit_imf.c_str()).GetFirstBinEdge();
    fitaxis_hi = v1imf.GetAxis(axistofit_imf.c_str()).GetLastBinEdge();

    Qn::DataContainerStatDiscriminator v1imf_rebinned = v1imf.Rebin({axistofit_imf, 1, fitaxis_lo, fitaxis_hi});
    Qn::DataContainerStatDiscriminator v1imf_reduced = v1imf_rebinned.Select({axistofit_imf, 1, fitaxis_lo, fitaxis_hi});

    Qn::DataContainerStatDiscriminator v1imf_slope;
    Qn::DataContainerStatDiscriminator v1imf_intercept;

    v1imf_slope.AddAxes(v1imf_reduced.GetAxes());
    v1imf_intercept.AddAxes(v1imf_reduced.GetAxes());

    GraphExtractor gex_rec;
    gex_rec.SetDataContainer(&v1imf);
    gex_rec.SetSelectAxis(axistofit_imf.c_str());

    for (int i = 0; i < v1mcmatch_reduced.size(); i++) {
      gex_rec.ReduceDataContainerToBin(v1imf_reduced.GetIndex(i));
      TGraphErrors* gr_rec = gex_rec.GetGraph();

      TF1* fimf = new TF1("fimf", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
      fimf->FixParameter(2, midrapidity);

      gr_rec->Fit(fimf, "0", "");

      v1imf_intercept[i].SetVEW(fimf->GetParameter(0), fimf->GetParError(0));
      v1imf_slope[i].SetVEW(fimf->GetParameter(1), fimf->GetParError(1));
      delete fimf;

      std::vector<TGraphErrors*> gr_recs = gex_rec.GetSamplesGraphs();
      std::vector<double> samples_weights = gex_rec.GetSamplesWeights();
      for(int isample = 0; isample<Nsamples; isample++) {
        fimf = new TF1("fimf", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
        fimf->FixParameter(2, midrapidity);

        gr_recs.at(isample)->Fit(fimf, "0");

        v1imf_intercept[i].AddSampleMean(fimf->GetParameter(0));
        v1imf_slope[i].AddSampleMean(fimf->GetParameter(1));
        v1imf_intercept[i].AddSampleWeight(samples_weights.at(isample));
        v1imf_slope[i].AddSampleWeight(samples_weights.at(isample));
        delete fimf;
        delete gr_recs.at(isample);

      }
    }
    v1imf_intercept.Write(("v1imf_intercept." + co).c_str());
    v1imf_slope.Write(("v1imf_slope." + co).c_str());
  }
  fileOut->Close();
}
