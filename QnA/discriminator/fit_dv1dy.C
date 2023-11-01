#include "fit_dv1dy.h"

void fit_dv1dy() {
  //**** aXmass
  excludedFolders = {"Chi2s", "Entries", "R1"};
  excludedObjects = {"bckgr"};
  axestofit = {"ReconstructedParticles_rapidity", "SimParticles_rapidity", "y"};
  axestoignore = {"ReconstructedParticles_mass"};
  axestoslice = {"ReconstructedParticles_pT", "SimParticles_pT", "pT"};


//   //**** stf
//   excludedFolders = {"R1", "QQ", "R1R1", "pipos", "pineg", "uQ"};
//   excludedObjects = {"x1y1", "y1x1", "_all_PLAIN", "_odd", "_even"};
//   axestofit = {"ReconstructedParticles_rapidity", "SimParticles_rapidity", "y"};
//   axestoslice = {"SimEventHeader_centrality_impactpar", "SimParticles_pT"};
//   rangestoslice.resize(axestoslice.size());
//   rangestoslice.at(0) = {0, 6, 10, 15, 20, 30, 40, 70};

//   std::string evegen = "dcmqgsm"; std::string pbeam = "12"; midrapidity = 1.6217901;
  std::string evegen = "dcmqgsm"; std::string pbeam = "3.3"; midrapidity = 0.985344;
//   std::string evegen = "urqmd"; std::string pbeam = "12"; midrapidity = 1.6217901;

//   std::string pdg = "3122"; std::string cuts = "lc1";
  std::string pdg = "310"; std::string cuts = "lc1";
//   std::string pdg = "3312"; std::string cuts = "dc";

  if(pdg!="3312" && pbeam=="3.3") cuts = "oc1";

  //   std::string is_fine_pt = "";
  std::string is_fine_pt = "_finept";

  //**** aXmass
  fileInPath = "/home/oleksii/cbmdir/working/qna/aXmass";
//   fileInName = "vR." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + is_fine_pt + ".root";
  fileInName = "of." + evegen + "." + pbeam + "agev." + cuts + "." + pdg + is_fine_pt + ".root";

//   //**** stf
//   fileInPath = "/home/oleksii/cbmdir/working/qna/simtracksflow/" + evegen  + "/" + pbeam + "agev";
//   fileInName = "v1andR1.stf." + evegen + "." + pbeam + "agev.root";

  TFile* fileIn = TFile::Open((fileInPath + "/" + fileInName).c_str(), "read");

  bool is_calculate_average = true;

  TDirectory* dir = fileIn->GetDirectory("");

  BuildObjectList(dir);

  fileOutName = fileInName;
  fileOutName.ReplaceAll("vR", "vR.dv1dy");
  fileOutName.ReplaceAll("v1andR1", "v1andR1.dv1dy");
  fileOutName.ReplaceAll("of", "of.dv1dy");
  TFile* fileOut = TFile::Open(fileOutName, "recreate");

//   std::cout << "object_paths.size() = " << object_paths.size() << "\n"; throw;
  for(int iObj=0; iObj<object_paths.size(); iObj++){
    Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>* obj{nullptr};
    if(stat_types.at(iObj) == 1) obj = new Qn::DataContainerStatDiscriminator(*(fileIn->Get<Qn::DataContainerStatCollect>(object_paths.at(iObj) + "/" + object_names.at(iObj))));
    if(stat_types.at(iObj) == 2) obj = new Qn::DataContainerStatDiscriminator(*(fileIn->Get<Qn::DataContainerStatCalculate>(object_paths.at(iObj) + "/" + object_names.at(iObj))));
    if(stat_types.at(iObj) == 3) obj = fileIn->Get<Qn::DataContainerStatDiscriminator>(object_paths.at(iObj) + "/" + object_names.at(iObj));

//     //**** stf
//     if(evegen == "dcmqgsm" && pbeam == "12") {
//       if(object_paths.at(iObj).Contains("lambda")) rangestoslice.at(1) = {0, 0.4, 1.0, 1.6};
//       if(object_paths.at(iObj).Contains("kshort")) rangestoslice.at(1) = {0, 0.4, 0.8, 1.6};
//       if(object_paths.at(iObj).Contains("xi")) rangestoslice.at(1) = {0, 0.4, 1.0, 1.6};
//     } else if(evegen == "urqmd" && pbeam == "12") {
//       if(object_paths.at(iObj).Contains("lambda")) rangestoslice.at(1) = {0, 0.8, 1.2, 1.6};
//       if(object_paths.at(iObj).Contains("kshort")) rangestoslice.at(1) = {0, 0.8, 1.2, 1.6};
//       if(object_paths.at(iObj).Contains("xi")) rangestoslice.at(1) = {0, 0.8, 1.2, 1.6};
//     } else if(evegen == "dcmqgsm" && pbeam == "3.3") {
//       if(object_paths.at(iObj).Contains("lambda")) rangestoslice.at(1) = {0, 0.4, 0.8, 1.2};
//       if(object_paths.at(iObj).Contains("kshort")) rangestoslice.at(1) = {0, 0.8, 1.2, 1.6};
//       if(object_paths.at(iObj).Contains("xi")) rangestoslice.at(1) = {0, 0.8, 1.2, 1.6};
//     }


    if(!CheckAxes(obj)) {
      delete obj;
      continue;
    }
    std::cout << object_paths.at(iObj) << " " << object_names.at(iObj) << "\n";

    auto fitted_obj = Fit(*obj);
    CD(fileOut, (object_paths.at(iObj) + "/slope").Data());
    fitted_obj.first.Write(object_names.at(iObj));
    CD(fileOut, (object_paths.at(iObj) + "/intercept").Data());
    fitted_obj.second.Write(object_names.at(iObj));

    if(is_calculate_average) {
      if(object_names.at(iObj).Contains("y1y1")) {
        delete obj;
        continue;
      }
      Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>* obj_2{nullptr};
      auto object_name_2 = object_names.at(iObj);
      object_name_2.ReplaceAll("x1x1", "y1y1");

      if(stat_types.at(iObj) == 1) obj_2 = new Qn::DataContainerStatDiscriminator(*(fileIn->Get<Qn::DataContainerStatCollect>(object_paths.at(iObj) + "/" + object_name_2)));
      if(stat_types.at(iObj) == 2) obj_2 = new Qn::DataContainerStatDiscriminator(*(fileIn->Get<Qn::DataContainerStatCalculate>(object_paths.at(iObj) + "/" + object_name_2)));
      if(stat_types.at(iObj) == 3) obj_2 = fileIn->Get<Qn::DataContainerStatDiscriminator>(object_paths.at(iObj) + "/" + object_name_2);

      object_name_2.ReplaceAll("y1y1", "ave");

      auto fitted_obj_2 = Fit(*obj, *obj_2);
      CD(fileOut, (object_paths.at(iObj) + "/slope").Data());
      fitted_obj_2.first.Write(object_name_2);
      CD(fileOut, (object_paths.at(iObj) + "/intercept").Data());
      fitted_obj_2.second.Write(object_name_2);

      delete obj_2;
    }
    delete obj;
  }
  fileOut->Close();
}

std::pair<Qn::DataContainerStatDiscriminator, Qn::DataContainerStatDiscriminator> Fit(Qn::DataContainerStatDiscriminator dcIn) {
  const int Nsamples = dcIn.At(0).GetSampleMeans().size();

//   dcIn = dcIn.Rebin({axestoslice.at(0).c_str(), rangestoslice.at(0)});
//   dcIn = dcIn.Rebin({axestoslice.at(1).c_str(), rangestoslice.at(1)});

  const double fitaxis_lo = dcIn.GetAxis(axistofit.c_str()).GetFirstBinEdge();
  const double fitaxis_hi = dcIn.GetAxis(axistofit.c_str()).GetLastBinEdge();

  Qn::DataContainerStatDiscriminator dc_reduced = dcIn.Rebin({axistofit, 1, fitaxis_lo, fitaxis_hi});
  dc_reduced = dc_reduced.Select({axistofit, 1, fitaxis_lo, fitaxis_hi});

  Qn::DataContainerStatDiscriminator slope;
  Qn::DataContainerStatDiscriminator intercept;

  slope.AddAxes(dc_reduced.GetAxes());
  intercept.AddAxes(dc_reduced.GetAxes());

  GraphExtractor gex;
  gex.SetDataContainer(&dcIn);
  gex.SetSelectAxis(axistofit.c_str());

  for (int i = 0; i < dc_reduced.size(); i++) {
    gex.ReduceDataContainerToBin(dc_reduced.GetIndex(i));
    TGraphErrors* gr = gex.GetGraph();

    TF1* ffit = new TF1("ffit", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
    ffit->FixParameter(2, midrapidity);

    gr->Fit(ffit, "0");

    intercept[i].SetVEW(ffit->GetParameter(0), ffit->GetParError(0));
    slope[i].SetVEW(ffit->GetParameter(1), ffit->GetParError(1));
    delete ffit;

    std::vector<TGraphErrors*> grs = gex.GetSamplesGraphs();
    std::vector<double> samples_weights = gex.GetSamplesWeights();
    for(int isample = 0; isample<Nsamples; isample++) {
      ffit = new TF1("ffit", "[0]+[1]*(x-[2])", fitaxis_lo, fitaxis_hi);
      ffit->FixParameter(2, midrapidity);

      grs.at(isample)->Fit(ffit, "0");

      intercept[i].AddSampleMean(ffit->GetParameter(0));
      slope[i].AddSampleMean(ffit->GetParameter(1));
      intercept[i].AddSampleWeight(samples_weights.at(isample));
      slope[i].AddSampleWeight(samples_weights.at(isample));
      delete ffit;
      delete grs.at(isample);
    }
  }

  return std::make_pair(slope, intercept);
}
