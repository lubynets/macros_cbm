//
// Created by mikhail on 8/1/21.
//
void example_v1(){
  auto file = TFile::Open( "~/Correlations/au123_p_pi_2021_09_10.root" );
  std::vector<std::string> ep_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED" };
  std::vector<std::string> res_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED", "Mf_RESCALED", "Mb_RESCALED" };
  std::vector<std::string> components{"x1x1", "y1y1"};

  auto file_out = TFile::Open( "au123_p_pi_2021_09_15.root", "recreate" );
  file_out->cd();
  file_out->mkdir("resolutions");
  file_out->mkdir("protons");
  file_out->mkdir("pi_pos");
  file_out->mkdir("pi_neg");

  for( auto qa : ep_vectors ){
    Correlation protons_qa( file, "SP/u1Q1", {"protons_RESCALED", qa}, {"x1x1", "y1y1"});
    Correlation pi_pos_qa( file, "SP/u1Q1", {"pi_pos_RESCALED", qa}, {"x1x1", "y1y1"});
    Correlation pi_neg_qa( file, "SP/u1Q1", {"pi_neg_RESCALED", qa}, {"x1x1", "y1y1"});
    auto res_v = Functions::VectorResolutions3S( file, "SP/Q1Q1", qa, res_vectors, components );
    for( auto R1 : res_v ) {
      auto v1_protons = protons_qa * 2 / R1;
      auto v1_pi_pos = pi_pos_qa * 2 / R1;
      auto v1_pi_neg = pi_neg_qa * 2 / R1;
      file_out->cd("resolutions");
      R1.Save("R1."+R1.Title());
      file_out->cd("protons");
      v1_protons.Save("v1."+R1.Title());
      file_out->cd("pi_pos");
      v1_pi_pos.Save("v1."+R1.Title());
      file_out->cd("pi_neg");
      v1_pi_neg.Save("v1."+R1.Title());
    }
  }
  file_out->Close();
  file->Close();
}