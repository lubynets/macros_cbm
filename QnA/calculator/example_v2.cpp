//
// Created by mikhail on 8/1/21.
//
void example_v2(){
  auto file = TFile::Open( "~/Correlations/au123_all_elliptic_2021_07_28.root" );
  std::vector<std::string> ep_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED" };
  std::vector<std::string> res_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED", "Mf_RESCALED", "Mb_RESCALED" };
  std::vector<std::string> components{"x1x1", "y1y1"};

  auto file_out = TFile::Open( "au123_all_elliptic_2021_07_28.root", "recreate" );
  file_out->cd();

  for( size_t i=0; i<ep_vectors.size(); ++i ){
    auto qa = ep_vectors.at(i);
    auto res_va = Functions::VectorResolutions3S( file, "SP/QQ", qa, res_vectors, components );
    for ( size_t j=i+1; j<ep_vectors.size(); ++j ){
      auto qb = ep_vectors.at(j);
      auto res_vb = Functions::VectorResolutions3S( file, "SP/QQ", qb, res_vectors, components );
      Correlation u_qa_qb( file, "SP/u2Q1Q1", {"protons_RESCALED", qa, qb}, {"x2x1x1", "y2x1y1", "y2y1x1" ,"x2y1y1"});
      for( auto Ra : res_va )
        for( auto Rb : res_vb ){
          auto RaxRb = MatrixMultiply( Ra, Rb );
          auto v2 = u_qa_qb/RaxRb*std::vector{4.0, 4.0, 4.0, -4.0};
          RaxRb.Save( "R2."+Ra.Title()+"."+Rb.Title() );
          v2.Save( "v2."+Ra.Title()+"."+Rb.Title() );
        }

    }
  }
  file_out->Close();
  file->Close();
}