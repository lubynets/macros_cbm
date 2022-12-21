//
// Created by mikhail on 8/1/21.
//
void v2_calculator(){
  auto file = TFile::Open( "~/Correlations/au123_p_pi_2021_09_10.root" );
//  std::vector<std::string> ep_vectors{ "W1_RESCALED", "W3_RESCALED" };
  std::vector<std::string> res_vectors{ "Mf_RESCALED", "Mb_RESCALED" };
  std::vector<std::string> components{"x1x1", "y1y1"};

  auto file_out = TFile::Open( "au123_pi_pos_2021_09_10.root", "recreate" );
  file_out->cd();

  Correlation u_w1_w3( file, "SP/u2Q1Q1", {"pi_pos_RESCALED", "W1_RESCALED", "W3_RESCALED"}, {"x2x1x1", "x2y1y1", "y2x1y1", "y2y1x1"});
  u_w1_w3[2] = u_w1_w3[0]*u_w1_w3[1]*-1;
  auto vector_names = u_w1_w3.GetComponentNames();
  u_w1_w3.GetComponentNames()[2] = std::string("y2x1y1.y2y1x1");
  u_w1_w3.RemoveComponent(3);
  Correlation w1_w3( file, "SP/Q1Q1", {"W1_RESCALED", "W3_RESCALED"}, {"x1x1", "y1y1"});
  w1_w3.AddComponent(w1_w3[0]*w1_w3[1], "x1x1.y1y1");
  auto v2 = u_w1_w3 / w1_w3 * std::vector{2.0, -2.0, 4.0};
  v2[2] = Sqrt(v2[2]);
  file_out->cd();
  v2.Save( "v2.pi_pos_RESCALED.W1_RESCALED.W3_RESCALED" );

  Correlation u_w1_w2( file, "SP/u2Q1Q1", {"pi_pos_RESCALED", "W1_RESCALED", "W3_RESCALED"}, {"x2x1x1", "x2y1y1"});
  for( auto mdc : res_vectors ){
    Correlation w2_mdc( file, "SP/Q1Q1", {"W2_RESCALED", mdc}, {"x1x1", "y1y1"});
    Correlation w3_mdc( file, "SP/Q1Q1", {"W3_RESCALED", mdc}, {"x1x1", "y1y1"});
    w1_w3 = Correlation( file, "SP/Q1Q1", {"W1_RESCALED", "W3_RESCALED"}, {"x1x1", "y1y1"});
    v2 = u_w1_w2 / ( w2_mdc * w1_w3 / w3_mdc ) * std::vector{2.0, -2.0} ;
    v2.Save( "v2.pi_pos_RESCALED.W1_RESCALED.W2_RESCALED("+mdc+")" );
  }
  Correlation u_w2_w3( file, "SP/u2Q1Q1", {"pi_pos_RESCALED", "W2_RESCALED", "W3_RESCALED"}, {"x2x1x1", "x2y1y1"});
  for( auto mdc : res_vectors ){
    Correlation w2_mdc( file, "SP/Q1Q1", {"W2_RESCALED", mdc}, {"x1x1", "y1y1"});
    Correlation w1_mdc( file, "SP/Q1Q1", {"W1_RESCALED", mdc}, {"x1x1", "y1y1"});
    w1_w3 = Correlation( file, "SP/Q1Q1", {"W1_RESCALED", "W3_RESCALED"}, {"x1x1", "y1y1"});
    v2 = u_w2_w3 / ( w2_mdc * w1_w3 / w1_mdc ) * std::vector{2.0, -2.0} ;
    v2.Save( "v2.pi_pos_RESCALED.W2_RESCALED.W3_RESCALED("+mdc+")" );
  }

  std::vector<std::string> ep_vectors{ "W1_RESCALED", "W2_RESCALED", "W3_RESCALED" };
  for( auto wall : ep_vectors ) {
    for (auto mdc : res_vectors) {
      Correlation u_wall_mdc( file, "SP/u2Q1Q1", {"pi_pos_RESCALED", wall, mdc}, {"x2x1x1", "x2y1y1"});
      Correlation wall_mdc(file, "SP/Q1Q1", {wall, mdc}, {"x1x1", "y1y1"});
      v2 = u_wall_mdc / wall_mdc * std::vector{2.0, -2.0};
      v2.Save("v2.pi_pos_RESCALED."+wall+"."+mdc);
    }
  }
  file_out->Close();
  file->Close();
}