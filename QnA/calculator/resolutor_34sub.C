void resolutor_34sub()
{
  std::string fileName = "/home/user/cbmdir/working/qna/resolutions/OlegSelection_sumw/psd.dcmqgsm.12agev.root";
  
  TFile* fileIn = TFile::Open(fileName.c_str());
  
//   std::vector<std::string> steps{"PLAIN", "RECENTERED", "TWIST", "RESCALED"};
  std::vector<std::string> steps{"RECENTERED"};
  std::vector<std::string> components{"x1x1", "y1y1"};
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  fileOut->cd();
  
  // --------------------- 3 sub -----------------------------------------------------------
  for(auto& st : steps)
  {
    std::vector<std::string> subevents{"psd1_" + st, "psd2_" + st, "psd3_" + st};
  
    for(auto& se : subevents)
    {
      auto res = Functions::VectorResolutions3S( fileIn, "QQ", se, subevents, components );
      for(auto& r : res)
        r.Save(se);
    }
  }
  // ----------------------------------------------------------------------------------------
  
  //--------------------- 4 sub, Formula 2.12 E.K. --------------------------------------------------------------------------------------------------------------
  auto res_psd1_via_sts_p = Functions::VectorResolutions3S( fileIn, "QQ", "psd1_RECENTERED", {"sts_p_RESCALED", "psd3_RECENTERED"}, components );
  res_psd1_via_sts_p.at(0).Save("res_psd1_via_sts_p");
  
  auto res_psd3_via_sts_p = Functions::VectorResolutions3S( fileIn, "QQ", "psd3_RECENTERED", {"sts_p_RESCALED", "psd1_RECENTERED"}, components );
  res_psd3_via_sts_p.at(0).Save("res_psd3_via_sts_p");
  
  auto res_psd1_via_sts_pipos = Functions::VectorResolutions3S( fileIn, "QQ", "psd1_RECENTERED", {"sts_pipos_RESCALED", "psd3_RECENTERED"}, components );
  res_psd1_via_sts_pipos.at(0).Save("res_psd1_via_sts_pipos");
  
  auto res_psd3_via_sts_pipos = Functions::VectorResolutions3S( fileIn, "QQ", "psd3_RECENTERED", {"sts_pipos_RESCALED", "psd1_RECENTERED"}, components );
  res_psd3_via_sts_pipos.at(0).Save("res_psd3_via_sts_pipos");
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
 
  //--------------------- 4 sub, Formula 2.13 E.K. --------------------------------------------------------------------------------------------------------------
  auto res_sts_p = Functions::VectorResolutions3S( fileIn, "QQ", "sts_p_RESCALED", {"psd1_RECENTERED", "psd3_RECENTERED"}, components );
  Correlation psd2_sts_p( fileIn, "QQ", {"psd2_RECENTERED", "sts_p_RESCALED"}, components);
  auto res_psd2_via_sts_p = psd2_sts_p/res_sts_p.at(0)*(-2.);
  res_psd2_via_sts_p.Save("res_psd2_via_sts_p");
  
  auto res_sts_pipos = Functions::VectorResolutions3S( fileIn, "QQ", "sts_pipos_RESCALED", {"psd1_RECENTERED", "psd3_RECENTERED"}, components );
  Correlation psd2_sts_pipos( fileIn, "QQ", {"psd2_RECENTERED", "sts_pipos_RESCALED"}, components);
  auto res_psd2_via_sts_pipos = psd2_sts_pipos/res_sts_pipos.at(0)*(-2.);
  res_psd2_via_sts_pipos.Save("res_psd2_via_sts_pipos");
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  fileOut->Close();
}