void TH3_Divider()
{
  TFile* file_fit = TFile::Open("/home/user/cbmdir/QnAnalysis/build/src/QnAnalysisDiscriminator/out.fitter.root");
  TFile* file_mc_fit = TFile::Open("/home/user/cbmdir/QnAnalysis/build/src/QnAnalysisDiscriminator/out.mcfitter.root");
  
  std::vector<std::string> histonames{"hsignal", "hbckgr_0", "hbckgr_1"};
  
  TFile* fileOut = TFile::Open("out.ratio.root", "recreate");
  TDirectory* dirRatio = fileOut->mkdir("ratios");
  TDirectory* dirDiff = fileOut->mkdir("diff_from_ones");
  
  for(auto& histoname : histonames)
  {
    TH3F* hfit = file_fit -> Get<TH3F>(("parameters/"+histoname).c_str());
    TH3F* hmcfit = file_mc_fit -> Get<TH3F>(("parameters/"+histoname).c_str());
    TH3F* hratio = (TH3F*) hfit -> Clone();
    TH3F* hsigma = (TH3F*) hfit -> Clone();
    
    hratio -> Sumw2();
    hratio -> Divide(hmcfit);
    hratio -> SetTitle((histoname+"_RATIO").c_str());
    hratio -> SetName((histoname+"_ratio").c_str());
    
    dirRatio -> cd();
    hratio -> Write();
    
    hsigma -> Reset();
    
    for(int i=0; i<=hratio->GetNbinsX(); i++)
      for(int j=0; j<=hratio->GetNbinsY(); j++)
        for(int k=0; k<=hratio->GetNbinsZ(); k++)
          hsigma -> SetBinContent(i, j, k, (hratio->GetBinContent(i, j, k) - 1.) / hratio->GetBinError(i, j, k));
    
    hsigma -> SetTitle((histoname+"_DIFF").c_str());
    hsigma -> SetName((histoname+"_diff_from_one").c_str());
    
    dirDiff -> cd();
    hsigma -> Write();
  }
  
  fileOut -> Close();
}