TH3F* Average(TH3F* h1, TH3F* h2);

void TH3_Averager()
{
//     TFile* file_XX = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.root");
//     TFile* file_YY = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.YY.root");
    
//     TFile* file_XX = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.third1.root");
//     TFile* file_YY = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.YY.third1.root");
    
    TFile* file_XX = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcv1.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.third1.root");
    TFile* file_YY = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcv1.apr20.dcmqgsm.nopid.lightcuts1.set4.YY.third1.root");
    
    std::vector<std::pair<std::string, std::string>> histonames
    {
      {"parameters", "hsignal"},
      {"parameters", "hbckgr_0"},
      {"parameters", "hbckgr_1"},
      {"",           "hv1_mc"}
    };
    
    TFile* fileOut = TFile::Open("out.average.root", "recreate");
    
    for(auto& hn : histonames)
    {
      std::string histoName;
      if(hn.first == "")
        histoName = hn.second;
      else
        histoName = hn.first + "/" + hn.second;
      
      TH3F* h_XX = file_XX->Get<TH3F>(histoName.c_str());
      if(h_XX == nullptr) continue;
      
      TH3F* h_YY = file_YY->Get<TH3F>(histoName.c_str());
      
      TH3F* h_ave = Average(h_XX, h_YY);
      
      if(fileOut->GetDirectory(hn.first.c_str()) == nullptr)
        fileOut -> mkdir(hn.first.c_str());
      fileOut -> cd(hn.first.c_str());
      
      h_ave->Write(hn.second.c_str());
    }
  
  fileOut -> Close();  
}

TH3F* Average(TH3F* h1, TH3F* h2)
{
  TH3F* haverage = (TH3F*) h1 -> Clone();
  
  haverage -> Reset();
  
  for(int i=0; i<=haverage->GetNbinsX(); i++)
    for(int j=0; j<=haverage->GetNbinsY(); j++)
      for(int k=0; k<=haverage->GetNbinsZ(); k++)
      {
        const float average = (h1->GetBinContent(i, j, k) + h2->GetBinContent(i, j, k))/2.;
        const float error1 = h1->GetBinError(i, j, k);
        const float error2 = h2->GetBinError(i, j, k);
        const float error_average = TMath::Sqrt(error1*error1 + error2*error2)/2.;
        
        haverage -> SetBinContent(i, j, k, average);
        haverage -> SetBinError(i, j, k, error_average);
      }
    
  return haverage;
}