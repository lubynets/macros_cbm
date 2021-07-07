TH3F* Subtract(TH3F* h1, TH3F* h2);

void TH3_Subtractor()
{
  TFile* file_fit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.AVE.third1.root");
  TFile* file_mcfit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.AVE.third1.root");
  TFile* file_mcv1 = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcv1.apr20.dcmqgsm.nopid.lightcuts1.set4.AVE.third1.root");
  
  struct infotype
  {
    std::string name_;
    std::string folder_;
    bool is_mcv1_;
  };
  
  std::vector<infotype> infotypes
  {
    {"hsignal",  "sgnl",            true},
    {"hbckgr_0", "bckgr/intercept", false},
    {"hbckgr_1", "bckgr/slope",     false}
  };
  
  TH3F* h_fit;
  TH3F* h_mcfit;
  TH3F* h_mcv1;
  
  TH3F* diff_fit_mcfit;
  TH3F* diff_fit_mcv1;
  TH3F* diff_mcfit_mcv1;
  
  TFile* fileOut = TFile::Open("out.difference.root", "recreate");
  
  for(auto it : infotypes)
  {
                    h_fit = file_fit -> Get<TH3F>(("parameters/" + it.name_).c_str());
                    h_mcfit = file_mcfit -> Get<TH3F>(("parameters/" + it.name_).c_str());
    if(it.is_mcv1_) h_mcv1 = file_mcv1 -> Get<TH3F>("hv1_mc"); 
    
    std::string dirname = it.folder_;
    fileOut -> mkdir(dirname.c_str());
    fileOut -> cd(dirname.c_str()); 
    
                      diff_fit_mcfit = Subtract(h_fit, h_mcfit);
    if(it.is_mcv1_) { diff_fit_mcv1 = Subtract(h_fit, h_mcv1);
                      diff_mcfit_mcv1 = Subtract(h_mcfit, h_mcv1); }
    
                      diff_fit_mcfit -> Write("diff_fit_mcfit");
    if(it.is_mcv1_) { diff_fit_mcv1 -> Write("diff_fit_mcv1");
                      diff_mcfit_mcv1 -> Write("diff_mcfit_mcv1"); }    
  }
    
  fileOut -> Close();
}

TH3F* Subtract(TH3F* h1, TH3F* h2)
{
  TH3F* hdifference = (TH3F*) h1 -> Clone();
  
  hdifference -> Reset();
  
  for(int i=0; i<=hdifference->GetNbinsX(); i++)
    for(int j=0; j<=hdifference->GetNbinsY(); j++)
      for(int k=0; k<=hdifference->GetNbinsZ(); k++)
      {
        const float difference = h1->GetBinContent(i, j, k) - h2->GetBinContent(i, j, k);
        const float error1 = h1->GetBinError(i, j, k);
        const float error2 = h2->GetBinError(i, j, k);
        const float error_difference = TMath::Sqrt(error1*error1 + error2*error2);
        
        hdifference -> SetBinContent(i, j, k, difference/error_difference);
      }
  
//   hdifference -> SetTitle((histoname+"_diff").c_str());
//   hdifference -> SetName((histoname+"_diff_interms_err").c_str());
  
  return hdifference;
}