TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis, int bin);
void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname);
std::string StringBinNumber(int number);

void TH3ToTH2Differs()
{
  TFile* fileIn = TFile::Open("/home/user/cbmdir/working/qna/fits/dirivatives/out.difference.set4.half2_half1.root");
  
  struct axis
  {
    std::string name_;
    std::string letter_;
    std::string title_;
    std::string id_;
    int another_first_;
    int another_second_;
    int nbins_;
  };
  
  std::vector<axis> axes
  {
    {"centrality", "C",  "Centrality, %", "x", 1, 2, -1},
    {"rapidity",   "y",  "y_{LAB}",       "y", 0, 2, -1},
    {"pT",         "pT", "p_{T}, GeV/c",  "z", 0, 1, -1}
  };
  
  struct infotype
  {
    std::string folder_;
    bool is_mcv1_;
  };
  
  std::vector<infotype> infotypes
  {
    {"sgnl",            true},
    {"bckgr/intercept", false},
    {"bckgr/slope",     false}
  }; 

  TH3F* h_fit_mcfit;
  TH3F* h_fit_mcv1;
  TH3F* h_mcfit_mcv1;
  
  TH2F* h2_fit_mcfit;
  TH2F* h2_fit_mcv1;
  TH2F* h2_mcfit_mcv1;  
  
  TFile* fileOut = TFile::Open("out.th3toth2differs.root", "recreate");
  
 for(auto it : infotypes)
  {
                      h_fit_mcfit = fileIn -> Get<TH3F>((it.folder_ + "/diff_fit_mcfit").c_str());
    if(it.is_mcv1_) { h_fit_mcv1 = fileIn -> Get<TH3F>((it.folder_ + "/diff_fit_mcv1").c_str());
                      h_mcfit_mcv1 = fileIn -> Get<TH3F>((it.folder_ + "/diff_mcfit_mcv1").c_str()); }   
                      
    axes.at(0).nbins_ = h_fit_mcfit->GetXaxis()->GetNbins();
    axes.at(1).nbins_ = h_fit_mcfit->GetYaxis()->GetNbins();
    axes.at(2).nbins_ = h_fit_mcfit->GetZaxis()->GetNbins();
    
    for(auto ax : axes)
    {
      for(int ibin=1; ibin<=ax.nbins_; ibin++)
      {
        std::string layername = StringBinNumber(ibin);
        std::string dirname = it.folder_ +"/" + ax.name_ + "/" + layername;
        fileOut -> mkdir(dirname.c_str());
        fileOut -> cd(dirname.c_str()); 
        
                          h2_fit_mcfit  = Get2DLayerFromTH3(h_fit_mcfit,  ax.id_, ibin);
        if(it.is_mcv1_) { h2_fit_mcv1   = Get2DLayerFromTH3(h_fit_mcv1,   ax.id_, ibin);
                          h2_mcfit_mcv1 = Get2DLayerFromTH3(h_mcfit_mcv1, ax.id_, ibin); }
        
                          SetAxesNames(h2_fit_mcfit , axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_);
        if(it.is_mcv1_) { SetAxesNames(h2_fit_mcv1  , axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_);
                          SetAxesNames(h2_mcfit_mcv1, axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_); }
        
                          h2_fit_mcfit  -> Write("fit_mcfit");
        if(it.is_mcv1_) { h2_fit_mcv1   -> Write("fit_mcv1");
                          h2_mcfit_mcv1 -> Write("mcfit_mcv1"); }
      }
    }
  } 
  fileOut -> Close();
}

TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis, int bin)
{
  if(axis=="X" || axis =="x" || axis=="YZ" || axis=="yz")
  {
    const int nbins_1 = histo3D->GetYaxis()->GetNbins();
    const int nbins_2 = histo3D->GetZaxis()->GetNbins();
    
    double* bin_edges_1 = new double[nbins_1+1];
    double* bin_edges_2 = new double[nbins_2+1];
    
    for(int i=0; i<=nbins_1; i++)
      bin_edges_1[i] = histo3D->GetYaxis()->GetBinLowEdge(i+1);
    for(int i=0; i<=nbins_2; i++)
      bin_edges_2[i] = histo3D->GetZaxis()->GetBinLowEdge(i+1);
    
    TH2F* histo2D = new TH2F("", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
    for(int i=1; i<=nbins_1; i++)
      for(int j=1; j<=nbins_2; j++)
      {
        histo2D -> SetBinContent(i, j, histo3D->GetBinContent(bin, i, j));
        histo2D -> SetBinError(i, j, histo3D->GetBinError(bin, i, j));
      }
      
    histo2D -> SetTitle((std::string(histo3D->GetXaxis()->GetTitle()) + " from " + std::to_string(histo3D->GetXaxis()->GetBinLowEdge(bin)) + " to " + std::to_string(histo3D->GetXaxis()->GetBinLowEdge(bin+1))).c_str());
    histo2D -> GetXaxis() -> SetTitle(histo3D->GetYaxis()->GetTitle());
    histo2D -> GetYaxis() -> SetTitle(histo3D->GetZaxis()->GetTitle());
    
    return histo2D;
  }
  
  if(axis=="Y" || axis =="y" || axis=="XZ" || axis=="xz")
  {
    const int nbins_1 = histo3D->GetXaxis()->GetNbins();
    const int nbins_2 = histo3D->GetZaxis()->GetNbins();
    
    double* bin_edges_1 = new double[nbins_1+1];
    double* bin_edges_2 = new double[nbins_2+1];
    
    for(int i=0; i<=nbins_1; i++)
      bin_edges_1[i] = histo3D->GetXaxis()->GetBinLowEdge(i+1);
    for(int i=0; i<=nbins_2; i++)
      bin_edges_2[i] = histo3D->GetZaxis()->GetBinLowEdge(i+1);
    
    TH2F* histo2D = new TH2F("", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
    for(int i=1; i<=nbins_1; i++)
      for(int j=1; j<=nbins_2; j++)
      {
        histo2D -> SetBinContent(i, j, histo3D->GetBinContent(i, bin, j));
        histo2D -> SetBinError(i, j, histo3D->GetBinError(i, bin, j));
      }
      
    histo2D -> SetTitle((std::string(histo3D->GetYaxis()->GetTitle()) + " from " + std::to_string(histo3D->GetYaxis()->GetBinLowEdge(bin)) + " to " + std::to_string(histo3D->GetYaxis()->GetBinLowEdge(bin+1))).c_str());
    histo2D -> GetXaxis() -> SetTitle(histo3D->GetXaxis()->GetTitle());
    histo2D -> GetYaxis() -> SetTitle(histo3D->GetZaxis()->GetTitle());

    return histo2D;
  }
  
  if(axis=="Z" || axis =="z" || axis=="XY" || axis=="xy")
  {
    const int nbins_1 = histo3D->GetXaxis()->GetNbins();
    const int nbins_2 = histo3D->GetYaxis()->GetNbins();
    
    double* bin_edges_1 = new double[nbins_1+1];
    double* bin_edges_2 = new double[nbins_2+1];
    
    for(int i=0; i<=nbins_1; i++)
      bin_edges_1[i] = histo3D->GetXaxis()->GetBinLowEdge(i+1);
    for(int i=0; i<=nbins_2; i++)
      bin_edges_2[i] = histo3D->GetYaxis()->GetBinLowEdge(i+1);
    
    TH2F* histo2D = new TH2F("", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
    for(int i=1; i<=nbins_1; i++)
      for(int j=1; j<=nbins_2; j++)
      {
        histo2D -> SetBinContent(i, j, histo3D->GetBinContent(i, j, bin));
        histo2D -> SetBinError(i, j, histo3D->GetBinError(i, j, bin));
      }
      
    histo2D -> SetTitle((std::string(histo3D->GetZaxis()->GetTitle()) + " from " + std::to_string(histo3D->GetZaxis()->GetBinLowEdge(bin)) + " to " + std::to_string(histo3D->GetZaxis()->GetBinLowEdge(bin+1))).c_str());
    histo2D -> GetXaxis() -> SetTitle(histo3D->GetXaxis()->GetTitle());
    histo2D -> GetYaxis() -> SetTitle(histo3D->GetYaxis()->GetTitle());

    return histo2D;
  }
  
}

void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname)
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}