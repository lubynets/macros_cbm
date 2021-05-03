TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis, int bin);
void SetAxesNames(TH2F* histo, TString xaxisname, TString yaxisname);
std::string StringBinNumber(int number);

void TH3ToTH2ShapeFits()
{
  TFile* fileIn = TFile::Open("shapetensor_fit.apr20.dcmqgsm.nopid.lightcuts1.set4.root");
  
  struct axis
  {
    std::string name_;
    std::string title_;
    std::string id_;
    int another_first_;
    int another_second_;
    int nbins_;
  };
  
  std::vector<axis> axes
  {
    {"centrality", "Centrality, %", "x", 1, 2, -1},
    {"rapidity",   "y_{LAB}",       "y", 0, 2, -1},
    {"pT",         "p_{T}, GeV/c",  "z", 0, 1, -1}
  };
  
  std::vector<std::string> histos{"chi2_bckgr_fit", "chi2_bckgr_func_histo"};
    
  TH3F* hIn;
  
  TH2F* h2Out;
      
  TFile* fileOut = TFile::Open("out.th3toth2chi2bckgrfit.root", "recreate");  
  
  for(auto& histo : histos)
  {
  
    hIn = fileIn -> Get<TH3F>(("h" + histo).c_str());
    
    axes.at(0).nbins_ = hIn->GetXaxis()->GetNbins();
    axes.at(1).nbins_ = hIn->GetYaxis()->GetNbins();
    axes.at(2).nbins_ = hIn->GetZaxis()->GetNbins();
    
    for(auto ax : axes)
    {     
      for(int ibin=1; ibin<=ax.nbins_; ibin++)
      {
        std::string layername = StringBinNumber(ibin);
        std::string dirname = ax.name_ + "/" + layername;
        if(fileOut->GetDirectory(dirname.c_str()) == nullptr)
          fileOut -> mkdir(dirname.c_str());
        fileOut -> cd(dirname.c_str()); 
        
        h2Out   = Get2DLayerFromTH3(hIn  , ax.id_, ibin);
        
        SetAxesNames(h2Out  , axes.at(ax.another_first_).title_, axes.at(ax.another_second_).title_);  
        
        h2Out   -> Write(("h2_" + histo).c_str());        
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