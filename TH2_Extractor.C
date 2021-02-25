TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis="Z", int bin=1);
std::vector<TH2F*> Get2DFromTH3(TH3F* histo3D, TString axis="Z");

// void TH2_Extractor()
// {
//   TFile* file_fit = TFile::Open("/home/user/cbmdir/QnAnalysis/build/src/QnAnalysisDiscriminator/out.fitter.root");
//   std::vector<std::string> histonames{"hsignal", "hbckgr_0", "hbckgr_1"};
//   
//   TFile* fileOut = TFile::Open("out.layer.root", "recreate");
//   
//   for(auto& histoname : histonames)
//   {
//     TH3F* hfit = file_fit -> Get<TH3F>(("parameters/"+histoname).c_str());
//     std::vector<TH2F*> vec = Get2DFromTH3(hfit, "X");
//     
//     int i=1;
//     
//     for(auto& v : vec)
//     {
//       v->SetName((histoname + "_layer" + std::to_string(i)).c_str());
//       v->Write();
//       i++;
//     }
//   }
//   
//   fileOut -> Close();  
// }

void TH2_Extractor()
{
  TFile* fileIn = TFile::Open("out.ratio.set2.root");
  
  std::vector<std::string> histotypes{"ratio", "diff_from_one"};
  std::vector<std::string> histonames{"hsignal", "hbckgr_0", "hbckgr_1"};
  
  TFile* fileOut = TFile::Open("out.ratio_layer.root", "recreate");
  TDirectory* dirRatio = fileOut->mkdir("ratios");
  TDirectory* dirDiff = fileOut->mkdir("diff_from_ones");
  
  for(auto& histotype : histotypes)
  {
    if(histotype == "ratio")
      dirRatio->cd();
    else
      dirDiff->cd();
    for(auto& histoname : histonames)
    {
      TH3F* h = fileIn -> Get<TH3F>((histotype + "s/" + histoname + "_" + histotype).c_str());
      std::cout << (histotype + "s/" + histoname + "_" + histotype).c_str() << "\n";
      std::vector<TH2F*> vec = Get2DFromTH3(h, "X");
      
      int i=1;
      
      for(auto& v : vec)
      {
        v->SetName((histoname + "_layer" + std::to_string(i)).c_str());
        v->Write();
        i++;
      }
    }
  }
  
  fileOut -> Close();  
}

std::vector<TH2F*> Get2DFromTH3(TH3F* histo3D, TString axis="Z")
{
  int nbins;
  if(axis=="X" || axis =="x" || axis=="YZ" || axis=="yz")
    nbins = histo3D->GetXaxis()->GetNbins();
  else if(axis=="Y" || axis =="y" || axis=="XZ" || axis=="xz")
    nbins = histo3D->GetYaxis()->GetNbins();
  else if(axis=="Z" || axis =="z" || axis=="XY" || axis=="xy")
    nbins = histo3D->GetZaxis()->GetNbins();
  else
    throw std::runtime_error("No axis with such name");
  
  std::vector<TH2F*> histo2D_vector;
    
  for(int i=1; i<=nbins; i++)
    histo2D_vector.push_back(Get2DLayerFromTH3(histo3D, axis, i));
  
  return histo2D_vector;
}

TH2F* Get2DLayerFromTH3(TH3F* histo3D, TString axis="Z", int bin=1)
{
//   auto getaxis = [histo3D] (TString ax) {
//       if ax == "X" return histo -> GetXax();
//       ...
//       ...
//   }
//   
//   
//   auto ax_ptr = getaxis("X")
//   
//   "xy"
  
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
    
    TH2F* histo2D = new TH2F("histo2D", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
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
    
    TH2F* histo2D = new TH2F("histo2D", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
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
    
    TH2F* histo2D = new TH2F("histo2D", "", nbins_1, bin_edges_1, nbins_2, bin_edges_2);
    
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