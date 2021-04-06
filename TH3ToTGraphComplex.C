TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2);
void CustomizeGraphsRanges(std::vector<TGraphErrors*> v_graph);
void SetGraphProperties(TGraphErrors* graph, int color, int line_width, int line_style, int marker_size, int marker_style);
void SetAxesNames(TGraphErrors* graph, TString xaxisname, TString yaxisname);
std::string StringBinNumber(int number);


void TH3ToTGraphComplex()
{  
  const int marker_size = 2;
  const int line_width = 2;
  
  TFile* file_fit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.root");
  TFile* file_mcfit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.root");
  TFile* file_mcv1 = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcv1.apr20.dcmqgsm.nopid.lightcuts1.set4.root");
   
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
    std::string name_;
    std::string folder_;
    bool is_mcfitter_;
    bool is_mcv1_;
  };
  
  std::vector<infotype> infotypes
  {
    {"hsignal",        "sgnl",            true,  true},
    {"hbckgr_0",       "bckgr/intercept", true,  false},
    {"hbckgr_1",       "bckgr/slope",     true,  false},
    {"hentries_sgnl",  "Nentries/sgnl",   false, false},
    {"hentries_bckgr", "Nentries/bckgr",  false, false}
  };
  
  TH3F* h_fit;
  TH3F* h_mcfit;
  TH3F* h_mcv1;
  
  TGraphErrors* graph_fit;
  TGraphErrors* graph_mcfit;
  TGraphErrors* graph_mcv1;
      
  TFile* fileOut = TFile::Open("out.th3totgrphcplx.root", "recreate");
  
  for(auto it : infotypes)
  {
  
                        h_fit = file_fit -> Get<TH3F>(("parameters/" + it.name_).c_str());
    if(it.is_mcfitter_) h_mcfit = file_mcfit -> Get<TH3F>(("parameters/" + it.name_).c_str());
    if(it.is_mcv1_)     h_mcv1 = file_mcv1 -> Get<TH3F>("hv1_mc");
    
    axes.at(0).nbins_ = h_fit->GetXaxis()->GetNbins();
    axes.at(1).nbins_ = h_fit->GetYaxis()->GetNbins();
    axes.at(2).nbins_ = h_fit->GetZaxis()->GetNbins();
                
    for(auto ax : axes)
    {
      const int first_nbins = axes.at(ax.another_first_).nbins_;
      const int second_nbins =  axes.at(ax.another_second_).nbins_;
          
      for(int i_first=1; i_first<=first_nbins; i_first++)
        for(int i_second=1; i_second<=second_nbins; i_second++)
        {
          std::string binname = axes.at(ax.another_first_).letter_ + StringBinNumber(i_first) + "_" + axes.at(ax.another_second_).letter_ + StringBinNumber(i_second);
          std::string dirname = it.folder_ +"/" + ax.name_ + "/" + binname;
          fileOut -> mkdir(dirname.c_str());
          fileOut -> cd(dirname.c_str());
          
                              graph_fit = GetGraph(h_fit, ax.id_, i_first, i_second);
          if(it.is_mcfitter_) graph_mcfit = GetGraph(h_mcfit, ax.id_, i_first, i_second);
          if(it.is_mcv1_)     graph_mcv1 = GetGraph(h_mcv1, ax.id_, i_first, i_second);
          
                              SetAxesNames(graph_fit, ax.title_, "v_{1x}");
          if(it.is_mcfitter_) SetAxesNames(graph_mcfit, ax.title_, "v_{1x}");
          if(it.is_mcv1_)     SetAxesNames(graph_mcv1, ax.title_, "v_{1x}"); 
          
                              SetGraphProperties(graph_fit,   kBlue,     line_width, 1, marker_size, 8);
          if(it.is_mcfitter_) SetGraphProperties(graph_mcfit, kGreen+2,  line_width, 1, marker_size, 8);
          if(it.is_mcv1_)     SetGraphProperties(graph_mcv1,  kRed,      line_width, 1, marker_size, 8);
          
                              graph_fit -> Write("fit");
          if(it.is_mcfitter_) graph_mcfit -> Write("mcfit");
          if(it.is_mcv1_)     graph_mcv1 -> Write("mcv1");
          
          if(! it.is_mcfitter_) continue;
          
          TCanvas cc("canvas", "canvas", 1500, 900);
          cc.cd();
          
                          graph_fit -> Draw("APL");
                          graph_mcfit -> Draw("PL");
          if(it.is_mcv1_) graph_mcv1 -> Draw("PL");
            
          TLegend* leg = new TLegend(0.35, 0.68, 0.55, 0.88);
          leg -> SetBorderSize(0);
          if(it.is_mcv1_) leg -> AddEntry(graph_mcv1,  "MC",                "PL");
                          leg -> AddEntry(graph_mcfit, "RECO, MC-match",    "PL");
                          leg -> AddEntry(graph_fit,   "RECO, invmass-fit", "PL");
          leg -> Draw("same");
          
          cc.Write();
            
        }
    }
  
  }
   
  fileOut -> Close();
}

TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2)
{
  if(axis=="x" || axis=="X")
  {
    TGraphErrors* graph = new TGraphErrors();
    for(int i=1; i<=histo3D->GetXaxis()->GetNbins(); i++)
    {
      graph -> SetPoint(i-1, histo3D->GetXaxis()->GetBinCenter(i), histo3D->GetBinContent(i, bin1, bin2));
      graph -> SetPointError(i-1, 0., histo3D->GetBinError(i, bin1, bin2));
    }
    
    return graph;
  }
  
  if(axis=="y" || axis=="Y")
  {
    TGraphErrors* graph = new TGraphErrors();
    for(int i=1; i<=histo3D->GetYaxis()->GetNbins(); i++)
    {
      graph -> SetPoint(i-1, histo3D->GetYaxis()->GetBinCenter(i), histo3D->GetBinContent(bin1, i, bin2));
      graph -> SetPointError(i-1, 0., histo3D->GetBinError(bin1, i, bin2));
    }
    
    return graph;
  }
  
  if(axis=="z" || axis=="Z")
  {
    TGraphErrors* graph = new TGraphErrors();
    for(int i=1; i<=histo3D->GetZaxis()->GetNbins(); i++)
    {
      graph -> SetPoint(i-1, histo3D->GetZaxis()->GetBinCenter(i), histo3D->GetBinContent(bin1, bin2, i));
      graph -> SetPointError(i-1, 0., histo3D->GetBinError(bin1, bin2, i));
    }
    
    return graph;
  }  
}

void SetGraphProperties(TGraphErrors* graph, int color, int line_width, int line_style, int marker_size, int marker_style)
{
  graph -> SetLineColor(color);
  graph -> SetMarkerColor(color);
  graph -> SetLineWidth(line_width);
  graph -> SetLineStyle(line_style);
  graph -> SetMarkerStyle(marker_style);
  graph -> SetMarkerSize(marker_size);
  float y_min = graph -> GetYaxis() -> GetXmin();
  float y_max = graph -> GetYaxis() -> GetXmax();
  
  if(y_min>0)
    y_min = 0.3*y_min;
  else
    y_min = 1.2*y_min;
  
  if(y_max>0)
    y_max = 1.2*y_max;
  else
    y_max=0.3*y_max;
  
  graph -> GetYaxis() -> SetRangeUser(y_min, y_max);  
}

void SetAxesNames(TGraphErrors* graph, TString xaxisname, TString yaxisname)
{
  graph -> GetXaxis() -> SetTitle(xaxisname);
  graph -> GetYaxis() -> SetTitle(yaxisname);
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}

void CustomizeGraphsRanges(std::vector<TGraphErrors*> v_graph)
{
  float min =  999.;
  float max = -999.;
  
  for(auto g : v_graph)
  {
    min = std::min(min, (float)g->GetYaxis()->GetXmin());
    max = std::max(max, (float)g->GetYaxis()->GetXmax());
  }
  
  if(min>0)
    min = 0.3*min;
  else
    min = 1.2*min;
  
  if(max>0)
    max = 1.2*max;
  else
    max=0.3*max;
  
  for(auto g : v_graph)
    g->GetYaxis()->SetRangeUser(min, max);
  
}