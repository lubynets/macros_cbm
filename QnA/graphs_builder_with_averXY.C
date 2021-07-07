void SetGraphProperties(TGraphAsymmErrors* graph, int color, int line_width, int line_style, int marker_size, int marker_style)
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

void ManageGraph(TGraphAsymmErrors* graph, std::string name, std::string title, std::string x_title)
{
  graph -> SetName(name.c_str());
  graph -> SetTitle(title.c_str());
  graph -> GetXaxis() -> SetTitle(x_title.c_str());
  graph -> GetYaxis() -> SetTitle("v_{1}");
  graph -> Write();
}

void MultGrapgh(TGraphAsymmErrors* graph, double value)
{
  const int N = graph -> GetN();
  double* y = graph -> GetY();
  double* ey = graph -> GetEY();
  
  for(int i=0; i<N; i++)
  {
    y[i] = value*y[i];
    ey[i] = value*ey[i];
  }  
}

void ShiftGraph(TGraphAsymmErrors* graph, double value)
{
  const int N = graph -> GetN();
  double* x = graph -> GetX();
  
  for(int i=0; i<N; i++)
    x[i] = x[i] + value;
}

TGraphAsymmErrors* AverGraph(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2)
{
  const int N1 = graph1 -> GetN();
  const int N2 = graph2 -> GetN();
  
  if(N1 != N2)
  {
    std::cout << "Error at AverGraph(): Graphs are not compatible!\n";
    assert(false);
  }

  TGraphAsymmErrors* graph = (TGraphAsymmErrors*) graph1->Clone();
  double* y1 = graph -> GetY();
  double* ey1l = graph -> GetEYlow();
  double* ey1h = graph -> GetEYhigh();
  double* y2 = graph2 -> GetY();
  double* ey2l = graph2 -> GetEYlow();
  double* ey2h = graph2 -> GetEYhigh();

  for(int i=0; i<N1; i++)
  {
    y1[i] = (y1[i]+y2[i])/2.;
    ey1l[i] = std::sqrt(ey1l[i]*ey1l[i] + ey2l[i]*ey2l[i])/2.;
    ey1h[i] = std::sqrt(ey1h[i]*ey1h[i] + ey2h[i]*ey2h[i])/2.;
  }  

  return graph;
}

TGraphAsymmErrors* AverGraphErrorAsDiff(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2)
{
// Error is calculated as difference between two summing values  
  const int N1 = graph1 -> GetN();
  const int N2 = graph2 -> GetN();
  
  if(N1 != N2)
  {
    std::cout << "Error at AverGraph(): Graphs are not compatible!\n";
    assert(false);
  }

  TGraphAsymmErrors* graph = (TGraphAsymmErrors*) graph1->Clone();
  double* y1 = graph -> GetY();
  double* ey1l = graph -> GetEYlow();
  double* ey1h = graph -> GetEYhigh();
  double* y2 = graph2 -> GetY();

  for(int i=0; i<N1; i++)
  {
    y1[i] = (y1[i]+y2[i])/2.;
    ey1l[i] = fabs(y1[i]-y2[i]);
    ey1h[i] = fabs(y1[i]-y2[i]);
  }  

  return graph;
}

void graphs_builder_with_averXY(const TString& fileName)
{
  const int marker_size = 2;
  const int line_width = 2;
  
  struct axis
  {
    std::string name_;
    std::string sim_name_;
    std::string reco_name_;
    std::string title_;
  };
  
  std::vector<axis> axes
  {
//     {"pT", "SimTracks_pT", "RecParticlesMcPid_pT", "p_{T}, GeV/c"},
    {"rapidity", "SimTracks_rapidity", "RecParticlesMcPid_rapidity", "y_{CM}"},
//     {"centrality", "Centrality", "Centrality", "Centrality, %"}
  };
  
  enum eComponent : short
  {
    kXX = 0,
    kYY,
    kNumberOfComponents
  };
    
  enum eStep : short
  {
    kMc = 0,
    kPlain,
    kRecentered,
    kTwist,
    kRescaled,
    kNumberOfSteps
  };
  
  std::vector<std::vector<Qn::DataContainer<Qn::Stats,Qn::Axis<double>>*>> uQ;
  uQ.resize(kNumberOfComponents);
  for(int iComponent=0; iComponent<kNumberOfComponents; iComponent++)
    uQ.at(iComponent).resize(kNumberOfSteps);
   
  std::vector<std::vector<TGraphAsymmErrors*>> graph;
  graph.resize(kNumberOfComponents);
  for(int iComponent=0; iComponent<kNumberOfComponents; iComponent++)
    graph.at(iComponent).resize(kNumberOfSteps);
  
  std::vector<std::string> component_name {"Q1x_Q1x", "Q1y_Q1y"};
  
  struct step
  {
    std::string name_;
    std::string path_;
  };
  
  std::vector<step> steps
  {
    {"mc",          "mc_lambda_PLAIN"},
    {"plain",       "lambda_PLAIN"},
    {"recentered",  "lambda_RECENTERED"},
    {"twist",       "lambda_TWIST"},
    {"rescaled",    "lambda_RESCALED"}
  };
   
  TFile* v1file = TFile::Open(fileName);
  
  for(int iComponent=0; iComponent<kNumberOfComponents; iComponent++)
    for(int iStep=0; iStep<kNumberOfSteps; iStep++)
      uQ.at(iComponent).at(iStep) = (Qn::DataContainer<Qn::Stats,Qn::Axis<double>>*)v1file -> Get((steps.at(iStep).path_ + "_psi_PLAIN_" + component_name.at(iComponent)).c_str());

  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  for(auto ax : axes)
  {
    TDirectory* dir = fileOut->mkdir(ax.name_.c_str());
    dir -> cd();
    
    for(int iStep=0; iStep<kNumberOfSteps; iStep++)
    {
      for(int iComponent=0; iComponent<kNumberOfComponents; iComponent++)
      {
        std::string proj_name;
        if(iStep == kMc)
          proj_name = ax.sim_name_;
        else
          proj_name = ax.reco_name_;
        auto projection = uQ.at(iComponent).at(iStep) -> Projection({proj_name.c_str()});
        graph.at(iComponent).at(iStep) = Qn::DataContainerHelper::ToTGraph(projection);
        const double y_beam = 1.62179;
        if(ax.name_ == "rapidity")
          ShiftGraph(graph.at(iComponent).at(iStep), -y_beam);
        
        if(ax.name_ == "centrality")
          graph.at(iComponent).at(iStep) -> RemovePoint(graph.at(iComponent).at(iStep)->GetN()-1);
        
        ManageGraph(graph.at(iComponent).at(iStep), (steps.at(iStep).name_ + "_" + component_name.at(iComponent)).c_str(), (steps.at(iStep).name_ + "_" + component_name.at(iComponent)).c_str(), ax.title_);
      }
      TGraphAsymmErrors* graph_average = AverGraph(graph.at(kXX).at(iStep), graph.at(kYY).at(iStep));
      ManageGraph(graph_average, (steps.at(iStep).name_ + "_average_Q1x_Q1y").c_str(), (steps.at(iStep).name_ + "_average_Q1x_Q1y").c_str(), ax.title_);
      
      TGraphAsymmErrors* graph_average_error_as_diff = AverGraphErrorAsDiff(graph.at(kXX).at(iStep), graph.at(kYY).at(iStep));
      ManageGraph(graph_average_error_as_diff, (steps.at(iStep).name_ + "_average_Q1x_Q1y_syserr").c_str(), (steps.at(iStep).name_ + "_average_Q1x_Q1y_syserr").c_str(), ax.title_);
    }
      
    //---------------- canvas 2 ------------------------------------------------------------------
    TCanvas c2("XX", "XX", 800, 600);
    c2.cd();
    
    SetGraphProperties(graph.at(kXX).at(kMc),          kRed,        line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kXX).at(kPlain),       kGreen+2,    line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kXX).at(kRecentered),  kMagenta,    line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kXX).at(kTwist),       kBlue,       line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kXX).at(kRescaled),    kBlack,      line_width, 1, marker_size, 8);
    
    graph.at(kXX).at(kMc)         -> Draw("APL");
    graph.at(kXX).at(kPlain)      -> Draw("PL");
    graph.at(kXX).at(kRecentered) -> Draw("PL");
    graph.at(kXX).at(kTwist)      -> Draw("PL");
    graph.at(kXX).at(kRescaled)   -> Draw("PL");
    
    TLegend* leg2 = new TLegend(0.35, 0.58, 0.55, 0.88);
    leg2 -> SetBorderSize(0);
    leg2 -> AddEntry(graph.at(kXX).at(kMc),         "mc, XX",         "PL");
    leg2 -> AddEntry(graph.at(kXX).at(kPlain),      "plain, XX",      "PL");
    leg2 -> AddEntry(graph.at(kXX).at(kRecentered), "recentered, XX", "PL");
    leg2 -> AddEntry(graph.at(kXX).at(kTwist),      "twist, XX",      "PL");
    leg2 -> AddEntry(graph.at(kXX).at(kRescaled),   "rescaled, XX",   "PL");
    leg2 -> Draw("same");
    
    c2.Write();
    //--------------------------------------------------------------------------------------------
    
    //---------------- canvas 3 ------------------------------------------------------------------
    TCanvas c3("YY", "YY", 800, 600);
    c3.cd();
    
    SetGraphProperties(graph.at(kYY).at(kMc),          kRed,        line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kYY).at(kPlain),       kGreen+2,    line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kYY).at(kRecentered),  kMagenta,    line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kYY).at(kTwist),       kBlue,       line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kYY).at(kRescaled),    kBlack,      line_width, 1, marker_size, 8);
    
    graph.at(kYY).at(kMc)         -> Draw("APL");
    graph.at(kYY).at(kPlain)      -> Draw("PL");
    graph.at(kYY).at(kRecentered) -> Draw("PL");
    graph.at(kYY).at(kTwist)      -> Draw("PL");
    graph.at(kYY).at(kRescaled)   -> Draw("PL");
    
    TLegend* leg3 = new TLegend(0.35, 0.58, 0.55, 0.88);
    leg3 -> SetBorderSize(0);
    leg3 -> AddEntry(graph.at(kYY).at(kMc),         "mc, YY",         "PL");
    leg3 -> AddEntry(graph.at(kYY).at(kPlain),      "plain, YY",      "PL");
    leg3 -> AddEntry(graph.at(kYY).at(kRecentered), "recentered, YY", "PL");
    leg3 -> AddEntry(graph.at(kYY).at(kTwist),      "twist, YY",      "PL");
    leg3 -> AddEntry(graph.at(kYY).at(kRescaled),   "rescaled, YY",   "PL");
    leg3 -> Draw("same");
    
    c3.Write();
    //--------------------------------------------------------------------------------------------
    
    //-------------- canvas 4 --------------------------------------------------------------------
    TCanvas c4("XX_and_YY", "XX_and_YY", 800, 600);
    c4.cd();
    
    SetGraphProperties(graph.at(kXX).at(kMc),       kRed,   line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kYY).at(kMc),       kRed,   line_width, 2, marker_size, 4);
    SetGraphProperties(graph.at(kXX).at(kRescaled), kBlack, line_width, 1, marker_size, 8);
    SetGraphProperties(graph.at(kYY).at(kRescaled), kBlack, line_width, 2, marker_size, 4);
    
    graph.at(kXX).at(kMc)       -> Draw("APL");
    graph.at(kYY).at(kMc)       -> Draw("PL");
    graph.at(kXX).at(kRescaled) -> Draw("PL");
    graph.at(kYY).at(kRescaled) -> Draw("PL");
    
    TLegend* leg4 = new TLegend(0.35, 0.64, 0.55, 0.88);
    leg4 -> SetBorderSize(0);
    leg4 -> AddEntry(graph.at(kXX).at(kMc),       "mc, XX",       "PL");
    leg4 -> AddEntry(graph.at(kYY).at(kMc),       "mc, YY",       "PL");
    leg4 -> AddEntry(graph.at(kXX).at(kRescaled), "rescaled, XX", "PL");
    leg4 -> AddEntry(graph.at(kYY).at(kRescaled), "rescaled, YY", "PL");
    leg4 -> Draw("same");
    
    c4.Write();
    //--------------------------------------------------------------------------------------------
  }
  
  fileOut -> Close();
}
    