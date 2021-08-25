TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2, float shift);
std::string GetGraphCell(TH3F* histo3D, std::string axis, std::string letter1, std::string letter2, int bin1, int bin2);
void CustomizeGraphsRanges(std::vector<TGraphErrors*> v_graph);
void SetGraphProperties(TGraphErrors* graph, int color, int line_width, int line_style, int marker_size, int marker_style);
void SetAxesNames(TGraphErrors* graph, TString xaxisname, TString yaxisname);
std::string StringBinNumber(int number);
void ResetSystErrors(TGraphErrors* graph_AVE, const TGraphErrors* graph_XX, const TGraphErrors* graph_YY, float left_edge);

void TH3ToTGraphVeryComplex()
{  
  const int marker_size = 2;
  const int line_width = 2;
  
  TFile* file_fit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.AVE.root");
  TFile* file_mcfit_withcorr_AVE = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.AVE.root");
  TFile* file_mcfit_withcorr_XX = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.root");
  TFile* file_mcfit_withcorr_YY = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.YY.root");
  TFile* file_mcfit_wocorr_XX = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.now.XX.root");
  TFile* file_mcfit_wocorr_YY = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcfitter.apr20.dcmqgsm.nopid.lightcuts1.set4.now.YY.root");
  TFile* file_mcv1 = TFile::Open("/home/user/cbmdir/working/qna/fits/out.mcv1.apr20.dcmqgsm.nopid.lightcuts1.set4.AVE.root");
     
  struct axis
  {
    std::string name_;
    std::string letter_;
    std::string title_;
    std::string id_;
    int another_first_;
    int another_second_;
    int nbins_;
    float shift_;
    float left_edge_;
    float right_edge_;
  };
  
  std::vector<axis> axes
  {
    {"centrality", "C",  "Centrality, %", "x", 1, 2, -1, 0,        0,    100},
    {"rapidity",   "y",  "y_{CM}",        "y", 0, 2, -1, -1.62179, -0.6, 1.0},
    {"pT",         "pT", "p_{T}, GeV/c",  "z", 0, 1, -1, 0,        0.2,  1.4}
  };
  
  TH3F* h_fit;
  TH3F* h_mcfit_withcorr_AVE;
  TH3F* h_mcfit_withcorr_XX;
  TH3F* h_mcfit_withcorr_YY;
  TH3F* h_mcfit_wocorr_XX;
  TH3F* h_mcfit_wocorr_YY;
  TH3F* h_mcv1;
  
  TGraphErrors* graph_fit;
  TGraphErrors* graph_mcfit_withcorr_AVE;
  TGraphErrors* graph_mcfit_withcorr_AVE_clone;
  TGraphErrors* graph_mcfit_withcorr_XX;
  TGraphErrors* graph_mcfit_withcorr_YY;
  TGraphErrors* graph_mcfit_wocorr_XX;
  TGraphErrors* graph_mcfit_wocorr_YY;
  TGraphErrors* graph_mcv1;
      
  TFile* fileOut = TFile::Open("out.tgraph.verycplx.root", "recreate");
  
  bool is_first_canvas = true;
  
  h_fit                = file_fit                -> Get<TH3F>("parameters/hsignal");
  h_mcfit_withcorr_AVE = file_mcfit_withcorr_AVE -> Get<TH3F>("parameters/hsignal");
  h_mcfit_withcorr_XX  = file_mcfit_withcorr_XX  -> Get<TH3F>("parameters/hsignal");
  h_mcfit_withcorr_YY  = file_mcfit_withcorr_YY  -> Get<TH3F>("parameters/hsignal");
  h_mcfit_wocorr_XX    = file_mcfit_wocorr_XX    -> Get<TH3F>("parameters/hsignal");
  h_mcfit_wocorr_YY    = file_mcfit_wocorr_YY    -> Get<TH3F>("parameters/hsignal");
  h_mcv1               = file_mcv1               -> Get<TH3F>("hv1_mc");
  
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
        std::string dirname = "sgnl/" + ax.name_ + "/" + binname;
        fileOut -> mkdir(dirname.c_str());
        fileOut -> cd(dirname.c_str());
        
        graph_fit                = GetGraph(h_fit               , ax.id_, i_first, i_second, ax.shift_);
        graph_mcfit_withcorr_AVE = GetGraph(h_mcfit_withcorr_AVE, ax.id_, i_first, i_second, ax.shift_);
        graph_mcfit_withcorr_XX  = GetGraph(h_mcfit_withcorr_XX , ax.id_, i_first, i_second, ax.shift_);
        graph_mcfit_withcorr_YY  = GetGraph(h_mcfit_withcorr_YY , ax.id_, i_first, i_second, ax.shift_);
        graph_mcfit_wocorr_XX    = GetGraph(h_mcfit_wocorr_XX   , ax.id_, i_first, i_second, ax.shift_);
        graph_mcfit_wocorr_YY    = GetGraph(h_mcfit_wocorr_YY   , ax.id_, i_first, i_second, ax.shift_);
        graph_mcv1               = GetGraph(h_mcv1              , ax.id_, i_first, i_second, ax.shift_);
        
        graph_mcfit_withcorr_AVE_clone = (TGraphErrors*)graph_mcfit_withcorr_AVE -> Clone();
        ResetSystErrors(graph_mcfit_withcorr_AVE_clone, graph_mcfit_withcorr_XX, graph_mcfit_withcorr_YY, ax.left_edge_);
        
        graph_mcv1 -> SetTitle(("hsignal, " + GetGraphCell(h_fit, ax.id_, axes.at(ax.another_first_).letter_, axes.at(ax.another_second_).letter_, i_first, i_second)).c_str());
        
        SetAxesNames(graph_mcv1,  ax.title_, "v_{1}"); 
        
        SetGraphProperties(graph_fit               , kBlack,    line_width, 1, marker_size, 8);
        SetGraphProperties(graph_mcfit_withcorr_AVE, kRed,      line_width, 1, marker_size, 8);
        SetGraphProperties(graph_mcfit_wocorr_XX   , kBlue,     line_width, 1, marker_size, 4);
        SetGraphProperties(graph_mcfit_wocorr_YY   , kGreen+2,  line_width, 1, marker_size, 4);
        SetGraphProperties(graph_mcv1              , kRed,      line_width, 1, marker_size, 8);
        
        graph_mcfit_withcorr_AVE_clone -> SetFillStyle(3001);
        graph_mcfit_withcorr_AVE_clone -> SetFillColor(kRed-4);
                  
        TCanvas cc("canvas", "canvas", 1500, 900);
        cc.cd();
        
        CustomizeGraphsRanges({graph_fit, graph_mcfit_withcorr_AVE, graph_mcfit_withcorr_AVE_clone, graph_mcfit_wocorr_XX, graph_mcfit_wocorr_YY, graph_mcv1});
        graph_mcv1->GetXaxis()->SetLimits(ax.left_edge_, ax.right_edge_);
        
        graph_mcv1 -> Draw("ALX");
        graph_mcfit_wocorr_XX -> Draw("P");
        graph_mcfit_wocorr_YY-> Draw("P");
        graph_mcfit_withcorr_AVE -> Draw("P");
        graph_mcfit_withcorr_AVE_clone -> Draw("E2");
        graph_fit -> Draw("P");
        
          
//         TLegend* leg = new TLegend(0.85, 0.85, 1, 1);
// //           leg -> SetBorderSize(0);
//         leg -> AddEntry(graph_mcv1,  "MC",                "PL");
//         leg -> AddEntry(graph_mcfit, "RECO, MC-match",    "PL");
//         leg -> AddEntry(graph_fit,   "RECO, invmass-fit", "PL");
//         leg -> Draw("same");
        
        cc.Write();
        if(is_first_canvas)
          cc.Print("out.tgraph.verycplx.pdf(", "pdf");
        else
          cc.Print("out.tgraph.verycplx.pdf", "pdf");

        is_first_canvas = false;
      }
  }
  
  
  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("out.tgraph.verycplx.pdf)", "pdf");
   
  fileOut -> Close();
}

TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2, float shift)
{
  if(axis=="x" || axis=="X")
  {
    TGraphErrors* graph = new TGraphErrors();
    for(int i=1; i<=histo3D->GetXaxis()->GetNbins(); i++)
    {
      graph -> SetPoint(i-1, histo3D->GetXaxis()->GetBinCenter(i) + shift, histo3D->GetBinContent(i, bin1, bin2));
      graph -> SetPointError(i-1, 0., histo3D->GetBinError(i, bin1, bin2));
    }
    
    return graph;
  }
  
  if(axis=="y" || axis=="Y")
  {
    TGraphErrors* graph = new TGraphErrors();
    for(int i=1; i<=histo3D->GetYaxis()->GetNbins(); i++)
    {
      graph -> SetPoint(i-1, histo3D->GetYaxis()->GetBinCenter(i) + shift, histo3D->GetBinContent(bin1, i, bin2));
      graph -> SetPointError(i-1, 0., histo3D->GetBinError(bin1, i, bin2));
    }
    
    return graph;
  }
  
  if(axis=="z" || axis=="Z")
  {
    TGraphErrors* graph = new TGraphErrors();
    for(int i=1; i<=histo3D->GetZaxis()->GetNbins(); i++)
    {
      graph -> SetPoint(i-1, histo3D->GetZaxis()->GetBinCenter(i) + shift, histo3D->GetBinContent(bin1, bin2, i));
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

std::string GetGraphCell(TH3F* histo3D, std::string axis, std::string letter1, std::string letter2, int bin1, int bin2)
{
  float lo1;
  float lo2;
  float hi1;
  float hi2;
  
  if(axis=="x" || axis=="X")
  {
    lo1 = histo3D -> GetYaxis() -> GetBinLowEdge(bin1);
    hi1 = histo3D -> GetYaxis() -> GetBinUpEdge(bin1);
    lo2 = histo3D -> GetZaxis() -> GetBinLowEdge(bin2);
    hi2 = histo3D -> GetZaxis() -> GetBinUpEdge(bin2);
  }
  
  if(axis=="y" || axis=="Y")
  {
    lo1 = histo3D -> GetXaxis() -> GetBinLowEdge(bin1);
    hi1 = histo3D -> GetXaxis() -> GetBinUpEdge(bin1);
    lo2 = histo3D -> GetZaxis() -> GetBinLowEdge(bin2);
    hi2 = histo3D -> GetZaxis() -> GetBinUpEdge(bin2);
  }
  
  if(axis=="z" || axis=="Z")
  {
    lo1 = histo3D -> GetXaxis() -> GetBinLowEdge(bin1);
    hi1 = histo3D -> GetXaxis() -> GetBinUpEdge(bin1);
    lo2 = histo3D -> GetYaxis() -> GetBinLowEdge(bin2);
    hi2 = histo3D -> GetYaxis() -> GetBinUpEdge(bin2);
  }
  
  std::string cell = letter1 + "#in[" +
                     std::to_string(lo1).substr(0, std::to_string(lo1).find(".")+3) + ", " +
                     std::to_string(hi1).substr(0, std::to_string(hi1).find(".")+3) + "], " +
                     letter2 + "#in[" +
                     std::to_string(lo2).substr(0, std::to_string(lo2).find(".")+3) + ", " +
                     std::to_string(hi2).substr(0, std::to_string(hi2).find(".")+3) + "]";
                                          
  return cell;
}

void ResetSystErrors(TGraphErrors* graph_AVE, const TGraphErrors* graph_XX, const TGraphErrors* graph_YY, float left_edge)
{
  const int N = graph_AVE -> GetN();
  float left_bin_edge = left_edge;
  
  for(int i=0; i<N; i++)
  {
    const float x_syst_err = graph_AVE -> GetPointX(i) - left_bin_edge;
    left_bin_edge += 2*x_syst_err;
    
    const float y_XX = graph_XX -> GetPointY(i);
    const float y_YY = graph_YY -> GetPointY(i);
    const float y_syst_err = std::fabs(y_YY - y_XX)/2;
    
    graph_AVE -> SetPointError(i, x_syst_err, y_syst_err);
  }
}