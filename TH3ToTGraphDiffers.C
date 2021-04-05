TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2);
void CustomizeGraphsRanges(std::vector<TGraphErrors*> v_graph);
void SetGraphProperties(TGraphErrors* graph, int color, int line_width, int line_style, int marker_size, int marker_style);
void SetAxesNames(TGraphErrors* graph, TString xaxisname, TString yaxisname);
std::string StringBinNumber(int number);


void TH3ToTGraphDiffers()
{
  const int marker_size = 2;
  const int line_width = 2;
  
  TFile* file_fit2mcfit = TFile::Open("/home/user/cbmdir/working/qna/fits/dirivatives/out.difference_fitter2mcfitter.set4.root");
  TFile* file_fit2mcv1 = TFile::Open("/home/user/cbmdir/working/qna/fits/dirivatives/out.difference_fitter2mcv1.set4.root");
  
  TH3F* h_fit2mcfit = file_fit2mcfit -> Get<TH3F>("hsignal_diff_interms_err");
  TH3F* h_fit2mcv1 = file_fit2mcv1 -> Get<TH3F>("hsignal_diff_interms_err");
  
  const int C_nbins = h_fit2mcfit->GetXaxis()->GetNbins();
  const int y_nbins = h_fit2mcfit->GetYaxis()->GetNbins();
  const int pT_nbins = h_fit2mcfit->GetZaxis()->GetNbins();
  
  TFile* fileOut = TFile::Open("out.th3totgrphdiffers.root", "recreate");
    
  TGraphErrors* graph_fit2mcfit;
  TGraphErrors* graph_fit2mcv1;
  
  //************ centrality ******************************************
  fileOut->mkdir("centrality/graphs");
  fileOut->mkdir("centrality/canvas");
  
  for(int i_y=1; i_y<=y_nbins; i_y++)
    for(int i_pT=1; i_pT<=pT_nbins; i_pT++)
    {
      graph_fit2mcfit = GetGraph(h_fit2mcfit, "x", i_y, i_pT);
      graph_fit2mcv1 = GetGraph(h_fit2mcv1, "x", i_y, i_pT);
      
      SetAxesNames(graph_fit2mcfit, "Centrality, %", "\" #chi^{2} \"");
      SetAxesNames(graph_fit2mcv1, "Centrality, %", "\" #chi^{2} \"");
      
      SetGraphProperties(graph_fit2mcfit,   kBlue, line_width, 1, marker_size, 8);
      SetGraphProperties(graph_fit2mcv1,    kRed,  line_width, 1, marker_size, 8);
      
      TString graphname = "y" + StringBinNumber(i_y) + "_pT" + StringBinNumber(i_pT);
      
      fileOut -> cd("centrality/graphs");
      graph_fit2mcfit -> Write(graphname + "_fit_mcfit");
      graph_fit2mcv1 -> Write(graphname + "_fit_mcv1");
      
      TCanvas cc(graphname, graphname, 1500, 900);
      cc.cd();
      
      CustomizeGraphsRanges({graph_fit2mcfit, graph_fit2mcv1});
      graph_fit2mcfit -> Draw("APL");
      graph_fit2mcv1 -> Draw("PL");
      
      TLegend* leg = new TLegend(0.25, 0.75, 0.65, 0.88);
      leg -> SetBorderSize(0);
      leg -> AddEntry(graph_fit2mcfit, "invmass fit - RECO with MC-match", "PL");
      leg -> AddEntry(graph_fit2mcv1,  "invmass fit - MC",                 "PL");
      leg -> Draw("same");
      
      fileOut -> cd("centrality/canvas");
      cc.Write();
    }
  //******************************************************************
  
  //************ rapidity ******************************************
  fileOut->mkdir("rapidity/graphs");
  fileOut->mkdir("rapidity/canvas");
  
  for(int i_C=1; i_C<=C_nbins; i_C++)
    for(int i_pT=1; i_pT<=pT_nbins; i_pT++)
    {
      graph_fit2mcfit = GetGraph(h_fit2mcfit, "y", i_C, i_pT);
      graph_fit2mcv1 = GetGraph(h_fit2mcv1, "y", i_C, i_pT);
      
      SetAxesNames(graph_fit2mcfit, "y_{LAB}", "\" #chi^{2} \"");
      SetAxesNames(graph_fit2mcv1, "y_{LAB}", "\" #chi^{2} \"");
      
      SetGraphProperties(graph_fit2mcfit,   kBlue, line_width, 1, marker_size, 8);
      SetGraphProperties(graph_fit2mcv1,    kRed,  line_width, 1, marker_size, 8);
      
      TString graphname = "C" + StringBinNumber(i_C) + "_pT" + StringBinNumber(i_pT);
      
      fileOut -> cd("rapidity/graphs");
      graph_fit2mcfit -> Write(graphname + "_fit_mcfit");
      graph_fit2mcv1 -> Write(graphname + "_fit_mcv1");
      
      TCanvas cc(graphname, graphname, 1500, 900);
      cc.cd();
      
      CustomizeGraphsRanges({graph_fit2mcfit, graph_fit2mcv1});
      graph_fit2mcfit -> Draw("APL");
      graph_fit2mcv1 -> Draw("PL");
      
      TLegend* leg = new TLegend(0.25, 0.75, 0.65, 0.88);
      leg -> SetBorderSize(0);
      leg -> AddEntry(graph_fit2mcfit, "invmass fit - RECO with MC-match", "PL");
      leg -> AddEntry(graph_fit2mcv1,  "invmass fit - MC",                 "PL");
      leg -> Draw("same");
      
      fileOut -> cd("rapidity/canvas");
      cc.Write();
    }
  //******************************************************************
  
  //************ pT **************************************************
  fileOut->mkdir("pT/graphs");
  fileOut->mkdir("pT/canvas");
  
  for(int i_C=1; i_C<=C_nbins; i_C++)
    for(int i_y=1; i_y<=y_nbins; i_y++)
    {
      graph_fit2mcfit = GetGraph(h_fit2mcfit, "z", i_C, i_y);
      graph_fit2mcv1 = GetGraph(h_fit2mcv1, "z", i_C, i_y);
      
      SetAxesNames(graph_fit2mcfit, "p_{T}, GeV/c", "\" #chi^{2} \"");
      SetAxesNames(graph_fit2mcv1, "p_{T}, GeV/c", "\" #chi^{2} \"");
      
      SetGraphProperties(graph_fit2mcfit,   kBlue, line_width, 1, marker_size, 8);
      SetGraphProperties(graph_fit2mcv1,    kRed,  line_width, 1, marker_size, 8);
      
      TString graphname = "C" + StringBinNumber(i_C) + "_y" + StringBinNumber(i_y);
      
      fileOut -> cd("pT/graphs");
      graph_fit2mcfit -> Write(graphname + "_fit_mcfit");
      graph_fit2mcv1 -> Write(graphname + "_fit_mcv1");
      
      TCanvas cc(graphname, graphname, 1500, 900);
      cc.cd();
      
      CustomizeGraphsRanges({graph_fit2mcfit, graph_fit2mcv1});
      graph_fit2mcfit -> Draw("APL");
      graph_fit2mcv1 -> Draw("PL");
      
      TLegend* leg = new TLegend(0.25, 0.75, 0.65, 0.88);
      leg -> SetBorderSize(0);
      leg -> AddEntry(graph_fit2mcfit, "invmass fit - RECO with MC-match", "PL");
      leg -> AddEntry(graph_fit2mcv1,  "invmass fit - MC",                 "PL");
      leg -> Draw("same");
      
      fileOut -> cd("pT/canvas");
      cc.Write();
    }
  //******************************************************************
  

  
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