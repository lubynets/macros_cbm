TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2);
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
  
  TH3F* h_fit = file_fit -> Get<TH3F>("parameters/hsignal");
  TH3F* h_mcfit = file_mcfit -> Get<TH3F>("parameters/hsignal");
  TH3F* h_mcv1 = file_mcv1 -> Get<TH3F>("hv1_mc");
  
  const int C_nbins = h_fit->GetXaxis()->GetNbins();
  const int y_nbins = h_fit->GetYaxis()->GetNbins();
  const int pT_nbins = h_fit->GetZaxis()->GetNbins();
  
  TFile* fileOut = TFile::Open("out.th3totgrphcplx.root", "recreate");
    
  TGraphErrors* graph_fit;
  TGraphErrors* graph_mcfit;
  TGraphErrors* graph_mcv1;
  
  //************ centrality ******************************************
  fileOut->mkdir("centrality/graphs");
  fileOut->mkdir("centrality/canvas");
  
  for(int i_y=1; i_y<=y_nbins; i_y++)
    for(int i_pT=1; i_pT<=pT_nbins; i_pT++)
    {
      graph_fit = GetGraph(h_fit, "x", i_y, i_pT);
      graph_mcfit = GetGraph(h_mcfit, "x", i_y, i_pT);
      graph_mcv1 = GetGraph(h_mcv1, "x", i_y, i_pT);
      
      SetAxesNames(graph_fit, "Centrality, %", "v_{1x}");
      SetAxesNames(graph_mcfit, "Centrality, %", "v_{1x}");
      SetAxesNames(graph_mcv1, "Centrality, %", "v_{1x}");
      
      SetGraphProperties(graph_fit,   kBlue,     line_width, 1, marker_size, 8);
      SetGraphProperties(graph_mcfit, kGreen+2,  line_width, 1, marker_size, 8);
      SetGraphProperties(graph_mcv1,  kRed,      line_width, 1, marker_size, 8);
      
      TString graphname = "y" + StringBinNumber(i_y) + "_pT" + StringBinNumber(i_pT);
      
      fileOut -> cd("centrality/graphs");
      graph_fit -> Write(graphname + "_fit");
      graph_mcfit -> Write(graphname + "_mcfit");
      graph_mcv1 -> Write(graphname + "_mcv1");
      
      TCanvas cc(graphname, graphname, 1500, 900);
      cc.cd();
      
      graph_fit -> Draw("APL");
      graph_mcfit -> Draw("PL");
      graph_mcv1 -> Draw("PL");
      
      TLegend* leg = new TLegend(0.35, 0.68, 0.55, 0.88);
      leg -> SetBorderSize(0);
      leg -> AddEntry(graph_mcv1,  "MC",                "PL");
      leg -> AddEntry(graph_mcfit, "RECO, MC-match",    "PL");
      leg -> AddEntry(graph_fit,   "RECO, invmass-fit", "PL");
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
      graph_fit = GetGraph(h_fit, "y", i_C, i_pT);
      graph_mcfit = GetGraph(h_mcfit, "y", i_C, i_pT);
      graph_mcv1 = GetGraph(h_mcv1, "y", i_C, i_pT);
      
      SetAxesNames(graph_fit, "y_{LAB}", "v_{1x}");
      SetAxesNames(graph_mcfit, "y_{LAB}", "v_{1x}");
      SetAxesNames(graph_mcv1, "y_{LAB}", "v_{1x}");
      
      SetGraphProperties(graph_fit,   kBlue,     line_width, 1, marker_size, 8);
      SetGraphProperties(graph_mcfit, kGreen+2,  line_width, 1, marker_size, 8);
      SetGraphProperties(graph_mcv1,  kRed,      line_width, 1, marker_size, 8);
      
      TString graphname = "C" + StringBinNumber(i_C) + "_pT" + StringBinNumber(i_pT);
      
      fileOut -> cd("rapidity/graphs");
      graph_fit -> Write(graphname + "_fit");
      graph_mcfit -> Write(graphname + "_mcfit");
      graph_mcv1 -> Write(graphname + "_mcv1");
      
      TCanvas cc(graphname, graphname, 1500, 900);
      cc.cd();
      
      graph_fit -> Draw("APL");
      graph_mcfit -> Draw("PL");
      graph_mcv1 -> Draw("PL");
      
      TLegend* leg = new TLegend(0.35, 0.68, 0.55, 0.88);
      leg -> SetBorderSize(0);
      leg -> AddEntry(graph_mcv1,  "MC",                "PL");
      leg -> AddEntry(graph_mcfit, "RECO, MC-match",    "PL");
      leg -> AddEntry(graph_fit,   "RECO, invmass-fit", "PL");
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
      graph_fit = GetGraph(h_fit, "z", i_C, i_y);
      graph_mcfit = GetGraph(h_mcfit, "z", i_C, i_y);
      graph_mcv1 = GetGraph(h_mcv1, "z", i_C, i_y);
      
      SetAxesNames(graph_fit, "p_{T}, GeV/c", "v_{1x}");
      SetAxesNames(graph_mcfit, "p_{T}, GeV/c", "v_{1x}");
      SetAxesNames(graph_mcv1, "p_{T}, GeV/c", "v_{1x}");
      
      SetGraphProperties(graph_fit,   kBlue,     line_width, 1, marker_size, 8);
      SetGraphProperties(graph_mcfit, kGreen+2,  line_width, 1, marker_size, 8);
      SetGraphProperties(graph_mcv1,  kRed,      line_width, 1, marker_size, 8);
      
      TString graphname = "C" + StringBinNumber(i_C) + "_y" + StringBinNumber(i_y);
      
      fileOut -> cd("pT/graphs");
      graph_fit -> Write(graphname + "_fit");
      graph_mcfit -> Write(graphname + "_mcfit");
      graph_mcv1 -> Write(graphname + "_mcv1");
      
      TCanvas cc(graphname, graphname, 1500, 900);
      cc.cd();
      
      graph_fit -> Draw("APL");
      graph_mcfit -> Draw("PL");
      graph_mcv1 -> Draw("PL");
      
      TLegend* leg = new TLegend(0.35, 0.68, 0.55, 0.88);
      leg -> SetBorderSize(0);
      leg -> AddEntry(graph_mcv1,  "MC",                "PL");
      leg -> AddEntry(graph_mcfit, "RECO, MC-match",    "PL");
      leg -> AddEntry(graph_fit,   "RECO, invmass-fit", "PL");
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