TGraphErrors* GetGraph(TH3F* histo3D, TString axis, int bin1, int bin2);

void TH3ToTGraph()
{
  TFile* file_fit = TFile::Open("/home/user/cbmdir/QnAnalysis/build/src/QnAnalysisDiscriminator/out.fitter.apr20.dcmqgsm.nopid.defcuts.set3.root");
//   TFile* file_fit = TFile::Open("/home/user/cbmdir/QnAnalysis/build/src/QnAnalysisDiscriminator/out.mcfitter.apr20.dcmqgsm.nopid.defcuts.set3.root");
//   TFile* file_fit = TFile::Open("/home/user/cbmdir/QnAnalysis/build/src/QnAnalysisDiscriminator/out.mcv1.apr20.dcmqgsm.nopid.defcuts.set3.root");
  
  TH3F* hfit = file_fit -> Get<TH3F>("parameters/hsignal");
//   TH3F* hfit = file_fit -> Get<TH3F>("hv1_mc");
  
  TGraphErrors* g = GetGraph(hfit, "z", 1, 4);
  
  TFile* fileOut = TFile::Open("out.th3totgraph.root", "recreate");
  g -> SetMarkerStyle(8);
  g -> SetMarkerSize(1.6);
  g -> Write("graph");
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