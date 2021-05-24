TCanvas* GraphToCanvas(TGraph* graph, std::string graphname);
std::string StringBinNumber(int number);

void TGraph_Fits_pdfer()
{
  TFile* file_fit = TFile::Open("/home/user/cbmdir/working/qna/fits/out.fitter.apr20.dcmqgsm.nopid.lightcuts1.set4.XX.root");
  const int C_nbins = 3;
  const int y_nbins = 4;
  const int pT_nbins = 4;
  
  bool is_first_canvas = true;
  
  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string graphname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
        TGraph* graph = file_fit -> Get<TGraph>(("fit/"+graphname).c_str());
        TCanvas* cc = GraphToCanvas(graph, graphname);
        
        if(is_first_canvas)
          cc->Print("out.fits.pdf(", "pdf");
        else
          cc->Print("out.fits.pdf", "pdf");

        is_first_canvas = false;
      }
      
  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("out.fits.pdf)", "pdf");  
}

TCanvas* GraphToCanvas(TGraph* graph, std::string graphname)
{
  TCanvas* cc = new TCanvas("canvas", "canvas", 1500, 900);
  graph -> SetTitle(graphname.c_str());
  graph -> Draw();
  
  return cc;
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}