#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TLegend.h>
#include <TPaveText.h>
// #include "plot_tools.C"

#include "AliHFInvMassFitter.h"
// #include "AliHFMassFitter.h"
#endif

Double_t double_gaus(Double_t *x, Double_t *par)
{
   Double_t f1 = par[0]*((1.-par[3])/TMath::Sqrt(2.*TMath::Pi())/par[2]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]));
   Double_t f2 = par[0]*(par[3]/TMath::Sqrt(2.*TMath::Pi())/par[4]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[4]/par[4]));
   Double_t f3=(f1+f2);
   return f3;
}

std::string StringBinNumber(int number)
{
  if(number<10)
    return "0" + std::to_string(number);
  else
    return std::to_string(number);
}

void SetAxesNames(TH3F* histo, TString xaxisname="centrality, %",
                  TString yaxisname="rapidity",
                  TString zaxisname="p_{T}, GeV")
{
  histo -> GetXaxis() -> SetTitle(xaxisname);
  histo -> GetYaxis() -> SetTitle(yaxisname);
  histo -> GetZaxis() -> SetTitle(zaxisname);
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int main()
{
//   TString sgnlfilename="/home/user/cbmdir/working/qna/OLD/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.sgnl_12.root";
  TString sgnlfilename="/home/user/cbmdir/working/qna/OLD/shapes/out.mass3D.apr20.dcmqgsm.nopid.lightcuts1.set4.sgnl_12.reduced.root";
  TFile* sgnlfile = TFile::Open(sgnlfilename, "read");
  
  TH1F* histosgnl = nullptr;
  
  const int C_nbins = 3;                                                    // TODO rm hardcoded nbins
  const int y_nbins = 4;
  const int pT_nbins = 4;
  
  double C_edges_array[] = {0, 20, 40, 100};                                // TODO rm hardcoded binranges
  double y_edges_array[] = {1.02179, 1.42179, 1.82179, 2.22179, 2.62179};
  double pT_edges_array[] = {0.2, 0.5, 0.8, 1.1, 1.4};
  
  double* C_edges = &C_edges_array[0];
  double* y_edges = &y_edges_array[0];
  double* pT_edges = &pT_edges_array[0];
  
//   TH3F hchi2_root("hchi2_root", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_root);
//   
//   TH3F hchi2_alice("hchi2_alice", "", C_nbins, C_edges, y_nbins, y_edges, pT_nbins, pT_edges);
//   SetAxesNames(&hchi2_alice);  
  
  bool is_first_canvas = true;
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  fileOut->mkdir("fit");
  fileOut->mkdir("ratio");

  for(int iC=0; iC<C_nbins; iC++)
    for(int iy=0; iy<y_nbins; iy++)
      for(int ipT=0; ipT<pT_nbins; ipT++)
      {
        std::string binname = "C" + StringBinNumber(iC+1) + "_y" + StringBinNumber(iy+1) + "_pT" + StringBinNumber(ipT+1);
                
//         histosgnl = (TH1F*)((TH1F*)sgnlfile->Get(binname.c_str()));
        histosgnl = (TH1F*)((TH1F*)sgnlfile->Get(binname.c_str())) -> Rebin(5);
        
        const double intest = (double)histosgnl->GetEntries()*histosgnl->GetBinWidth(1);
                
        Double_t minMassForFit0 = histosgnl->GetMean()-4*(histosgnl->GetRMS());
        Double_t maxMassForFit0 = histosgnl->GetMean()+4*(histosgnl->GetRMS());        
        
        TF1 *f2 = new TF1("double_gaus",double_gaus,minMassForFit0,maxMassForFit0,5);
        f2->SetParameters(intest,1.1156,0.001,0.2,0.0009);
        
        std::cout << "\n\n\t\t------------------------- Started ROOT fitting bin\t" << binname << "\t\t------------------------\n";
//         histosgnl->Fit(f2,"MNRI","L",minMassForFit0,maxMassForFit0);     
        histosgnl->Fit(f2,"R,L,E,+,0", "",minMassForFit0,maxMassForFit0);     
        std::cout << "\n\n\t\t------------------------- Finished ROOT fitting bin\t" << binname << "\t\t------------------------\n";
        
        AliHFInvMassFitter* fitter0 = new AliHFInvMassFitter(histosgnl, minMassForFit0, maxMassForFit0, 3, 1);
        fitter0->SetInitialGaussianMean(1.1156);
        fitter0->SetInitialGaussianSigma(0.001);
        fitter0->SetInitialSecondGaussianSigma(0.0009);
        
        std::cout << "\n\n\t\t------------------------- Started ALICE fitting bin\t" << binname << "\t\t------------------------\n";
        Bool_t ok = fitter0->MassFitter(kFALSE);
        std::cout << "\n\n\t\t------------------------- Finished ALICE fitting bin\t" << binname << "\t\t------------------------\n";
        
        if(!ok)
          std::cout << binname << " is not ok\n";
        
        TF1 *sig_func = fitter0->GetSignalFunc();
//         TF1 *tot_func = fitter0->GetMassFunc();
        
        f2->SetLineColor(kRed);
        sig_func->SetLineColor(kGreen+2);
        sig_func->SetLineStyle(2);
//         tot_func->SetLineColor(kMagenta);
        
        TLegend leg(0.15, 0.75, 0.3, 0.85);
        leg.AddEntry(f2, "ROOT", "L");
        leg.AddEntry(sig_func, "ALICE", "L");
        
        const float chi2root = f2->GetChisquare()/f2->GetNDF();
        const float chi2alice = fitter0->GetReducedChiSquare();
        
        TPaveText pt(0.15, 0.55, 0.3, 0.65, "brNDC");
        pt.AddText(("#chi^{2}/_{ndf} root = " + to_string_with_precision(chi2root, 2)).c_str());
        pt.AddText(("#chi^{2}/_{ndf} alice = " + to_string_with_precision(chi2alice, 2)).c_str());
        pt.SetFillColor(0);
        pt.SetTextSize(0.03);
        pt.SetTextFont(22);
        
        TCanvas cc("", "", 1500, 900);
        cc.cd();
        histosgnl->GetXaxis()->SetRangeUser(minMassForFit0-histosgnl->GetRMS(), maxMassForFit0+histosgnl->GetRMS());
        histosgnl->GetYaxis()->SetRangeUser(0, 1.2*std::max(histosgnl->GetYaxis()->GetXmax(), std::max(f2->GetMaximum(), sig_func->GetMaximum())));
        histosgnl->SetMarkerColor(kBlack);
        histosgnl->SetMarkerSize(0.8);
        histosgnl->SetMarkerStyle(8);
        histosgnl->Draw("PE,X0");
        f2->Draw("same");
        sig_func->Draw("same");
//         tot_func->Draw("same");
        leg.Draw("same");
        pt.Draw("same");
        fileOut->cd("fit");
        cc.Write(binname.c_str());
        
        TCanvas cratio("", "", 1500, 900);
        cratio.cd();
        histosgnl->Sumw2();
        histosgnl->Divide(f2);
        histosgnl->GetYaxis()->SetRangeUser(0.5,2);
//         cratio.SetLogy();
        histosgnl->Draw("PE,X0");
        TF1 one("one", "1", minMassForFit0-histosgnl->GetRMS(), maxMassForFit0+histosgnl->GetRMS());
        one.Draw("same");
        fileOut->cd("ratio");
        cratio.Write(binname.c_str());
        
        if(is_first_canvas) {
          cc.Print("fit.pdf(", "pdf");
          cratio.Print("ratio.pdf(", "pdf");
        }
        else {
          cc.Print("fit.pdf", "pdf");
          cratio.Print("ratio.pdf", "pdf");
        }
        
        is_first_canvas = false;
        
//         hchi2_root.SetBinContent(iC+1, iy+1, ipT+1, chi2root);
//         hchi2_alice.SetBinContent(iC+1, iy+1, ipT+1, chi2alice);
      }
      
  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("fit.pdf)", "pdf");
  emptycanvas.Print("ratio.pdf)", "pdf");
      
//   hchi2_root.Write();
//   hchi2_alice.Write();
  fileOut -> Close();
  
  return 0;
}