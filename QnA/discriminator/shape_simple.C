#include "Helper.hpp"

void shape_simple() {
  gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style_1.cc" );

  std::string evegen = "dcmqgsm"; std::string pbeam = "3.3";

  std::string fileName = "/home/oleksii/cbmdir/working/qna/aXmass/shapes/massDC." + evegen + "." + pbeam + "agev.oc1.310.root";
  const float mu = 0.497611;
  const float shift = 0.0476;
  const float tail_width = 0.02;
  const float mid_width = 0.027;
  std::string particle = "K^{0}_{S}";

  TFile* fileIn = TFile::Open(fileName.c_str(), "read");
  Qn::DataContainer<TH1F, Qn::Axis<double>> dcprimary = *(fileIn->Get<Qn::DataContainer<TH1F, Qn::Axis<double>>>("dcmass_all"));
//   auto dcIn = dcprimary.Rebin({"centrality", {0, 10, 20, 40, 70}});
  auto dcIn = dcprimary;

  Qn::DataContainerStatDiscriminator dcPurity;
  dcPurity.AddAxes(dcIn.GetAxes());

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");

  for (int i = 0; i < dcIn.size(); i++) {

    std::vector<size_t> indices = dcIn.GetIndex(i);
    std::string binname = "C" + StringBinNumber(indices.at(0) + 1) + "_pT" + StringBinNumber(indices.at(1) + 1) + "_y" + StringBinNumber(indices.at(2) + 1);
    const float C_lo = dcIn.GetAxis("centrality").GetLowerBinEdge(indices.at(0));
    const float C_hi = dcIn.GetAxis("centrality").GetUpperBinEdge(indices.at(0));
    const float pT_lo = dcIn.GetAxis("pT").GetLowerBinEdge(indices.at(1));
    const float pT_hi = dcIn.GetAxis("pT").GetUpperBinEdge(indices.at(1));
    const float y_lo = dcIn.GetAxis("y").GetLowerBinEdge(indices.at(2));
    const float y_hi = dcIn.GetAxis("y").GetUpperBinEdge(indices.at(2));

    std::cout << "DataContainer cell # " << i << "\n";
    std::cout << binname << "\n";
    std::cout << ("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str() << "\n";
    std::cout << ("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str() << "\n";
    std::cout << ("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str() << "\n";

    Sig2BckgrSimple sbs(&dcIn[i]);
    sbs.SetMu(mu);
    sbs.SetLeftShift(shift);
    sbs.SetRightShift(shift);
    sbs.SetLeftWidth(tail_width);
    sbs.SetRightWidth(tail_width);
    sbs.SetMidWidth(mid_width);

    sbs.Calculate();
    dcPurity[i].SetValue(sbs.GetPurity());

    CD(fileOut, "Histos");
    TCanvas cc("", "", 1500, 900);
    cc.cd();

    dcIn[i].SetTitle(binname.c_str());
    dcIn[i].GetXaxis()->SetRangeUser(mu - 2*shift, mu + 2*shift);
    dcIn[i].SetLineWidth(2);
    dcIn[i].SetStats(kFALSE);
    dcIn[i].Draw("");

    TPaveText binedges(0.12, 0.75, 0.27, 0.92, "brNDC");
    binedges.AddText(particle.c_str());
    binedges.AddText(("C: " + to_string_with_precision(C_lo, 2) + " - " + to_string_with_precision(C_hi, 2) + " %").c_str());
    binedges.AddText(("p_{T}: " + to_string_with_precision(pT_lo, 2) + " - " + to_string_with_precision(pT_hi, 2) + " GeV/c").c_str());
    binedges.AddText(("y_{LAB}: " + to_string_with_precision(y_lo, 2) + " - " + to_string_with_precision(y_hi, 2)).c_str());
    binedges.SetFillColor(0);
    binedges.SetTextSize(0.03);
    binedges.SetTextFont(22);
    binedges.Draw("same");

    TLine line1(mu-shift-tail_width/2, 0, mu-shift-tail_width/2, dcIn[i].GetMaximum()); line1.Draw("same");
    TLine line2(mu-shift+tail_width/2, 0, mu-shift+tail_width/2, dcIn[i].GetMaximum()); line2.Draw("same");
    TLine line3(mu+shift-tail_width/2, 0, mu+shift-tail_width/2, dcIn[i].GetMaximum()); line3.Draw("same");
    TLine line4(mu+shift+tail_width/2, 0, mu+shift+tail_width/2, dcIn[i].GetMaximum()); line4.Draw("same");
    TLine line5(mu-mid_width/2, 0, mu-mid_width/2, dcIn[i].GetMaximum()); line5.Draw("same");
    TLine line6(mu+mid_width/2, 0, mu+mid_width/2, dcIn[i].GetMaximum()); line6.Draw("same");

    TPaveText ptyield(0.12, 0.40, 0.27, 0.50, "brNDC");
    ptyield.AddText(("Purity = " + to_string_with_precision(sbs.GetPurity(), 3)).c_str());
    ptyield.SetFillColor(0);
    ptyield.SetTextSize(0.03);
    ptyield.SetTextFont(22);
    ptyield.Draw("same");

    cc.Write(binname.c_str());
  }

  fileOut->cd();
  dcPurity.Write("Purities");

  fileOut->Close();
  fileIn->Close();

  return 0;
}
