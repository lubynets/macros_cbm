float DetermineCategorySlopeRatio(float value, float error);
float DetermineCategoryIntercept(float value, float error);
bool OverlapLineSegments(float A1, float A2, float B1, float B2);

void TxtToTH2_mod(std::string fileInName) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow); // FIXME if there is no this line, for some reasons doesn't work custom palette

  bool is_slope{true};
//   bool is_slope{false};

  Int_t palette[6];
  palette[0] = kBlack;
  palette[1] = kGray;
  palette[2] = kGreen+2;
  palette[3] = kRed-7;
  palette[4] = kRed-9;
  palette[5] = kRed-4;
  gStyle->SetPalette(6, palette);

  std::ifstream fileIn;
  fileIn.open(fileInName);

  constexpr int size_x = 9;
  constexpr int size_y = 6;

  std::array<std::array<float, size_y>, size_x> values;
  std::array<std::array<float, size_y>, size_x> errors;

  TH2F* histo = new TH2F("histo", "", size_x, 0, size_x, size_y, 0, size_y);

  std::string line;

  int I;
  int j;

  while(!fileIn.eof()) {
    std::getline(fileIn, line);
    std::istringstream iss(line);
    std::vector<std::string> words((std::istream_iterator<std::string>(iss)),
                                    std::istream_iterator<std::string>());
    if(words.size() != 10) continue;
    if(is_slope && words.at(1) != "slope") continue;
    if(!is_slope && words.at(1) != "intercept") continue;
//     if(words.at(3) == "X") continue;  // TODO how to deal with components

    if(words.at(0) == "etacut_1_charged") j = 0;
    if(words.at(0) == "etacut_1_all")     j = 1;
    if(words.at(0) == "etacut_2_charged") j = 2;
    if(words.at(0) == "etacut_2_all")     j = 3;
    if(words.at(0) == "etacut_3_charged") j = 4;
    if(words.at(0) == "etacut_3_all")     j = 5;

    if(words.at(2) == "pt1") I = 0;
    if(words.at(2) == "pt2") I = 1;
    if(words.at(2) == "pt3") I = 2;

    values.at(0+I).at(j) = std::stof(words.at(4));
    values.at(3+I).at(j) = std::stof(words.at(6));
    values.at(6+I).at(j) = std::stof(words.at(8));
    errors.at(0+I).at(j) = std::stof(words.at(5));
    errors.at(3+I).at(j) = std::stof(words.at(7));
    errors.at(6+I).at(j) = std::stof(words.at(9));
  }
  fileIn.close();

  for(int ih=0; ih<size_x; ih++) {
    for(int jh=0; jh<size_y; jh++) {
      if(is_slope) histo->SetBinContent(ih+1, size_y-jh, DetermineCategorySlopeRatio(values.at(ih).at(jh), errors.at(ih).at(jh)));
      if(!is_slope) histo->SetBinContent(ih+1, size_y-jh, DetermineCategoryIntercept(values.at(ih).at(jh), errors.at(ih).at(jh)));
    }
  }

  TCanvas* cc = new TCanvas("cc", "", 1500, 500);
  cc->SetLeftMargin(0);
  cc->SetRightMargin(0);
  cc->SetTopMargin(0);
  cc->SetBottomMargin(0);
  histo->GetXaxis()->SetTickSize(0);
  histo->GetYaxis()->SetTickSize(0);
  histo->SetContour(6);
  histo->SetContourLevel(0, -0.5);
  histo->SetContourLevel(1, 0.5);
  histo->SetContourLevel(2, 1.5);
  histo->SetContourLevel(3, 2.5);
  histo->SetContourLevel(4, 3.5);
  histo->SetContourLevel(5, 4.5);
  histo->Draw("col");

  cc->SaveAs("cc.png");

//   TFile* fileOut = TFile::Open("fileOut.root", "recreate");
//   histo->Write();
//   cc->Write();
//   fileOut->Close();
}

float DetermineCategorySlopeRatio(float value, float error) {
  const float eps = 0.1;
  const float n_sigma_sgnfc = 3.;
  const float n_sigma_eps = 0.7;

  const float upper = 1+eps;
  const float lower = 1./(1+eps);

  if(value-n_sigma_eps*error > lower && value+n_sigma_eps*error < upper) return 2.;                          // fAcceptable, kGreen+2
  if(OverlapLineSegments(value - n_sigma_sgnfc*error, value + n_sigma_sgnfc*error, lower, upper)) return 1.; // fNotGoodNotBad kGray
  if(value<0) return 5.;                                                                                     // fWrongSign kRed-4
  if(value - n_sigma_sgnfc*error >= upper) return 3.;                                                        // fOverEstimate kRed-3
  if(value + n_sigma_sgnfc*error <= lower) return 4.;                                                        // fUnderEstimate kRed-9

  std::cout << value << "\t" << error << "\n";

  return 0.01;  // not legal
}

float DetermineCategoryIntercept(float value, float error) {
  const float n_sigma_sgnfc = 3.;
  if(std::fabs(value-0.) < n_sigma_sgnfc*error)  return 1.; // fCoincidence
  if(value > 0.) return 3.;                                 // fOverEstimate
  if(value < 0.) return 4.;                                 // fUnderEstimate
  return 0.01;  // not legal
}

bool OverlapLineSegments(float A1, float A2, float B1, float B2) {
  if(B1>A1 && B1<A2) return true;
  if(B2>A1 && B2<A2) return true;
  if(A1>B1 && A1<B2) return true;
  if(A2>B1 && A2<B2) return true;

  return false;
}
