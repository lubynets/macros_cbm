float DetermineCategorySlopeRatio(float value, float error);
float DetermineCategoryIntercept(float value, float error);
bool OverlapLineSegments(float A1, float A2, float B1, float B2);

void TxtToTH2(std::string fileInName) {
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

  std::array<std::array<float, 9>, 9> values;
  std::array<std::array<float, 9>, 9> errors;

  TH2F* histo = new TH2F("histo", "", 9, 0, 9, 9, 0, 9);

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
    if(words.at(3) == "X") continue;  // TODO how to deal with components

    if(words.at(0) == "etacut_1_charged") j = 0;
    if(words.at(0) == "etacut_1_all")     j = 1;
    if(words.at(0) == "psd1")             j = 2;
    if(words.at(0) == "etacut_2_charged") j = 3;
    if(words.at(0) == "etacut_2_all")     j = 4;
    if(words.at(0) == "psd2")             j = 5;
    if(words.at(0) == "etacut_3_charged") j = 6;
    if(words.at(0) == "etacut_3_all")     j = 7;
    if(words.at(0) == "psd3")             j = 8;

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

  for(int ih=0; ih<9; ih++) {
    for(int jh=0; jh<9; jh++) {
//       histo->SetBinContent(ih+1, 9-jh, values.at(ih).at(jh));
      if(is_slope) histo->SetBinContent(ih+1, 9-jh, DetermineCategorySlopeRatio(values.at(ih).at(jh), errors.at(ih).at(jh)));
      if(!is_slope) histo->SetBinContent(ih+1, 9-jh, DetermineCategoryIntercept(values.at(ih).at(jh), errors.at(ih).at(jh)));
    }
  }

  TCanvas* cc = new TCanvas("cc", "", 1500, 448);
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
  histo->Draw("col text");

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

//   if(std::fabs(value-1.) < n_sigma*error) {
//     if(value-n_sigma*error > lower && value+n_sigma*error < upper) return 2.;  // fAcceptable
//     else return 1.; // fCoincidence
//   }
//   if(value<0) return 5.;                      // fWrongSign
//   if(value > lower && value < upper) return 2.;  // fAcceptable
//   if(value >=upper) return 3.;                  // fOverEstimate
//   if(value <= lower) return 4.;                // fUnderEstimate

  if(value-n_sigma_eps*error > lower && value+n_sigma_eps*error < upper) return 2.; // fAcceptable, kGreen+2
  if(OverlapLineSegments(value - n_sigma_sgnfc*error, value + n_sigma_sgnfc*error, lower, upper)) return 1.; // fNotGoodNotBad kGray
  if(value<0) return 5.;                      // fWrongSign
  if(value - n_sigma_sgnfc*error >= upper) return 3.;                  // fOverEstimate kRed-...
  if(value + n_sigma_sgnfc*error <= lower) return 4.;                // fUnderEstimate kRed-...

  std::cout << value << "\t" << error << "\n";

  return 0.01;  // not legal
}

float DetermineCategoryIntercept(float value, float error) {
  if(std::fabs(value-0.) < 3*error)  return 1.; // fCoincidence
  if(value > 0.) return 3.;                  // fOverEstimate
  if(value < 0.) return 4.;                // fUnderEstimate
  return 0.;  // not legal
}

bool OverlapLineSegments(float A1, float A2, float B1, float B2) {
  if(B1>A1 && B1<A2) return true;
  if(B2>A1 && B2<A2) return true;
  if(A1>B1 && A1<B2) return true;
  if(A2>B1 && A2<B2) return true;

  return false;
}
