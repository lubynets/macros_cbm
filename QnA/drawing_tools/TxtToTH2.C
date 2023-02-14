void TxtToTH2(std::string fileInName) {
  std::ifstream fileIn;
  fileIn.open(fileInName);

  std::array<std::array<float, 9>, 9> values;
  std::array<std::array<float, 9>, 9> errors;

  TH2F histovalues("histovalues", "", 9, 0, 9, 9, 0, 9);
  TH2F histoerrors("histoerrors", "", 9, 0, 9, 9, 0, 9);

  std::string line;

  int I;
  int j;

  while(!fileIn.eof()) {
    std::getline(fileIn, line);
    std::istringstream iss(line);
    std::vector<std::string> words((std::istream_iterator<std::string>(iss)),
                                    std::istream_iterator<std::string>());
    if(words.size() != 10) continue;
    if(words.at(1) != "slope") continue;
    if(words.at(3) == "Y") continue;  // TODO how to deal with components

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
      histovalues.SetBinContent(ih+1, 9-jh, values.at(ih).at(jh));
      histoerrors.SetBinContent(ih+1, 9-jh, errors.at(ih).at(jh));
    }
  }

  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  histovalues.Write();
  histoerrors.Write();
  fileOut->Close();
}
