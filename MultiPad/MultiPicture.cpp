//
// Created by oleksii on 05.09.23.
//

#include "MultiPicture.h"

#include <iostream>
#include <sstream>
#include <utility>

MultiPicture::MultiPicture(int nx, int ny) : nx_{nx},
                                             ny_{ny} {
  pad_names_.resize(nx_ * ny_);
  left_margins_.resize(nx_);
  right_margins_.resize(nx_);
  top_margins_.resize(ny_);
  bottom_margins_.resize(ny_);
}

void MultiPicture::Run() {
  for(int j=0; j<ny_; j++) {
    for(int i=0; i<nx_; i++) {
      CropPicture(i, j);
    }
    MergeLine(j);
  }
  MergeAll();
  if(!save_intermediate_pictures_) {
    ExeBash("rm cropped_* merged_*");
  }
  if(remove_original_pictures_) {
    for(auto& pn : pad_names_) {
      ExeBash(("rm " + pn).c_str());
    }
  }
}

std::pair<int, int> MultiPicture::DeterminePicturesWH(const std::string& name) {
  FILE* pFile;
  char line[1024];

  pFile = popen(("identify -format '%w ' " + name).c_str(), "r");
  fgets(line, sizeof(line), pFile);
  int W;
  sscanf(line, "%d", &W);

  pFile = popen(("identify -format '%h ' " + name).c_str(), "r");
  fgets(line, sizeof(line), pFile);
  int H;
  sscanf(line, "%d", &H);

  return std::make_pair(W, H);
}

void MultiPicture::CropPicture(std::string inname, float left, float right, float bottom, float top, const std::string& outname) const {
  bool need_remove_png{false};
  if(inname.substr(inname.size()-3) == "pdf") {
    Pdf2Png(inname);
    inname.erase(inname.size()-3);
    inname += "png";
    need_remove_png = true;
  }
  const std::pair<int, int> aspect_ratio = DeterminePicturesWH(inname);
  const float fNpixL = left * (float)aspect_ratio.first;
  const float fNpixB = bottom * (float)aspect_ratio.second;
  const float fNpixR = right * (float)aspect_ratio.first;
  const float fNpixT = top * (float)aspect_ratio.second;
  const std::string sNpixL = to_string_with_precision(fNpixL, 0);
  const std::string sNpixB = to_string_with_precision(fNpixB, 0);
  const std::string sNpixR = to_string_with_precision(fNpixR, 0);
  const std::string sNpixT = to_string_with_precision(fNpixT, 0);

  std::string command = (std::string)"convert " + inname +
                              " -gravity north -chop 0x" + sNpixT +
                              " -gravity east -chop " + sNpixR + "x0" +
                              " -gravity south -chop 0x" + sNpixB +
                              " -gravity west -chop " + sNpixL + "x0" +
                              " " + outname;

  ExeBash(command);
  if(need_remove_png) {
    command = "rm " + inname;
    ExeBash(command);
  }
}

void MultiPicture::CropPicture(int i, int j) const {
  const std::string outname = "cropped_" + std::to_string(i) + "_" + std::to_string(j) + ".png";
  const float l = i == 0 ? 0 : left_margins_.at(i);
  const float r = i == nx_-1 ? 0 : right_margins_.at(i);
  const float t = j == 0 ? 0 : top_margins_.at(j);
  const float b = j == ny_-1 ? 0 : bottom_margins_.at(j);
  CropPicture(pad_names_.at(TransformCoordinates(i, j)), l, r, b, t, outname);
}

void MultiPicture::Pdf2Png(const std::string& inname) const {
  if(inname.substr(inname.size()-3) != "pdf") {
    throw std::runtime_error("MultiPicture::Pdf2Png() - requires pdf file as input");
  }
  const std::string command = "pdftoppm -png -cropbox -singlefile " + inname + " " + inname.substr(0, inname.size()-4);

  ExeBash(command);
}

void MultiPicture::MergeLine(int j) const {
  const std::string command = "convert +append cropped_*_" + std::to_string(j) + ".png merged_" + std::to_string(j) + ".png";

  ExeBash(command);
}

void MultiPicture::MergeAll() const {
  std::string command = "convert -append ";
  for(int j=0; j<ny_; j++) {
    command += "merged_" + std::to_string(j) + ".png ";
  }
  command += "out.png";

  ExeBash(command);
}

int MultiPicture::TransformCoordinates(int i, int j) const {
  if(i>=nx_ || j>=ny_ || i<0 || j<0) {
    std::cout << "Error: MultiPicture::TransformCoordinates() - i>=nx_ || j>=ny_ || i<0 || j<0\n";
    std::cout << "i = " << i << ", j = " << j << ", nx_ = " << nx_ << ", ny_ = " << ny_ << "\n";
    throw;
  }
  const int k = nx_*j + i;
  return k;
}

void MultiPicture::SetPictureName(int i, int j, std::string name) {
  pad_names_.at(TransformCoordinates(i, j)) = std::move(name);
}

void MultiPicture::SetLeftMargins(const std::vector<float>& margins) {
  if(margins.size() != nx_) {
    throw std::runtime_error("MultiPicture::SetLeftMargins() - margins.size() != nx_");
  }
  left_margins_ = margins;
}

void MultiPicture::SetLeftMargins(float margins) {
  std::vector<float> vec = std::vector<float>(nx_, margins);
  SetLeftMargins(vec);
}

void MultiPicture::SetRightMargins(const std::vector<float>& margins) {
  if(margins.size() != nx_) {
    throw std::runtime_error("MultiPicture::SetRightMargins() - margins.size() != nx_");
  }
  right_margins_ = margins;
}

void MultiPicture::SetRightMargins(float margins) {
  std::vector<float> vec = std::vector<float>(nx_, margins);
  SetRightMargins(vec);
}

void MultiPicture::SetTopMargins(const std::vector<float>& margins) {
  if(margins.size() != ny_) {
    throw std::runtime_error("MultiPicture::SetTopMargins() - margins.size() != nx_");
  }
  top_margins_ = margins;
}

void MultiPicture::SetTopMargins(float margins) {
  std::vector<float> vec = std::vector<float>(ny_, margins);
  SetTopMargins(vec);
}

void MultiPicture::SetBottomMargins(const std::vector<float>& margins) {
  if(margins.size() != ny_) {
    throw std::runtime_error("MultiPicture::SetBottomMargins() - margins.size() != nx_");
  }
  bottom_margins_ = margins;
}

void MultiPicture::ZeroAllMargins() {
  SetBottomMargins(0.f);
  SetTopMargins(0.f);
  SetLeftMargins(0.f);
  SetRightMargins(0.f);
}

void MultiPicture::SetBottomMargins(float margins) {
  std::vector<float> vec = std::vector<float>(ny_, margins);
  SetBottomMargins(vec);
}

void MultiPicture::ExeBash(const std::string& command) const {
  if(verbose_) std::cout << command << "\n";
  std::system(command.c_str());
}

template<typename T>
std::string MultiPicture::to_string_with_precision(const T a_value, int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
