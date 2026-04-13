//
// Created by oleksii on 05.09.23.
//

#include "MultiPicture.h"

#include <algorithm>
#include <iostream>
#include <utility>

MultiPicture::MultiPicture(int nx, int ny) : nx_{nx},
                                             ny_{ny} {
  pad_names_.resize(nx_ * ny_);
  left_margins_.resize(nx_, 0.f);
  right_margins_.resize(nx_, 0.f);
  top_margins_.resize(ny_, 0.f);
  bottom_margins_.resize(ny_, 0.f);
  pixel_sizes_x_.resize(nx_, 0);
  pixel_sizes_y_.resize(ny_, 0);
  left_margin_shifts_.resize(ny_-1, 0.f);
  bottom_margin_shifts_.resize(nx_-1, 0.f);
}

void MultiPicture::Run() {
  for(int j=0; j<ny_; j++) {
    for(int i=0; i<nx_; i++) {
      CropPicture(i, j);
    }
    MergeLine(j);
  }
  MergeAllWoMargins();
  if(!(left_margins_.at(0) == 0.f && bottom_margins_.at(ny_-1) == 0.f)) {
    for (int j = 0; j < ny_; j++) {
      if(left_margins_.at(0) != 0.f) CropLeftMargin(j);
    }
    for (int i = 0; i < nx_; i++) {
      if(bottom_margins_.at(ny_-1) != 0.f) CropBottomMargin(i);
    }
    if(left_margins_.at(0) != 0) {
      MergeLeftMargins();
    } else {
      ExeBash("mv out.womargins.png out.withleftmargins.png");
    }
    if(bottom_margins_.at(ny_-1) != 0.f) MergeBottomMargins();
    if(left_margins_.at(0) != 0) MergeLeftMarginsToAllWoMargins();
    MergeAll();
    ExeBash("mv out.png " + outfile_name_);
  } else {
    ExeBash("mv out.womargins.png " + outfile_name_);
  }
  if(!save_intermediate_pictures_) {
    ExeBash("rm cropped_* merged_* bottom_margin* left_margin* out.*margins.png");
  }
  if(remove_original_pictures_) {
    for(auto& pn : pad_names_) {
      ExeBash(("rm " + pn).c_str());
    }
  }
}

std::pair<int, int> MultiPicture::DeterminePicturesWH(const std::string& name) const {
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

  const int k = std::distance(pad_names_.begin(), std::find_if(pad_names_.begin(), pad_names_.end(), [&name] (const std::string& s) { return s.substr(0, s.size()-4) == name.substr(0, name.size()-4); }));
  auto [i, j] = TransformCoordinates(k);
  if(pixel_sizes_x_.at(i) != 0 && pixel_sizes_x_.at(i) != W) std::cout << "Warning: DeterminePicturesWH() - pixel_sizes_x_ != W\n";
  if(pixel_sizes_y_.at(j) != 0 && pixel_sizes_y_.at(j) != H) std::cout << "Warning: DeterminePicturesWH() - pixel_sizes_y_ != H\n";

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
  const std::string outname = "cropped_" + to_string_with_forward_zero(i) + "_" + to_string_with_forward_zero(j) + ".png";
  const float l = left_margins_.at(i);
  const float r = i == nx_-1 ? 0 : right_margins_.at(i);
  const float t = j == 0 ? 0 : top_margins_.at(j);
  const float b = bottom_margins_.at(j);
  CropPicture(pad_names_.at(TransformCoordinates(i, j)), l, r, b, t, outname);
}

void MultiPicture::CropLeftMargin(int j) const {
  const std::string outname = "left_margin_" + to_string_with_forward_zero(j) + ".png";
  if(j != 0 && top_margins_.at(j) + left_margin_shifts_.at(j-1) < 0) throw std::runtime_error("CropLeftMargin() - top_margins_.at(j) + left_margin_shifts_.at(j-1) < 0");
  if(j != ny_-1 && bottom_margins_.at(j) - left_margin_shifts_.at(j) < 0) throw std::runtime_error("CropLeftMargin() - bottom_margins_.at(j) - left_margin_shifts_.at(j) < 0");
  const float l = 0;
  const float r = 1 - left_margins_.at(0);
  const float t = j == 0 ? 0 : top_margins_.at(j) + left_margin_shifts_.at(j-1);
  const float b = j == ny_-1 ? bottom_margins_.at(j) : bottom_margins_.at(j) - left_margin_shifts_.at(j);
  CropPicture(pad_names_.at(TransformCoordinates(0, j)), l, r, b, t, outname);
}

void MultiPicture::CropBottomMargin(int i) const {
  const std::string outname = "bottom_margin_" + to_string_with_forward_zero(i) + ".png";
  if(i != 0 && left_margins_.at(i) - bottom_margin_shifts_.at(i-1) < 0) throw std::runtime_error("CropBottomMargin() - left_margins_.at(i) - bottom_margin_shifts_.at(i-1) < 0");
  if(i != nx_-1 && right_margins_.at(i) + bottom_margin_shifts_.at(i) < 0) throw std::runtime_error("CropBottomMargin() - right_margins_.at(i) + bottom_margin_shifts_.at(i) < 0");
  const float l = i == 0 ? 0 : left_margins_.at(i) - bottom_margin_shifts_.at(i-1);
  const float r = i == nx_-1 ? 0 : right_margins_.at(i) + bottom_margin_shifts_.at(i);
  const float t = 1 - bottom_margins_.at(ny_-1);
  const float b = 0;
  CropPicture(pad_names_.at(TransformCoordinates(i, ny_-1)), l, r, b, t, outname);
}

void MultiPicture::Pdf2Png(const std::string& inname) const {
  if(inname.substr(inname.size()-3) != "pdf") {
    throw std::runtime_error("MultiPicture::Pdf2Png() - requires pdf file as input");
  }
  const int k = std::distance(pad_names_.begin(), std::find(pad_names_.begin(), pad_names_.end(), inname.c_str()));
  auto [i, j] = TransformCoordinates(k);
  const std::string scaleTo = (pixel_sizes_x_.at(i) != 0 && pixel_sizes_y_.at(j) != 0) ? "-scale-to-x " + std::to_string(pixel_sizes_y_.at(j)) + " -scale-to-y " + std::to_string(pixel_sizes_x_.at(i)) + " " : "";
  const std::string command = "pdftoppm -png " + scaleTo + " -cropbox -singlefile " + inname + " " + inname.substr(0, inname.size()-4);

  ExeBash(command);
}

void MultiPicture::MergeLine(int j) const {
  if(merge_one_by_one_) {
    ExeBash("mv cropped_00_" + to_string_with_forward_zero(j) + ".png cropped_merge_base.png");
    for(int i=1; i<nx_; ++i) {
      const std::string command = "convert +append cropped_merge_base.png cropped_" + to_string_with_forward_zero(i) + "_" + to_string_with_forward_zero(j) + ".png cropped_merge_result.png";
      ExeBash(command);
      ExeBash("mv cropped_merge_result.png cropped_merge_base.png");
    }
    ExeBash("mv cropped_merge_base.png merged_" + to_string_with_forward_zero(j) + ".png");
  } else {
    const std::string command = "convert +append cropped_*_" + to_string_with_forward_zero(j) + ".png merged_" + to_string_with_forward_zero(j) + ".png";
    ExeBash(command);
  }
}

void MultiPicture::MergeAllWoMargins() const {
  if(merge_one_by_one_) {
    ExeBash("mv merged_00.png merged_merge_base.png");
    for(int j=1; j<ny_; j++) {
      ExeBash("convert -append merged_merge_base.png merged_" + to_string_with_forward_zero(j) + ".png merged_merge_result.png");
      ExeBash("mv merged_merge_result.png merged_merge_base.png");
    }
    ExeBash("mv merged_merge_base.png out.womargins.png");
  } else {
    std::string command = "convert -append ";
    for(int j=0; j<ny_; j++) {
      command += "merged_" + to_string_with_forward_zero(j) + ".png ";
    }
    command += "out.womargins.png";
    ExeBash(command);
  }
}

void MultiPicture::MergeLeftMargins() const {
  if(merge_one_by_one_) {
    ExeBash("mv left_margin_00.png left_margin_merge_base.png");
    for(int j=1; j<ny_; j++) {
      ExeBash("convert -append left_margin_merge_base.png left_margin_" + to_string_with_forward_zero(j) + ".png left_margin_merge_result.png");
      ExeBash("mv left_margin_merge_result.png left_margin_merge_base.png");
    }
    ExeBash("mv left_margin_merge_base.png out.leftmargins.png");
  } else {
    std::string command = "convert -append ";
    for(int j=0; j<ny_; j++) {
      command += "left_margin_" + to_string_with_forward_zero(j) + ".png ";
    }
    command += "out.leftmargins.png";
    ExeBash(command);
  }
}

void MultiPicture::MergeBottomMargins() const {
  if(merge_one_by_one_) {
    ExeBash("mv bottom_margin_00.png bottom_margin_merge_base.png");
    for(int i=1; i<nx_; i++) {
      ExeBash("convert +append bottom_margin_merge_base.png bottom_margin_" + to_string_with_forward_zero(i) + ".png bottom_margin_merge_result.png");
      ExeBash("mv bottom_margin_merge_result.png bottom_margin_merge_base.png");
    }
    ExeBash("mv bottom_margin_merge_base.png out.bottommargins.png");
  } else {
    std::string command = "convert +append ";
    for(int i=0; i<nx_; i++) {
      command += "bottom_margin_" + to_string_with_forward_zero(i) + ".png ";
    }
    command += "out.bottommargins.png";
    ExeBash(command);
  }
}

void MultiPicture::MergeLeftMarginsToAllWoMargins() const {
  const std::string command = "convert +append out.leftmargins.png out.womargins.png out.withleftmargins.png";
  ExeBash(command);
}

void MultiPicture::MergeAll() const {
  const std::string command = "convert -append out.withleftmargins.png out.bottommargins.png out.png";
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

std::pair<int, int> MultiPicture::TransformCoordinates(int k) const {
  if(k<0 || k>=nx_*ny_) {
    std::cout << "Error: MultiPicture::TransformCoordinates() - k<0 || k>=nx_*ny_\n";
    std::cout << "k = " << k << ", nx_ = " << nx_ << ", ny_ = " << ny_ << "\n";
    throw;
  }
  const int I = k%nx_;
  const int J = k/nx_;

  return std::make_pair(I, J);
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

void MultiPicture::SetBottomMargins(float margins) {
  std::vector<float> vec = std::vector<float>(ny_, margins);
  SetBottomMargins(vec);
}

void MultiPicture::SetPixelSizesX(const std::vector<int>& sizes) {
  if(sizes.size() != nx_) {
    throw std::runtime_error("MultiPicture::SetPixelSizesX() - sizes.size() != nx_");
  }
  pixel_sizes_x_ = sizes;
}

void MultiPicture::SetPixelSizesX(int sizes) {
  std::vector<int> vec = std::vector<int>(nx_, sizes);
  SetPixelSizesX(vec);
}

void MultiPicture::SetPixelSizesY(const std::vector<int>& sizes) {
  if(sizes.size() != ny_) {
    throw std::runtime_error("MultiPicture::SetPixelSizesY() - sizes.size() != ny_");
  }
  pixel_sizes_y_ = sizes;
}

void MultiPicture::SetPixelSizesY(int sizes) {
  std::vector<int> vec = std::vector<int>(ny_, sizes);
  SetPixelSizesY(vec);
}

void MultiPicture::SetLeftMarginShifts(const std::vector<float>& shifts) {
  if(shifts.size() != ny_-1) {
    throw std::runtime_error("MultiPicture::SetLeftMarginShifts() - shifts.size() != ny_-1");
  }
  left_margin_shifts_ = shifts;
}

void MultiPicture::SetBottomMarginShifts(const std::vector<float>& shifts) {
  if(shifts.size() != nx_-1) {
    throw std::runtime_error("MultiPicture::SetBottomMarginShifts() - shifts.size() != nx_-1");
  }
  bottom_margin_shifts_ = shifts;
}

void MultiPicture::SetLeftMarginShifts(float shifts) {
  std::vector<float> vec = std::vector<float>(ny_-1, shifts);
  SetLeftMarginShifts(vec);
}

void MultiPicture::SetBottomMarginShifts(float shifts) {
  std::vector<float> vec = std::vector<float>(nx_-1, shifts);
  SetBottomMarginShifts(vec);
}

void MultiPicture::ExeBash(const std::string& command) const {
  if(verbose_) std::cout << command << "\n";
  std::system(command.c_str());
}

std::string MultiPicture::to_string_with_forward_zero(int value) {
  return (value > 9 ? "" : "0") + std::to_string(value);
}

void MultiPicture::SetOutputFileName(std::string name) {
  outfile_name_ = std::move(name);
  if (outfile_name_.size() < 4 || outfile_name_.compare(outfile_name_.size() - 4, 4, ".png") != 0) outfile_name_ += ".png";
}
