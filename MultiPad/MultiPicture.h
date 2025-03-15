//
// Created by oleksii on 05.09.23.
//

#ifndef MULTIPAD_MULTIPICTURE_H
#define MULTIPAD_MULTIPICTURE_H

#include <utility>
#include <vector>
#include <string>

class MultiPicture {
public:
  MultiPicture() = delete;
  MultiPicture(int nx, int ny);
  ~MultiPicture() = default;

  void SetPictureName(int i, int j, std::string name);

  void SetVerbose(bool value = true) { verbose_ = value; }
  void SetSaveIntermediatePictures(bool value = true) { save_intermediate_pictures_ = value; }
  void SetRemoveOriginalPictures(bool value = true) { remove_original_pictures_ = value; }

  void SetLeftMargins(const std::vector<float>& margins);
  void SetLeftMargins(float margins);
  void SetRightMargins(const std::vector<float>& margins);
  void SetRightMargins(float margins);
  void SetTopMargins(const std::vector<float>& margins);
  void SetTopMargins(float margins);
  void SetBottomMargins(const std::vector<float>& margins);
  void SetBottomMargins(float margins);
  void SetLeftMarginShifts(const std::vector<float>& shifts);
  void SetLeftMarginShifts(float shifts);
  void SetBottomMarginShifts(const std::vector<float>& shifts);
  void SetBottomMarginShifts(float shifts);
  void ZeroAllMargins();

  void Run();

protected:
  const int nx_{};
  const int ny_{};
  std::vector<std::string> pad_names_{};
  std::vector<float> left_margins_{};
  std::vector<float> right_margins_{};
  std::vector<float> top_margins_{};
  std::vector<float> bottom_margins_{};
  std::vector<float> left_margin_shifts_{};
  std::vector<float> bottom_margin_shifts_{};
  bool verbose_{false};
  bool save_intermediate_pictures_{false};
  bool remove_original_pictures_{false};

  static std::pair<int, int> DeterminePicturesWH(const std::string& name) ;
  void CropPicture(std::string inname, float left, float right, float bottom, float top, const std::string& outname) const;
  void CropPicture(int i, int j) const;
  void CropLeftMargin(int j) const;
  void CropBottomMargin(int i) const;
  void Pdf2Png(const std::string& inname) const;
  void MergeLine(int j) const;
  void MergeAllWoMargins() const;
  void MergeLeftMargins() const;
  void MergeBottomMargins() const;
  void MergeLeftMarginsToAllWoMargins() const;
  void MergeAll() const;
  void ExeBash(const std::string& command) const;

  int TransformCoordinates(int i, int j) const;

  template <typename T>
  static std::string to_string_with_precision(T a_value, int n = 6);

};


#endif //MULTIPAD_MULTIPICTURE_H
