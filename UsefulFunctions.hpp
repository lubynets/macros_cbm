template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

std::pair<float, float> GetHistoYRanges(TH1* h) {
  float min = 1.e9;
  float max = -1.e9;
  for(int iBin=1; iBin<=h->GetNbinsX(); iBin++) {
    const float up = h->GetBinContent(iBin) + h->GetBinError(iBin);
    const float lo = h->GetBinContent(iBin) - h->GetBinError(iBin);
    max = std::max(max, up);
    min = std::min(min, lo);
  }
  const float diff = max - min;
  max = max + diff/10;
  min = min - diff/10;
  return std::make_pair(min, max);
}

template<typename T>
void RemoveFirstLableFromAxis(T* obj, const std::string& axisletter) {
  TAxis* axis{nullptr};
  if(axisletter == "x") {
    axis = obj->GetXaxis();
  } else if (axisletter == "y") {
    axis = obj->GetYaxis();
  } else throw;
  if(axis->IsVariableBinSize()) {
    TArrayD* array = new TArrayD(*axis->GetXbins());
    array->SetAt(array->At(0)+1e-3, 0);
    axis->Set(axis->GetNbins(), array->GetArray());
  } else {
    axis->Set(axis->GetNbins(), axis->GetXmin()+1e-3, axis->GetXmax());
  }
}
