template <typename T>
void SetAxesNames(T* draw_obj, TString xaxisname, TString yaxisname, TString zaxisname = "") {
  draw_obj -> GetXaxis() -> SetTitle(xaxisname);
  draw_obj -> GetYaxis() -> SetTitle(yaxisname);
  if(zaxisname != "")
    draw_obj -> GetZaxis() -> SetTitle(zaxisname);
}

float MidRapidityByPbeam(float pbeam) {
  const float m = 0.938;
  const float E = std::sqrt(m*m + pbeam*pbeam);
  TLorentzVector v(0, 0, pbeam, E);
  
  return v.Rapidity()/2;
}