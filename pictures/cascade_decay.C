TGraph* CircleFromCenter(float a, float b, float r, float phi1, float phi2);
TGraph* CircleFromArc(float x, float y, float r, float alpha, float L);
TGraph* Line(float x1, float y1, float x2, float y2);

void cascade_decay() {
  TMultiGraph* mgr = new TMultiGraph();

  const float lambda_decay_x = 10;
  const float lambda_decay_y = 5;

  const float proton_r = 16;  float proton_alpha = 0.08; float proton_L = 10;
  const float pion_r = 12;    const float pion_alpha = -0.17;  const float pion_L = 10;

  const float lambda_alpha = (10*proton_alpha + pion_alpha) / 11;
  const float lambda_L = 5;

  const float xi_decay_x = lambda_decay_x - lambda_L*TMath::Cos(lambda_alpha);
  const float xi_decay_y = lambda_decay_y - lambda_L*TMath::Sin(lambda_alpha);

  const float xi_r = proton_r+3; const float xi_alpha = (10*lambda_alpha + pion_alpha)/11; float xi_L = 8;

  proton_alpha += TMath::Pi();  // because proton has charge opposite to pion
  proton_L = -proton_L;
  xi_L = -xi_L;

  //#################### MC ############################
  auto proton_sim = CircleFromArc(lambda_decay_x, lambda_decay_y, proton_r, proton_alpha, proton_L);
  proton_sim->SetLineColor(kBlue);
  proton_sim->SetLineWidth(3);
  proton_sim->SetLineStyle(1);
  mgr->Add(proton_sim, "L");

  auto pion_sim = CircleFromArc(lambda_decay_x, lambda_decay_y, pion_r, pion_alpha, pion_L);
  pion_sim->SetLineColor(kRed);
  pion_sim->SetLineWidth(3);
  pion_sim->SetLineStyle(1);
  mgr->Add(pion_sim, "L");

  auto lambda_sim = Line(lambda_decay_x, lambda_decay_y, xi_decay_x, xi_decay_y);
  lambda_sim->SetLineColor(kOrange-7);
  lambda_sim->SetLineWidth(3);
  lambda_sim->SetLineStyle(1);
  mgr->Add(lambda_sim, "L");

  auto xi_sim = CircleFromArc(xi_decay_x, xi_decay_y, xi_r, xi_alpha, xi_L);
  xi_sim->SetLineColor(kGreen+2);
  xi_sim->SetLineWidth(3);
  xi_sim->SetLineStyle(1);
  mgr->Add(xi_sim, "L");

  auto pion_fromxi_sim = CircleFromArc(xi_decay_x, xi_decay_y, pion_r, pion_alpha, pion_L);
  pion_fromxi_sim->SetLineColor(kViolet+1);
  pion_fromxi_sim->SetLineWidth(3);
  pion_fromxi_sim->SetLineStyle(1);
  mgr->Add(pion_fromxi_sim, "L");
  //######################################################

  const float xi_start_x = xi_sim->GetPointX(xi_sim->GetN()-1);
  const float xi_start_y = xi_sim->GetPointY(xi_sim->GetN()-1);

  //################### REC ##############################
  const float shift = 0.25;

  auto proton_rec = CircleFromArc(lambda_decay_x, lambda_decay_y+shift, proton_r, proton_alpha, proton_L);
  proton_rec->SetLineColor(kBlue);
  proton_rec->SetLineWidth(3);
  proton_rec->SetLineStyle(7);
  mgr->Add(proton_rec, "L");

  auto pion_rec = CircleFromArc(lambda_decay_x, lambda_decay_y-shift, pion_r, pion_alpha, pion_L);
  pion_rec->SetLineColor(kRed);
  pion_rec->SetLineWidth(3);
  pion_rec->SetLineStyle(7);
  mgr->Add(pion_rec, "L");

  auto lambda_rec = Line(lambda_decay_x, lambda_decay_y+shift/2, xi_decay_x, xi_decay_y+shift/2);
  lambda_rec->SetLineColor(kOrange-7);
  lambda_rec->SetLineWidth(3);
  lambda_rec->SetLineStyle(7);
  mgr->Add(lambda_rec, "L");

  auto xi_rec = CircleFromArc(xi_decay_x, xi_decay_y-shift/2, xi_r, xi_alpha, xi_L);
  xi_rec->SetLineColor(kGreen+2);
  xi_rec->SetLineWidth(3);
  xi_rec->SetLineStyle(7);
  mgr->Add(xi_rec, "L");

  auto pion_fromxi_rec = CircleFromArc(xi_decay_x, xi_decay_y-shift, pion_r, pion_alpha, pion_L);
  pion_fromxi_rec->SetLineColor(kViolet+1);
  pion_fromxi_rec->SetLineWidth(3);
  pion_fromxi_rec->SetLineStyle(7);
  mgr->Add(pion_fromxi_rec, "L");
  //######################################################

  TGraph* verteces = new TGraph();
  verteces->AddPoint(lambda_decay_x, lambda_decay_y);
  verteces->AddPoint(xi_decay_x, xi_decay_y);
  verteces->AddPoint(xi_start_x, xi_start_y);
  verteces->SetMarkerStyle(kFullCircle);
  verteces->SetMarkerSize(1);
  mgr->Add(verteces, "P");

  mgr->Draw("A");
}

TGraph* Line(float x1, float y1, float x2, float y2) {
  const int N = 1000;
  TGraph* gr = new TGraph();
  for (int i=0; i<N; i++) {
    const float x = x1 + (x2-x1)*i/N;
    const float y = y1 + (y2-y1)*i/N;
    gr->AddPoint(x,y);
  }

  return gr;
}

TGraph* CircleFromCenter(float a, float b, float r, float phi1, float phi2) {
  const int N = 1000;
  TGraph* gr = new TGraph();
  for (int i=0; i<N; i++) {
    const float phi = phi1 + (phi2-phi1)*i/N;
    const float x = a + r*TMath::Cos(phi);
    const float y = b + r*TMath::Sin(phi);
    gr->AddPoint(x,y);
  }

  return gr;
}

TGraph* CircleFromArc(float x, float y, float r, float alpha, float L) {
  const float a = x + r*TMath::Sin(alpha);
  const float b = y - r*TMath::Cos(alpha);
  const float phi1 = TMath::Pi()/2 + alpha;
  const float phi2 = phi1 - L/r;
  return CircleFromCenter(a, b, r, phi1, phi2);
}
