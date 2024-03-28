void Phasediagram() {
  gROOT->Macro( "phasediagram_style.cc" );
  const float marker_size(1.8);
  const int line_width = 3;
  const int text_font = 62;
  const float text_size = 0.035;
  auto crossover_color = kRed;
  auto phase_transition_color = kBlue;

  const float critical_point_x = 475;
  const float phase_transition_cross_Ox_x = 1200;
  const float critical_temperature = 156.5;
  const float phase_transition_interrupt_y = 0;
  const float alpha = 5;

  const float end = phase_transition_cross_Ox_x; // T = a*pow((end - mu_b), 1/alpha);
  const float a = critical_temperature / std::pow(end, 1/alpha);
  const float phase_transition_from = critical_point_x;
  const float phase_transition_to = end - std::pow(phase_transition_interrupt_y/a, alpha);

  TH1F* axes = new TH1F("ax", "ax", 1000, 0, 1600);
  axes->GetXaxis()->SetTitle("Baryon chemical potential, MeV");
  axes->GetYaxis()->SetTitle("Temperature, MeV");
  axes->GetXaxis()->CenterTitle();
  axes->GetYaxis()->CenterTitle();
  axes->GetXaxis()->SetRangeUser(0, 1600);
  axes->GetYaxis()->SetRangeUser(0, 250);
  axes->Draw("axis");

  TF1* phase_transition = new TF1("f1", "[0]*pow(([1]-x), 1/[2])", phase_transition_from, phase_transition_to);
  phase_transition->SetParameters(a, end, alpha);
  phase_transition->SetLineColor(phase_transition_color);
  phase_transition->SetLineWidth(line_width);
  phase_transition->Draw("same");

  TF1* crossover = new TF1("f2", "[0]*pow(([1]-x), 1/[2])", 0, phase_transition_from);
  crossover->SetParameters(a, end, alpha);
  crossover->SetLineStyle(7);
  crossover->SetLineColor(crossover_color);
  crossover->SetLineWidth(line_width);
  crossover->Draw("same");

  TF1* csc_line = new TF1("f3", "[0]*sqrt(x-[1])", 1200, 1590);
  csc_line->SetParameters(4.5, 1000);
  csc_line->SetLineColor(kMagenta);
  csc_line->SetLineWidth(line_width);
  csc_line->Draw("same");

  TGraph* lQCD = new TGraph();
  lQCD->AddPoint(10, critical_temperature);
  lQCD->SetMarkerStyle(kFullCircle);
  lQCD->SetMarkerSize(marker_size);
  lQCD->SetMarkerColor(kRed);
  lQCD->Draw("same P");

  const float critical_point_y = phase_transition->Eval(critical_point_x);
  TGraph* CP = new TGraph();
  CP->AddPoint(critical_point_x, critical_point_y);
  CP->SetMarkerStyle(kFullCircle);
  CP->SetMarkerSize(marker_size);
  CP->SetMarkerColor(phase_transition_color);
  CP->Draw("same P");

  TGraph* nucleus = new TGraph();
  nucleus->AddPoint(923, 5);
  nucleus->SetMarkerStyle(kFullCircle);
  nucleus->SetMarkerSize(marker_size);
  nucleus->SetMarkerColor(kGreen+2);
  nucleus->Draw("same P");

  TText* crossover_text = new TText(120, critical_temperature+5, "Crossover");
  crossover_text->SetTextColor(crossover_color);
  crossover_text->SetTextAngle(-7);
  crossover_text->SetTextFont(text_font);
  crossover_text->SetTextSize(text_size);
  crossover_text->Draw("same");

  TText* phase_transition_text = new TText(650, critical_temperature-15, "1-st order phase transition?");
  phase_transition_text->SetTextColor(phase_transition_color);
  phase_transition_text->SetTextAngle(-15);
  phase_transition_text->SetTextFont(text_font);
  phase_transition_text->SetTextSize(text_size);
  phase_transition_text->Draw("same");

  TText* CP_text = new TText(critical_point_x-120, critical_point_y-20, "Critical point?");
  CP_text->SetTextColor(kBlack);
  CP_text->SetTextFont(text_font);
  CP_text->SetTextSize(text_size);
  CP_text->Draw("same");

  TText* nucleus_text = new TText(923-120, 12, "Nuclear matter");
  nucleus_text->SetTextColor(kGreen+2);
  nucleus_text->SetTextFont(text_font);
  nucleus_text->SetTextSize(text_size);
  nucleus_text->Draw("same");

  TText* hadron_gas_text = new TText(300, 70, "Hadron gas");
  hadron_gas_text->SetTextColor(kBlack);
  hadron_gas_text->SetTextFont(text_font);
  hadron_gas_text->SetTextSize(text_size*1.5);
  hadron_gas_text->Draw("same");

  TText* QGP_text = new TText(300, 220, "Quark-gluon plasma");
  QGP_text->SetTextColor(kBlack);
  QGP_text->SetTextFont(text_font);
  QGP_text->SetTextSize(text_size*1.5);
  QGP_text->Draw("same");

  TText* LHC_text = new TText(200, 180, "LHC");
  LHC_text->SetTextColor(kGray+2);
  LHC_text->SetTextFont(text_font);
  LHC_text->SetTextSize(text_size*1.5);
  LHC_text->SetTextAngle(80);
  LHC_text->Draw("same");

  TText* RHIC_text = new TText(500, 160, "RHIC");
  RHIC_text->SetTextColor(kGray+2);
  RHIC_text->SetTextFont(text_font);
  RHIC_text->SetTextSize(text_size*1.5);
  RHIC_text->SetTextAngle(75);
  RHIC_text->Draw("same");

  TText* FAIR_text = new TText(800, 140, "FAIR");
  FAIR_text->SetTextColor(kGray+2);
  FAIR_text->SetTextFont(text_font);
  FAIR_text->SetTextSize(text_size*1.5);
  FAIR_text->SetTextAngle(70);
  FAIR_text->Draw("same");

  TText* CSC_text_1 = new TText(1220, 50, "Color super-");
  CSC_text_1->SetTextColor(kBlack);
  CSC_text_1->SetTextFont(text_font);
  CSC_text_1->SetTextSize(text_size*1.5);
  CSC_text_1->Draw("same");

  TText* CSC_text_2 = new TText(1220, 30, "conductor?");
  CSC_text_2->SetTextColor(kBlack);
  CSC_text_2->SetTextFont(text_font);
  CSC_text_2->SetTextSize(text_size*1.5);
  CSC_text_2->Draw("same");

}
