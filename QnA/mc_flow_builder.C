void mc_flow_builder()
{
  int num_init = 1;
  int num_fin = 5000;
  
  const float y_min = -0.5;
  const float y_max = 0.9;
  const int y_nbins = 7;
  
  const float flow_min = -1;
  const float flow_max = 1;
  const int flow_nbins = 2000;
  
  TH2D* histoflow = new TH2D("histoflow", "histoflow", y_nbins, y_min, y_max, flow_nbins, flow_min, flow_max);
  histoflow->GetXaxis()->SetTitle("y_{CM}");
  histoflow->GetYaxis()->SetTitle("cos(#varphi - #Psi_{RP})");
  
  TChain t("events");
  const TString dir = "/lustre/cbm/users/ogolosov/mc/generators/dcmqgsm_smm/auau/pbeam12agev/mbias/root/";
  for(int i=num_init; i<=num_fin; i++)
  {
    const TString num = std::to_string(i);
    const TString file = dir + "dcmqgsm_" + num + ".root";
    t.Add(file);
  }
  
  UEvent* event = new UEvent();
	t.SetBranchAddress("event", &event);
    
  const long nevents = t.GetEntries();
  
  for (long ievent=0; ievent<nevents; ievent++)
  {
    if ((ievent + 1) % 10000 == 0) 
        std::cout << ievent + 1 << " out of " << nevents << "\n";
    
    t.GetEntry(ievent);
    const int nparticles = event->GetNpa();
    const float psi_RP = event->GetPhi();
    
//     int multiplicity = 0;
    
    for (int iparticle=0; iparticle<nparticles; iparticle++)
    {
      const auto particle = event->GetParticle(iparticle);
      const int PDG = particle->GetPdg();
      if (PDG!=3122) continue;
      
      TLorentzVector momentum_cm = particle->GetMomentum();
      
      const float pt = momentum_cm.Pt();
      if(pt<0.3 || pt>1.35) continue;
      const float y_cm = momentum_cm.Rapidity();
      const float phi = momentum_cm.Phi();
      
      const float flow = TMath::Cos(phi-psi_RP);

      histoflow->Fill(y_cm, flow);
    }
  }
  
  TFile* f = TFile::Open("fileOut.root", "recreate");
  histoflow->Write();
  auto* profileflow = histoflow->ProfileX();
  profileflow->GetYaxis()->SetTitle("#LT cos(#varphi - #Psi_{RP}) #GT");
  profileflow->Write();
  f->Close();  
}
