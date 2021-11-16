void efficiency_calc()
{
  std::string sim_file_name = "/home/user/cbmdir/working/effmap/simmap_pt_y_phi.dcmqgsm.apr20.xi.root";
  std::string rec_file_name = "/home/user/cbmdir/working/effmap/recmap_pt_y.dcmqgsm.apr20.xi.defcuts.root";
  
  std::string sim_histo_name = "hsim_y_pt";
  std::string rec_histo_name = "hrec_y_pt";
  
  TFile* fileSim = TFile::Open(sim_file_name.c_str(), "read");
  TFile* fileRec = TFile::Open(rec_file_name.c_str(), "read");
  
  TH2F* histoSim = fileSim->Get<TH2F>(sim_histo_name.c_str());
  TH2F* histoRec = fileRec->Get<TH2F>(rec_histo_name.c_str());
  
  histoSim->Sumw2();
  histoRec->Sumw2();
  
  TH2F* histoEff = (TH2F*)histoSim->Clone();
  TH2F* histoErr = (TH2F*)histoSim->Clone();
  histoEff -> Reset();
  histoErr -> Reset();
  
  for(int i=0; i<=histoSim->GetNbinsX(); i++)
    for(int j=0; j<=histoSim->GetNbinsY(); j++)
    {
      const float Nsim = histoSim->GetBinContent(i, j);
      const float Nrec = histoRec->GetBinContent(i, j);
      if(Nsim == 0) continue;
      histoEff->SetBinContent(i, j, Nrec/Nsim);
      if(Nrec == 0) continue;
      histoErr->SetBinContent(i, j, TMath::Sqrt(1/Nsim + 1/Nrec));
    }
    
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  histoSim -> Write("sim");
  histoRec -> Write("rec");
  histoEff -> Write("eff");
  histoErr -> Write("err");
  fileOut -> Close();
}