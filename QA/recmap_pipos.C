#include "/lustre/cbm/users/lubynets/QA/macro/MacroHelper.h"

void recmap_pipos(std::string filelist_rec, float pbeam=12.) {

  std::cout << "Macro started\n";

  const TString y_name = "y_{CM}";
  const TString pT_name = "p_{T}, GeV/c";

  const int y_nbins = 1200;
  const int pT_nbins = 600;

  const float midrapidity = MidRapidityByPbeam(pbeam);
  const float y_low = -1;
  const float y_up = 2;
  const float pT_low = 0.;
  const float pT_up = 2.6;

  AnalysisTree::Chain* treeIn = new AnalysisTree::Chain(std::vector<std::string>({filelist_rec}), std::vector<std::string>({"aTree"}));
  auto* reco_tracks = new AnalysisTree::Particles();
  treeIn -> SetBranchAddress("RecParticles.", &reco_tracks);

  TH2F* hpty = new TH2F("hpty", "", y_nbins, y_low, y_up, pT_nbins, pT_low, pT_up);

  const int Nentries = treeIn->GetEntries();

  SetAxesNames(hpty, y_name, pT_name);

  for(int iEvent=0; iEvent<Nentries; iEvent++) {
    treeIn->GetEntry(iEvent);
    for(const auto& rectrack : *reco_tracks) {
      if(rectrack.GetPid() != 211) continue;

      const float pt = rectrack.GetPt();
      const float ycm = rectrack.GetRapidity() - midrapidity;

      hpty->Fill(ycm, pt);
    }
  }

  TFile fileOut("recmap_pipos.root", "recreate");
  hpty->Write();
  fileOut.Close();

  std::cout << "Macro finished successfully\n";
}
