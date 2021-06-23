void histo_extract(TString fileInName)
{
  TFile* fileIn = TFile::Open(fileInName, "read");
  TFile* fileOut = TFile::Open("CbmKFQA.root", "recreate");
  
  std::string path_to_histos = "KFTopoReconstructor/KFParticlesFinder/Particles/";
  std::vector<std::string> particles {"Lambda", "Xi-"};
  std::vector<std::string> components {"M", "p", "p_{t}", "y", "phi", "X", "Y", "Z"};
  
  for(auto& particle : particles)
    for(auto& component : components)
    {
      std::string histoname = path_to_histos + particle + "/Parameters/" + component;
      TH1F* histo = (TH1F*)fileIn->Get(histoname.c_str());
      
      if(fileOut->GetDirectory(particle.c_str()) == nullptr)
        fileOut -> mkdir(particle.c_str());
      fileOut -> cd(particle.c_str());

      histo->Write(component.c_str());
    }  
}