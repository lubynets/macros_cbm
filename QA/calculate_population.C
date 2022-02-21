void calculate_population(std::string filelist)
{
  std::cout << "Macro started\n";
  
  int Lambda=0;
  int Xi=0;
  int Omega=0;
  int AntiLambda=0;
  int AntiXi=0;
  int AntiOmega=0;
  
  AnalysisTree::Chain* sim_tree = new AnalysisTree::Chain(std::vector<std::string>({filelist}), std::vector<std::string>({"rTree"}));
  
  auto* sim_tracks = new AnalysisTree::Particles();
  
  sim_tree -> SetBranchAddress("SimParticles", &sim_tracks);

//   const int sim_entries = sim_tree->GetEntries();
  const int sim_entries = 50000;
  
  for(int iEvent=0; iEvent<sim_entries; iEvent++)
  {
    sim_tree->GetEntry(iEvent);
    if(iEvent%100==0)
      std::cout << iEvent << std::endl;
        
    for(const auto& simtrack : *(sim_tracks->GetChannels()))
    {
      const int mother_id = simtrack.GetField<int>(sim_tree->GetConfiguration()->GetBranchConfig("SimParticles").GetFieldId("mother_id"));
//       if(mother_id != -1) continue;
      
      const int pid = simtrack.GetPid();
      if(pid == 3122)
        Lambda++;
      if(pid == 3312)
        Xi++;
      if(pid == 3334)
        Omega++;
      if(pid == -3122)
        AntiLambda++;
      if(pid == -3312)
        AntiXi++;
      if(pid == -3334)
        AntiOmega++;
    }    
  }
  std::cout << "\n";
  std::cout << "Lambda = " << Lambda << "\n";
  std::cout << "Xi = " << Xi << "\n";
  std::cout << "Omega = " << Omega << "\n";
  std::cout << "AntiLambda = " << AntiLambda << "\n";
  std::cout << "AntiXi = " << AntiXi << "\n";
  std::cout << "AntiOmega = " << AntiOmega << "\n";  
}