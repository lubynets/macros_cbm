void CustomizeHistoRanges(std::vector<TH1F*> v_histo);

void pfsqa_kfpf_compare()
{
  gStyle -> SetOptStat(0);
  
  TFile* fileInKFPF = TFile::Open("/home/user/cbmdir/working/pfsqa/cascade/CbmKFQA.root");
  TFile* fileInPFS = TFile::Open("/home/user/cbmdir/working/pfsqa/cascade/pfsqa.casc.mcpid_alt.def2.xi.root");
  
  struct Particle
  {
    std::string kfpf_name_;
    std::string pfs_name_;
  };
  
  std::vector<Particle> particles
  {
//     {"Lambda", "Lambda"},
    {"Xi-", "Ksi"}
  };
  
  struct Component
  {
    std::string kfpf_name_;
    std::string pfs_name_;
  };
  
  std::vector<Component> components
  {
    {"M",     "mass"    },
    {"p",     "p"       },
    {"p_{t}", "pT"      },
    {"y",     "rapidity"},
    {"phi",   "phi"     },
    {"X",     "x"       },
    {"Y",     "y"       },
    {"Z",     "z"       }
  };
  
  std::string kfpf_path = "KFTopoReconstructor/KFParticlesFinder/Particles/";
  std::string pfs_path  = "Candidates_All_";
  
  TH1F* histo_KFPF;
  TH1F* histo_PFS;
  
  bool is_first_canvas = true;
    
  for(auto& pa : particles)
    for(auto& co : components)
    {
//       std::string histo_kfpf_name = kfpf_path + pa.kfpf_name_ + "/Parameters/Secondary/MassConstraint (2C-Fit)/" + co.kfpf_name_;
//       std::string histo_kfpf_name = kfpf_path + pa.kfpf_name_ + "/Parameters/" + co.kfpf_name_;
      std::string histo_kfpf_name = pa.kfpf_name_ + "/" + co.kfpf_name_;
     
      std::string histo_pfs_name  = pfs_path + pa.pfs_name_ + "/Candidates_" + co.pfs_name_ + "_All_" + pa.pfs_name_;
      
      histo_KFPF = fileInKFPF -> Get<TH1F>(histo_kfpf_name.c_str());
      histo_PFS  = fileInPFS  -> Get<TH1F>(histo_pfs_name.c_str());

      histo_PFS -> SetLineColor(kRed);
      
      CustomizeHistoRanges({histo_PFS, histo_KFPF});
      
      TLegend* leg = new TLegend(0.85, 0.85, 1, 1);
      leg -> AddEntry(histo_PFS, "PFSimple", "L");
      leg -> AddEntry(histo_KFPF, "KFPF", "L");
      
      TCanvas cc("canvas", "canvas", 1500, 900);
      cc.cd();
      histo_PFS -> Draw();
      histo_KFPF -> Draw("same");
      leg->Draw("same");
            
      if(is_first_canvas)
        cc.Print("out.pfs_kfpf.cmpr.pdf(", "pdf");
      else
        cc.Print("out.pfs_kfpf.cmpr.pdf", "pdf");        
    }
    
  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("out.pfs_kfpf.cmpr.pdf)", "pdf");
}

void CustomizeHistoRanges(std::vector<TH1F*> v_histo)
{
  float min =   1e9;
  float max = - 1e9;
  
  for(auto& h : v_histo)
  {
    min = std::min(min, (float)h->GetMinimum());
    max = std::max(max, (float)h->GetMaximum());
  }
  
  if(min>0)
    min = 0.9*min;
  else
    min = 1.1*min;
  
  if(max>0)
    max = 1.1*max;
  else
    max=0.9*max;
  
  for(auto h : v_histo)
    h->GetYaxis()->SetRangeUser(min, max);  
}