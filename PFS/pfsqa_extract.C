void pfsqa_extract()
{
  gStyle -> SetOptStat(0);
  
  TFile* fileIn = TFile::Open("/home/user/cbmdir/working/pfsqa/cascade/pfsqa.casc.mcpid_alt.def2.lambda.noextrapol.root");
    
  struct Particle
  {
    std::string name_;
  };
  
  std::vector<Particle> particles
  {
    {"Lambda"},
//     {"Ksi"}
  };
  
  struct SgnlBckgr
  {
    std::string name_;
  };
  
  std::vector<SgnlBckgr> sgnlbckgr
  {
//     {"All"},
    {"Sgnl"},
//     {"Bckgr"}
  };
  
  struct Component
  {
    std::string name_;
    int dim_;
  };
  
  std::vector<Component> components
  {
//     {"mass"            },
// //     {"p"               },
// //     {"pT"              },
// //     {"rapidity"        },
// //     {"phi"             },
// //     {"x"               },
// //     {"y"               },
// //     {"z"               },
//     {"generation"      },
//     {"chi2_prim_first" },
//     {"chi2_prim_second"},
//     {"distance"        },        
//     {"cosine_first"    },    
//     {"cosine_second"   },   
//     {"chi2_geo"        },             
//     {"l_over_dl"       },      
//     {"chi2_topo"       },       
//     {"chi2prim_mother" }, 
//     {"invmass_discr"   },
    {"Simulated_p_Candidates_p",               2},
    {"Simulated_pT_Candidates_pT",             2},
    {"Simulated_px_Candidates_px",             2},
    {"Simulated_py_Candidates_py",             2},
    {"Simulated_pz_Candidates_pz",             2},
    {"Simulated_rapidity_Candidates_rapidity", 2},
    {"Simulated_phi_Candidates_phi",           2},
    {"diff_p",                                 1},
    {"diff_px",                                1},
    {"diff_py",                                1},
    {"diff_pz",                                1},
    {"diff_pT",                                1},
    {"diff_rapidity",                          1},
    {"diff_phi" ,                              1},
  };
    
  TH1* histo{nullptr};
  
  bool is_first_canvas = true;
    
  for(auto& sb : sgnlbckgr)
    for(auto& pa : particles)
      for(auto& co : components)
      {
//         std::string histo_name  = "Candidates_" + sb.name_ + "_" + pa.name_ + "/Candidates_" + co.name_ + "_" + sb.name_ + "_" + pa.name_ + "";
        std::string histo_name  = "Candidates_Simulated_" + sb.name_ + "_" + pa.name_ + "/" + co.name_ + "_" + sb.name_ + "_" + pa.name_;
                
        if(co.dim_ == 1)
          histo  = fileIn  -> Get<TH1F>(histo_name.c_str());
        else if(co.dim_ == 2)
          histo  = fileIn  -> Get<TH2F>(histo_name.c_str());
        
        histo -> SetTitle(co.name_.c_str());
        
  //       CustomizeHistoRanges({histo});
        
        TCanvas cc("canvas", "canvas", 1500, 900);
        cc.cd();
        cc.SetLogz();
        histo -> Draw("colz");
              
        if(is_first_canvas)
          cc.Print("out.pfs.cplx.pdf(", "pdf");
        else
          cc.Print("out.pfs.cplx.pdf", "pdf");        
      }
    
  TCanvas emptycanvas("", "", 1500, 900);
  emptycanvas.Print("out.pfs.cplx.pdf)", "pdf");
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