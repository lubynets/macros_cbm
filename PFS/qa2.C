void qa2()
{  
  gStyle -> SetOptStat(0);
  
  TFile* fileIn = TFile::Open("/home/user/cbmdir/working/pfsqa/cascade/pfsqa.casc.mcpid_alt.L_optimcuts1.Xi_nocuts.cuts1.root");
  
  struct Particle
  {
    std::string name_;
  };
  
  std::vector<Particle> particles
  {
//     {"Lambda"},
    {"Ksi"}
  };
  
  struct Component
  {
    std::string name_;
  };
  
  std::vector<Component> components
  {
    {"mass"            },
//     {"p"               },
//     {"pT"              },
//     {"rapidity"        },
//     {"phi"             },
//     {"x"               },
//     {"y"               },
//     {"z"               },
    {"generation"      },
    {"chi2_prim_first" },
    {"chi2_prim_second"},
    {"distance"        },        
    {"cosine_first"    },    
    {"cosine_second"   },   
    {"chi2_geo"        },             
    {"l_over_dl"       },      
    {"chi2_topo"       },       
  };
  
  TFile* fileOut = TFile::Open("fileOut.root", "recreate");
  
  TH1F* histoSgnl;
  TH1F* histoBckgr;
  
  for(auto& pa : particles)
    for(auto& co : components)
    {
      std::string histo_sgnl_name  = "Candidates_Sgnl_" + pa.name_ + "/Candidates_" + co.name_ + "_Sgnl_" + pa.name_;
      std::string histo_bckgr_name  = "Candidates_Bckgr_" + pa.name_ + "/Candidates_" + co.name_ + "_Bckgr_" + pa.name_;
      
      histoSgnl = fileIn->Get<TH1F>(histo_sgnl_name.c_str());
      histoBckgr = fileIn->Get<TH1F>(histo_bckgr_name.c_str());
      
      histoBckgr -> Scale(1. * histoSgnl->GetEntries() / histoBckgr->GetEntries());
          
      const float maxSignal = histoSgnl -> GetBinContent(histoSgnl -> GetMaximumBin());
      const float maxBackground = histoBckgr -> GetBinContent(histoBckgr -> GetMaximumBin());
      const float maxGeneral = std::max(maxSignal, maxBackground);
      
      histoSgnl -> SetLineColor(kRed);
      histoSgnl -> SetLineWidth(2);
      histoSgnl -> GetYaxis() -> SetRangeUser(0.1, 1.1*maxGeneral);
      histoSgnl -> GetYaxis() -> SetTitle("Entries");
      histoSgnl -> GetXaxis() -> SetNdivisions(510); // 205
      histoSgnl -> GetYaxis() -> SetNdivisions(205);
      histoSgnl -> GetXaxis() -> SetLabelSize(0.05);
      histoSgnl -> GetXaxis() -> SetTitleSize(0.05);
      histoSgnl -> GetXaxis() -> SetTitleOffset(0.8);
      histoSgnl -> GetYaxis() -> SetLabelSize(0.05);
      histoSgnl -> GetYaxis() -> SetTitleSize(0.05);
      histoSgnl -> GetYaxis() -> SetTitleOffset(0.8);
      histoSgnl -> SetStats(false);
      
      histoBckgr -> SetLineColor(kBlue);
      histoBckgr -> SetLineWidth(2);
      histoBckgr -> GetYaxis() -> SetRangeUser(0.1, 1.1*maxGeneral);
      
      TCanvas* c1 = new TCanvas(co.name_.c_str(), "C1", 1500, 900);
      c1 -> cd();
      
      histoSgnl -> Draw();  
      histoBckgr -> Draw("same");
      
      TLegend *leg = new TLegend(0.75, 0.75, 0.875, 0.875);
      leg -> AddEntry(histoSgnl, "Signal", "L");
      leg -> AddEntry(histoBckgr, "Background", "L");
      leg -> Draw();
      
      c1 -> Write();
      delete c1;      
    }  
  fileOut -> Close();
}
