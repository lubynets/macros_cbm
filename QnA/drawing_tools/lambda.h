std::string fileName_;
int kSlice=-1;
int kProjection=-1;
int kSelect=-1;

struct Axis
{
  std::string name_;
  std::string sim_name_;
  std::string reco_name_;
  std::string reco_fit_name_;
  std::string title_;
  std::string unit_;
  std::vector<double> bin_edges_;
  float shift_;
  int precision_;
};

std::vector<Axis> axes // 3122, 310
{
  {"pT", "SimParticles_pT", "ReconstructedParticles_pT", "pT", "p_{T}", ", GeV/c", {}, 0., 1},
  {"rapidity", "SimParticles_rapidity", "ReconstructedParticles_rapidity", "y", "y_{CM}", "", {}, -1.62179, 1},
//   {"rapidity", "SimParticles_rapidity", "ReconstructedParticles_rapidity", "y", "y_{CM}", "", {}, -0.985344, 1},
  {"centrality", "RecEventHeader_centrality_tracks", "RecEventHeader_centrality_tracks", "centrality", "centrality", ", %", {}, 0., 0}
};
  
// std::vector<Axis> axes  // 3312
// {
//   {"pT", "SimParticles_pT", "ReconstructedParticles_pT", "ReconstructedParticles_pT", "p_{T}", ", GeV/c", {}, 0., 1},
//   {"rapidity", "SimParticles_rapidity", "ReconstructedParticles_rapidity", "ReconstructedParticles_rapidity", "y_{LAB}", "", {}, -1.62179, 1},
// //   {"rapidity", "SimParticles_rapidity", "ReconstructedParticles_rapidity", "y", "y_{CM}", "", {}, -0.985344, 1},
//   {"centrality", "AnaEventHeader_centrality_tracks", "AnaEventHeader_centrality_tracks", "AnaEventHeader_centrality_tracks", "centrality", ", %", {}, 0., 0}
// };
  
void SetFileName(const std::string& filename){
  fileName_ = filename;
}

void SetAxis(const std::string& name, const std::string& type) {
  int i=0;
  for(auto& ax : axes){
    if(ax.name_ == name){
      if(type == "select") {
        kSelect = i;
      }
      else if(type == "slice") {
        kSlice = i;
      }
      else if(type == "projection") {
        kProjection = i;
      }
      else
        std::cout << "No axis with type " << type << "\n";
      break;
    }
    i++;
  }
  if(i == axes.size())
    std::cout << "No axis with name " << name << "\n";
}

void IntegrateSelectAxis() {
  axes.at(kSelect).bin_edges_ = {axes.at(kSelect).bin_edges_.front(), axes.at(kSelect).bin_edges_.back()};
}

void IntegrateSliceAxis() {
  axes.at(kSlice).bin_edges_ = {axes.at(kSlice).bin_edges_.front(), axes.at(kSlice).bin_edges_.back()};
}

void SetSelectAxisBinEdges(std::vector<double> binedges){ axes.at(kSelect).bin_edges_ = binedges; };
void SetSliceAxisBinEdges(std::vector<double> binedges){ axes.at(kSlice).bin_edges_ = binedges; };
void SetProjectionAxisBinEdges(std::vector<double> binedges){ axes.at(kProjection).bin_edges_ = binedges; };

enum DrawOption {
  kPlain,
  kDifference,
  kChi2,
  kRatio
};

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

namespace Helper {
std::vector<int> palette1{kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed,
                           kOrange+1, kBlue, kGreen+2, kAzure-4, kGray+2, kViolet, kRed};
};

typedef Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>> QnDcSD;
