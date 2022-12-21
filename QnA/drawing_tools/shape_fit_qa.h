int kSlice=-1;
int kProjection=-1;
int kSelect=-1;

struct Axis {
  std::string name_;
  std::string title_;
  std::string unit_;
  std::vector<double> bin_edges_;
  float shift_;
  int precision_;
};

std::vector<Axis> axes // 3122, 310
{
  {"pT", "p_{T}", ", GeV/c", {}, 0., 1},
  {"y", "y_{CM}", "", {}, -1.62179, 1},
  {"centrality", "centrality", ", %", {}, 0., 0}
};

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

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}
