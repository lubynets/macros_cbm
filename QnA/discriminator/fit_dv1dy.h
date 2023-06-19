#include "Helper.hpp"

std::string fileInPath;
std::string fileInName;
std::vector <TString> object_paths;
std::vector <TString> object_names;
std::vector <int> stat_types;
std::vector <TString> excludedFolders;
std::vector <TString> excludedObjects;
int maxDepth{100};
std::vector<std::string> axestofit;
std::string axistofit;
std::vector<std::string> axestoslice;
std::string axistoslice;
std::vector<std::string> axestoignore;
float midrapidity;
TString fileOutName;

void BuildObjectList(TDirectory *folder, int depth = 0);
bool CheckAxes(Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>* obj);
std::pair<Qn::DataContainerStatDiscriminator, Qn::DataContainerStatDiscriminator> Fit(Qn::DataContainerStatDiscriminator dcIn);
// std::pair<Qn::DataContainerStatDiscriminator, Qn::DataContainerStatDiscriminator> Fit(Qn::DataContainerStatDiscriminator dcIn_1, Qn::DataContainerStatDiscriminator dcIn_2);

std::pair<Qn::DataContainerStatDiscriminator, Qn::DataContainerStatDiscriminator> Fit(Qn::DataContainerStatDiscriminator dcIn_1, Qn::DataContainerStatDiscriminator dcIn_2) {
  Qn::DataContainerStatDiscriminator dcIn = (dcIn_1 + dcIn_2)/2;

  return Fit(dcIn);
}

void BuildObjectList (TDirectory *folder, int depth) {
  TString folder_path = folder -> GetPath();
  folder_path.Remove (0, folder_path.Last (':') + 2);
  TList *keys = folder -> GetListOfKeys ();
  for (auto key : *keys) {
    TString object_name = key -> GetName();
    bool ok = true;
    for(auto& eo : excludedObjects) {
      if(object_name.Contains(eo)) ok = false;
    }
    if(!ok) continue;
    auto object = dynamic_cast <TKey*> (key) -> ReadObj ();
    TString className = object -> ClassName();

    if ( className.Contains ("StatCollect") || className.Contains ("StatCalculate") || className.Contains ("StatDiscriminator") ) {
      object_paths.push_back (folder_path);
      object_names.push_back (object_name);
      if(className.Contains ("StatCollect")) stat_types.push_back(1);
      else if(className.Contains ("StatCalculate")) stat_types.push_back(2);
      else if(className.Contains ("StatDiscriminator")) stat_types.push_back(3);
    }
    else if (className == "TDirectoryFile" && depth < maxDepth) {
      bool skipFolder = false;
      for (auto exFolder : excludedFolders) {
        if (object_name == exFolder) skipFolder = true;
      }
      if (!skipFolder) BuildObjectList (dynamic_cast <TDirectoryFile*> (object), depth + 1);
    }
    delete object;
  }
}

bool CheckAxes(Qn::DataContainer<Qn::StatDiscriminator,Qn::Axis<double>>* obj) {
  bool ok{false};
  std::vector<Qn::AxisD>& axes = obj->GetAxes();
  for(Qn::AxisD& ax : axes) {
    auto name = ax.Name();
    for(auto& atf : axestofit) {
      if(name == atf) {
        axistofit = atf;
        ok = true;
      }
    }
    for(auto& ats : axestoslice) {
      if(name == ats) {
        axistoslice = ats;
      }
    }
    for(auto& ati : axestoignore) {
      if(name == ati) return false;
    }
  }
  return ok;
}
