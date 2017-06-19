//
// author: Wu Shensen
//

#include "include/all_in_one.h"
#include "include/utils.h"
#include "include/json.hpp"

using namespace std;
using namespace cv;
using Json = nlohmann::json;

typedef struct color_palette {
  size_t K;
  string paletteID;
  vector<uint8_t> rgbList;
} Palette_s;

void sortByLAB(vector<uint8_t> &rgbList, vector<uint8_t> &mess_list) {
  rgbList.clear();
  vector<Pixel_D> unsort_list;
  for (size_t i = 0; i < mess_list.size(); i = i + 3) {
    Pixel_D t(mess_list[i], mess_list[i + 1], mess_list[i + 2]);
    unsort_list.push_back(t / 255.0);
  }

  vector<Pixel_D> unsort_list_lab;
  for (size_t i = 0; i < unsort_list.size(); ++i) {
    Pixel_D t;
    RGB2LAB(unsort_list[i], t);
    unsort_list_lab.push_back(t);
  }

  //sort by L component by LAB
  for (size_t i = 0; i < unsort_list_lab.size(); ++i) {
    for (size_t j = i + 1; j < unsort_list_lab.size(); ++j) {
      if (unsort_list_lab[i].x > unsort_list_lab[j].x) {
        Pixel_D tmp = unsort_list_lab[j];
        unsort_list_lab[j] = unsort_list_lab[i];
        unsort_list_lab[i] = tmp;
      }
    }
  }

  for (size_t i = 0; i < unsort_list_lab.size(); ++i) {
    Pixel_D rgb_double;
    LAB2RGB(unsort_list_lab[i], rgb_double);
    rgbList.push_back((uint8_t) round(rgb_double.x * 255));
    rgbList.push_back((uint8_t) round(rgb_double.y * 255));
    rgbList.push_back((uint8_t) round(rgb_double.z * 255));
  }
}

int main(int argc, char **argv) {
  // read a JSON file
  ifstream input(argv[1]);
  Json palettes_in;
  input >> palettes_in;

  Json palettes_out = Json::array();

  Palette_s palette_s;
  Json palette_s_json;
  for (Json &palette : palettes_in) {
    palette_s.paletteID = palette["paletteID"];
    vector<uint8_t> mess_list = palette["rgbList"];

    sortByLAB(palette_s.rgbList, mess_list);

    palette_s_json["paletteID"] = palette_s.paletteID;
    palette_s_json["rgbList"] = palette_s.rgbList;
    palette_s_json["K"] = palette_s.rgbList.size() / 3;
    palettes_out.push_back(palette_s_json);
  }

  cout << setw(4) << palettes_out << endl;
// write prettified JSON to another file
//  std::ofstream o("pretty.json");
//  o << std::setw(4) << j << std::endl;
//  Json colorList;


  return 0;
}
