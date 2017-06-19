//
// author: Wu Shensen
//

#include <boost/filesystem.hpp>

#include "include/all_in_one.h"
#include "include/json.hpp"
#include "include/utils.h"

using namespace std;
using namespace cv;
using namespace boost::filesystem;
using Json = nlohmann::json;

void sortByLAB(vector<Pixel_255> &rgbList, const vector<Pixel_255> &mess_list) {
  rgbList.clear();
  vector<Pixel_D> unsort_list;
  for (size_t i = 0; i < mess_list.size(); ++i) {
    unsort_list.push_back((Pixel_D) mess_list[i] / 255.0);
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
    Pixel_255 rgb((uint8_t) round(rgb_double.x * 255),
                  (uint8_t) round(rgb_double.y * 255),
                  (uint8_t) round(rgb_double.z * 255));
    rgbList.push_back(rgb);
  }
}

int main(int argc, char **argv) {
  uint64 K = (uint64) stoi(argv[1]);
  string inputImagePath = argv[2];

  Mat input = imread(inputImagePath);

  if (input.empty()) {
    cerr << "Error >> " << inputImagePath << " cannot be loaded!!" << endl;
    return -1;
  }

  int height = input.rows;
  int width = input.cols;

//  cout << height << "x" << width << endl;

  // padding + square_size = (width - padding) / 10
  // padding is 1/6 of square_size
  // line_thickness is 1/10 of padding.
  int padding = width / 71;
//  int line_thickness = padding / 10;
  int square_size = 6 * padding;
  int canvas_size = padding + square_size;

  int multi_canvas = (int) ceil(K / 10.0);
  int total_canvas_size = multi_canvas * canvas_size + padding;

  height = height - total_canvas_size;
//  cout << "original:" << height << "x" << width << endl;
//
//  cout << input.at<Vec3b>(20, 420) << endl;

  vector<Pixel_255> palettes_rgb(K);

  for (uint col = 0; col < K; ++col) {
    int row = col / 10;
    int mod_col = col % 10;
    int x_start = padding + mod_col * canvas_size;
    int y_start = height + padding + row * canvas_size;
    int x = x_start + square_size / 2;
    int y = y_start + square_size / 2;
    //R
    palettes_rgb[col].x = input.at<Vec3b>(Point(x, y))[2];
    //G
    palettes_rgb[col].y = input.at<Vec3b>(Point(x, y))[1];
    //B
    palettes_rgb[col].z = input.at<Vec3b>(Point(x, y))[0];
//    cout << x << "_" << y << endl;
//    cout << +palettes_rgb[col].x << "_"
//         << +palettes_rgb[col].y << "_"
//         << +palettes_rgb[col].z << endl << endl;
  }

  vector<Pixel_255> sorted_palettes_rgb;

  sortByLAB(sorted_palettes_rgb, palettes_rgb);

  //{"rgbList":[{"rgb":[105,127,164]},{"rgb":[154,189,221]},{"rgb":[152,194,192]},
  // {"rgb":[201,216,209]},{"rgb":[234,231,226]}],"rgbListSize":5}
  Json palette;

  palette["rgbListSize"] = K;
  for (uint i = 0; i < K; ++i) {
    vector<uint8_t> rgb3Point{sorted_palettes_rgb[i].x,
                              sorted_palettes_rgb[i].y,
                              sorted_palettes_rgb[i].z};
    Json rbg;
    rbg["rgb"] = rgb3Point;
    palette["rgbList"][i] = rbg;
  }

  cout << palette << endl;
}