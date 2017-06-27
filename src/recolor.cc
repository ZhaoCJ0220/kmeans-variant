//
// author: Wu Shensen
//
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include "include/all_in_one.h"
#include "include/palette-based-recolor.h"

#include "include/utils.h"
#include "include/json.hpp"

using Json = nlohmann::json;
using namespace boost::filesystem;

typedef struct kmeans_pack {
  double MDL = 0;//?
  uint k = 0;
  vector<uint> counts;
  vector<Pixel_D> palettes_rgb_d;
} KmeansPack;

int getPaletteInfo(const vector<Pixel_D> &pixels_rgb_d,
                   KmeansPack &paletteInfo, uint k) {

  paletteInfo.palettes_rgb_d.resize(k);
  paletteInfo.counts.resize(k);

  doKmeans(pixels_rgb_d, paletteInfo.palettes_rgb_d, paletteInfo.counts);

  return 0;
}

int getPalette(string input_image_path) {

  uint MIN_K_VALUE = 2;
  uint MAX_K_VALUE = 10;

  cout.precision(DBL::max_digits10);

  cv::Mat input_image = cv::imread(input_image_path);

  if (input_image.empty()) {
    cerr << "Error >> " << input_image_path << " cannot be loaded!!" << endl;
    return -1;
  }

  uint image_area = (uint) input_image.size().area();

//  int height = input_image.rows;
//  int width = input_image.cols;

  vector<Pixel_D> pixels_rgb_d(image_area);

  //read b-g-r three channels from opencv Mat
  for (uint i = 0; i < image_area; ++i) {
    //r
    pixels_rgb_d[i].x = input_image.data[3 * i + 2] / 255.0;
    //g
    pixels_rgb_d[i].y = input_image.data[3 * i + 1] / 255.0;
    //b
    pixels_rgb_d[i].z = input_image.data[3 * i] / 255.0;
  }

  vector<KmeansPack> res_kmeans;

  bool stop = false;
  uint limitation = image_area / 100;
  for (uint k_value = MIN_K_VALUE; k_value <= MAX_K_VALUE; ++k_value) {
    KmeansPack paletteInfo;
    getPaletteInfo(pixels_rgb_d, paletteInfo, k_value);
    for (auto &cnt : paletteInfo.counts) {
      if (res_kmeans.size() != 0 && cnt < limitation) {
        stop = true;
        break;
      }
    }
    if (stop) {
      break;
    }
    paletteInfo.k = k_value;
    res_kmeans.push_back(paletteInfo);
  }

  if (res_kmeans.size() == 0) {
    cerr << "Error >> " << input_image_path << " cannot be kmeans!!" << endl;
    return -1;
  }

  KmeansPack &best_PaletteInfo = res_kmeans[res_kmeans.size() - 1];

  string filename = basename(input_image_path);
  string file_type = extension(input_image_path);
  file_type.erase(
      remove(file_type.begin(), file_type.end(), '.'), file_type.end());

  //json output
  Json palette;
  palette["colorNum"] = best_PaletteInfo.k;
//  palette["image_serial"] = "PT********";
//  palette["image_filename"] = filename + '.' + file_type;
//  palette["image_type"] = file_type;
//  palette["image_height"] = height;
//  palette["image_width"] = width;
//  palette["image_area"] = image_area;
  for (size_t i = 0; i < best_PaletteInfo.palettes_rgb_d.size(); ++i) {
    Pixel_D rgb = 255 * best_PaletteInfo.palettes_rgb_d[i];
    vector<uint8_t> rgbList{(uint8_t) round(rgb.x), (uint8_t) round(rgb.y),
                            (uint8_t) round(rgb.z)};
    Json rbgAndCount;
    rbgAndCount["rgb"] = rgbList;
    rbgAndCount["count"] = best_PaletteInfo.counts[i];
    palette["mainColors"][i] = rbgAndCount;
  }
  cout << "@@@" << palette.dump() << "@@@" << endl;
  return 0;
}

int main(int argc, char **argv) {

  cout.precision(DBL::max_digits10);//输出精度

  //input path
  string input_image_path = argv[1];
  //output path
  string output_image_path = argv[2];
  //like "[72,83,148,124,135,177,171,176,206,220,187,164,232,237,240]"
  //sorted by illuminance in LAB color space
  string desired_palette_args = argv[3];

  int with_palette = stoi(argv[4]);//string to int

  vector<Pixel_D> desired_palette_list;

  boost::regex expr{"[0-9]{1, 3}"};
  boost::sregex_token_iterator
      parse(desired_palette_args.begin(), desired_palette_args.end(), expr, 0);
  boost::sregex_token_iterator end;

  for (; parse != end; ++parse) {
    string r(parse->first, parse->second);
    string g((++parse)->first, (parse)->second);
    string b((++parse)->first, (parse)->second);
    Pixel_D tmp(stoi(r), stoi(g), stoi(b));
    desired_palette_list.push_back(tmp / 255.0);
  }
#ifndef NDEBUG
  cout << "[main] >> input image: " << input_image_path << endl;
  cout << "[main] >> output path: " << output_image_path << endl;
  cout << "[main] >> desired palette:\n";
  for (size_t i = 0; i < desired_palette_list.size(); ++i) {
    cout << desired_palette_list[i].x << " - "
         << desired_palette_list[i].y << " - "
         << desired_palette_list[i].z << endl;
  }//desired_palette_list是目标色list
  cout << endl;
#endif

  int status = recolor5(input_image_path, output_image_path, desired_palette_list);
  if (status != 0) {
    return -1;
  }
  if (with_palette == 1) {
    cout << "start get palette:" << endl;
    return getPalette(output_image_path);
  }
  return 0;
}