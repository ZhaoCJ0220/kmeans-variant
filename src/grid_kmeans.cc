//
// author: Wu Shensen
//

#include <boost/filesystem.hpp>

#include "include/utils.h"
#include "include/json.hpp"

using Json = nlohmann::json;
using namespace boost::filesystem;

typedef struct kmeans_pack {
  double MDL = 0;
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

int main(int argc, char **argv) {

  uint MIN_K_VALUE = 2;
  uint MAX_K_VALUE = 10;

  cout.precision(DBL::max_digits10);

  string input_image_path = argv[1];

  int K = -1;

  if (argc == 3) {
    K = stoi(argv[2]);
  }

  cv::Mat input_image = cv::imread(input_image_path);

  if (input_image.empty()) {
    cerr << "Error >> " << input_image_path << " cannot be loaded!!" << endl;
    return -1;
  }

  uint image_area = (uint) input_image.size().area();

  int height = input_image.rows;
  int width = input_image.cols;

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


  /**
   * run (MAX_K_VALUE - MIN_K_VALUE) time K-means,
   * store them in vector<kmeans_pack>.
   * calculate the Minimum description length (MDL) on each kmeans_pack.
   * choose the best K_VALUE with the least MDL.
   *
   * reference blog:
   * "http://erikerlandson.github.io/blog/2016/08/03/
   * x-medoids-using-minimum-description-length-to-identify-the-k-in-k-medoids/"
   */
  bool stop = false;
  uint limitation = image_area / 100;
  if (K != -1) {
    MIN_K_VALUE = (uint) K;
    MAX_K_VALUE = MIN_K_VALUE;
  }
  for (uint k_value = MIN_K_VALUE; k_value <= MAX_K_VALUE; ++k_value) {
    stop = false;
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

//    for (size_t i = 0; i < res_kmeans.size(); ++i) {
//      cout << "K_VALUE: " << res_kmeans[i].k << endl;
//      for (size_t j = 0; j < res_kmeans[i].palettes_rgb_d.size(); ++j) {
//        cout << round(res_kmeans[i].palettes_rgb_d[j].x * 255) << " - "
//             << round(res_kmeans[i].palettes_rgb_d[j].y * 255) << " - "
//             << round(res_kmeans[i].palettes_rgb_d[j].z * 255) << " with "
//             << res_kmeans[i].counts[j] << " pixels\n";
//      }
//      cout << endl;
//    }

  if (res_kmeans.size() == 0) {
    cerr << "Error >> " << input_image_path << " cannot be kmeans!!" << endl;
    return -1;
  }

  KmeansPack &best_PaletteInfo = res_kmeans[res_kmeans.size() - 1];
//  for (size_t j = 0; j < best_PaletteInfo.palettes_rgb_d.size(); ++j) {
//    cout << round(best_PaletteInfo.palettes_rgb_d[j].x * 255) << " - "
//         << round(best_PaletteInfo.palettes_rgb_d[j].y * 255) << " - "
//         << round(best_PaletteInfo.palettes_rgb_d[j].z * 255) << " with "
//         << best_PaletteInfo.counts[j] << " pixels\n";
//  }
//  cout << endl;

  string filename = basename(input_image_path);
  string file_type = extension(input_image_path);
  file_type.erase(
      remove(file_type.begin(), file_type.end(), '.'), file_type.end());

  //json output
  Json palette;
  palette["colorNum"] = best_PaletteInfo.k;
//  palette["image_serial"] = "PT********";
  palette["image_filename"] = filename + '.' + file_type;
  palette["image_type"] = file_type;
  palette["image_height"] = height;
  palette["image_width"] = width;
  palette["image_area"] = image_area;
  for (size_t i = 0; i < best_PaletteInfo.palettes_rgb_d.size(); ++i) {
    Pixel_D rgb = 255 * best_PaletteInfo.palettes_rgb_d[i];
    vector<uint8_t> rgbList{(uint8_t) round(rgb.x), (uint8_t) round(rgb.y),
                            (uint8_t) round(rgb.z)};
    Json rbgAndCount;
    rbgAndCount["rgb"] = rgbList;
    rbgAndCount["count"] = best_PaletteInfo.counts[i];
    palette["mainColors"][i] = rbgAndCount;
  }
  cout << palette.dump() << endl;

}