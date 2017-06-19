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

static inline bool isNearWight(const Pixel_255 &rgb) {
  Pixel_D lab;
  RGB2LAB((Pixel_D) rgb / 255.0, lab);
  return lab.x > 80;
}

int main(int argc, char **argv) {

  Json palette;
  //json pipe input without extra single quote
  cin >> palette;

  string inputDir = argv[1];
  string outputDir = argv[2];
  cout << outputDir << endl;

  string image_filename = palette["image_filename"];

  string inputImagePath = inputDir + image_filename;

  if (exists(inputImagePath) && exists(outputDir)) {
    if (!is_regular_file(inputImagePath) || is_directory(inputImagePath)) {
      cerr << inputImagePath << " has something wrong!\n";
      return -1;
    }
  } else {
    cerr << inputImagePath << " or " << outputDir << " does not exist!\n";
    return -1;
  }

  uint K = palette["rgbListSize"];
  vector<Pixel_255> palettes_rgb(K);
  vector<uint> counts(K);
  for (uint i = 0; i < K; ++i) {
    Json rbgAndCount = palette["rgbList"][i];
    palettes_rgb[i].x = rbgAndCount["rgb"][0];
    palettes_rgb[i].y = rbgAndCount["rgb"][1];
    palettes_rgb[i].z = rbgAndCount["rgb"][2];
    counts[i] = rbgAndCount["count"];
  }

  for (size_t j = 0; j < K; ++j) {
    cout << +palettes_rgb[j].x << " - "
         << +palettes_rgb[j].y << " - "
         << +palettes_rgb[j].z << " with "
         << counts[j] << " pixels\n";
  }
  cout << endl;

  Mat input = imread(inputImagePath);

  int height = input.rows;
  int width = input.cols;

  // padding + square_size = (width - padding) / 10
  // padding is 1/6 of square_size
  // line_thickness is 1/10 of padding.
  int padding = width / 71;
  int line_thickness = padding / 10;
  int square_size = 6 * padding;
  int canvas_size = padding + square_size;

  int multi_canvas = (int) ceil(K / 10.0);
  int total_canvas_size = multi_canvas * canvas_size + padding;

  Mat output(Size(height + total_canvas_size, width), input.type());
  copyMakeBorder(input, output, 0, total_canvas_size, 0, 0, BORDER_CONSTANT,
                 Scalar(255, 255, 255));

  for (uint col = 0; col < K; ++col) {
    int row = col / 10;
    int mod_col = col % 10;
    Point start(padding + mod_col * canvas_size,
                height + padding + row * canvas_size);
    rectangle(output, start, start + Point(square_size, square_size),
              Scalar(palettes_rgb[col].z, palettes_rgb[col].y,
              palettes_rgb[col].x), cv::FILLED);

    if (isNearWight(palettes_rgb[col])) {
      rectangle(output, start, start + Point(square_size, square_size),
                Scalar(0, 0, 0), line_thickness);
    }
  }

  //output path, directory must be existed.
  imwrite(outputDir + image_filename, output);
}