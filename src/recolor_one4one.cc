//
// author: Wu Shensen
//
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include "include/all_in_one.h"
#include "include/palette-based-recolor.h"

int main(int argc, char **argv) {

  cout.precision(DBL::max_digits10);

  //input path
  string input_image_path = argv[1];
  //output path
  string output_image_path = argv[2];

  int transfer_id = stoi(argv[3]);

  string target_color_args = argv[4];

  Pixel_D target_color;

  boost::regex expr_{"[0-9]{1, 3}"};
  boost::sregex_token_iterator
      parse_(target_color_args.begin(), target_color_args.end(), expr_, 0);
  boost::sregex_token_iterator end_;

  for (; parse_ != end_; ++parse_) {
    string r(parse_->first, parse_->second);
    string g((++parse_)->first, (parse_)->second);
    string b((++parse_)->first, (parse_)->second);
    target_color.x = stoi(r) / 255.0;
    target_color.y = stoi(g) / 255.0;
    target_color.z = stoi(b) / 255.0;
  }



  //like "[72,83,148,124,135,177,171,176,206,220,187,164,232,237,240]"
  //sorted by illuminance in LAB color space
  string current_palette_args = argv[5];

  vector<Pixel_D> current_palette_list;

  boost::regex expr{"[0-9]{1, 3}"};
  boost::sregex_token_iterator
      parse(current_palette_args.begin(), current_palette_args.end(), expr, 0);
  boost::sregex_token_iterator end;

  for (; parse != end; ++parse) {
    string r(parse->first, parse->second);
    string g((++parse)->first, (parse)->second);
    string b((++parse)->first, (parse)->second);
    Pixel_D tmp(stoi(r), stoi(g), stoi(b));
    current_palette_list.push_back(tmp / 255.0);
  }
#ifndef NDEBUG
  cout << "[main] >> input image: " << input_image_path << endl;
  cout << "[main] >> output path: " << output_image_path << endl;
  cout << "[main] >> desired palette:\n";
  for (size_t i = 0; i < current_palette_list.size(); ++i) {
    cout << current_palette_list[i].x << " - "
         << current_palette_list[i].y << " - "
         << current_palette_list[i].z << endl;
  }
  cout << endl;
#endif

//  string input_image_path,
//  string output_image_path,
//  vector<Pixel_D> current_palette_rgb_d,
//  int transfer_id,
//  Pixel_D target_color_rgb_d
  int status = recolor4one(input_image_path, output_image_path,
                           current_palette_list, transfer_id, target_color);
  if (status != 0) {
    return -1;
  }
  return 0;
}