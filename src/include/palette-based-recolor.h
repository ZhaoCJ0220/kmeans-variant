//
// author: Wu Shensen
//

#ifndef INTELLID_RECOLOR_PALETTE_BASED_RECOLOR_H
#define INTELLID_RECOLOR_PALETTE_BASED_RECOLOR_H

#include "all_in_one.h"

int recolor5(string input_image_path,
             string output_image_path,
             vector<Pixel_D> desired_palette_list);

int recolor4one(string input_image_path,
                string output_image_path,
                vector<Pixel_D> current_palette_rgb_d,
                int transfer_id,
                Pixel_D target_color_rgb_d);

#endif //INTELLID_RECOLOR_PALETTE_BASED_RECOLOR_H
