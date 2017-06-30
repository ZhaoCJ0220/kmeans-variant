//
// author: Wu Shensen
//

#ifndef INTELLID_RECOLOR_UTILS_H
#define INTELLID_RECOLOR_UTILS_H

#include "all_in_one.h"

double LAB2Distance(const Pixel_D &a_point,
                    const Pixel_D &b_point);

double smoothL(double delta_L_absolute,
               double delta_L_relative);

double Norm2Distance(const Pixel_D &a_point,
                     const Pixel_D &b_point);

void LAB2RGB(const Pixel_D &lab,
             Pixel_D &pixel_rgb_tmp);

void RGB2LAB(const Pixel_D &rgb,
             Pixel_D &pixel_lab_tmp);

double RGB2CMYK(const Pixel_D &rgb,
              CMYK_Pix_D &pixel_cmyk_tmp,double rate);

int gridKmeans(const vector<Pixel_D> &pixels_rgb_d,
               vector<Pixel_D> &palette_lab,
               vector<uint> &palette_count,
               uint grid_num);

int doKmeans(const vector<Pixel_D> &pixels_rgb_d,
             vector<Pixel_D> &palette_rgb_d,
             vector<uint> &palette_count,double rate);

#endif //INTELLID_RECOLOR_UTILS_H
