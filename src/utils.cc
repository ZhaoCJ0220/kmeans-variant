//
// author: Wu Shensen
//
#include <math.h>
#include "include/all_in_one.h"
#include "include/utils.h"

double LAB2Distance(const Pixel_D &a_point, const Pixel_D &b_point) {
  // 以下是 LCH 颜色空间的距离计算方式
  //  for textiles k1 0.048， k2 0.014
  double K1 = 0.048;
  double K2 = 0.014;
  double del_L = a_point.x - b_point.x;
  double c1 = sqrt(a_point.y * a_point.y + a_point.z * a_point.z);
  double c2 = sqrt(b_point.y * b_point.y + b_point.z * b_point.z);
  double c_ab = c1 - c2;
  double h_ab = (a_point.y - b_point.y) * (a_point.y - b_point.y)
      + (a_point.z - b_point.z) * (a_point.z - b_point.z)
      - c_ab * c_ab;
  return del_L * del_L + c_ab * c_ab / (1 + K1 * c1) / (1 + K1 * c1)
      + h_ab / (1 + K2 * c1) / (1 + K2 * c1);
}

double smoothL(double delta_L_absolute, double delta_L_relative) {

  double lambda = 0.2 * log(2);
  double smooth = log(exp(lambda * delta_L_absolute)
                          + exp(lambda * delta_L_relative) - 1);
  return smooth / lambda - delta_L_absolute;
}

double Norm2Distance(const Pixel_D &a_point, const Pixel_D &b_point) {
  double delta_x = a_point.x - b_point.x;
  double delta_y = a_point.y - b_point.y;
  double delta_z = a_point.z - b_point.z;
  return (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
}

void LAB2RGB(const Pixel_D &lab, Pixel_D &pixel_rgb_tmp) {

  double r, g, b, x, y, z, fx, fy, fz;
  fy = (lab.x + 16.0f) / 116.0f;
  fx = lab.y / 500.0f + fy;
  fz = fy - lab.z / 200.0f;

  x = pow(fx, 3.0) > 0.008856f ? pow(fx, 3.0) : fx / 7.787f - 0.017713f;
  y = pow(fy, 3.0) > 0.008856f ? pow(fy, 3.0) : fy / 7.787f - 0.017713f;
  z = pow(fz, 3.0) > 0.008856f ? pow(fz, 3.0) : fz / 7.787f - 0.017713f;
  x = 0.95047f * x;
  z = 1.08883f * z;

  r = x * 3.2406f + y * -1.5372f + z * -0.4986f;
  g = x * -0.9689f + y * 1.8758f + z * 0.0415f;
  b = x * 0.0557f + y * -0.2040f + z * 1.0570f;

  r = r > 0.0031308 ? (1.055f * pow(r, 0.416667) - 0.055f) : 12.92f * r;
  g = g > 0.0031308 ? (1.055f * pow(g, 0.416667) - 0.055f) : 12.92f * g;
  b = b > 0.0031308 ? (1.055f * pow(b, 0.416667) - 0.055f) : 12.92f * b;

  //the value could be out of the RGB gamut
  //you must post-process the value validation outside.
  pixel_rgb_tmp.x = r;
  pixel_rgb_tmp.y = g;
  pixel_rgb_tmp.z = b;

}

void RGB2LAB(const Pixel_D &rgb, Pixel_D &pixel_lab_tmp) {

  double r, g, b, x, y, z, fx, fy, fz;
  r = rgb.x;
  g = rgb.y;
  b = rgb.z;

  r = r > 0.04045 ? pow((r + 0.055) / 1.055f, 2.4) : r / 12.92f;
  g = g > 0.04045 ? pow((g + 0.055) / 1.055f, 2.4) : g / 12.92f;
  b = b > 0.04045 ? pow((b + 0.055) / 1.055f, 2.4) : b / 12.92f;

  //Observer. = 2°, Illuminant = D65
  x = (r * 0.4124f + g * 0.3576f + b * 0.1805f) * 100.0f / 95.047f;
  y = r * 0.2126f + g * 0.7152f + b * 0.0722f;
  z = (r * 0.0193f + g * 0.1192f + b * 0.9505f) * 100.0f / 108.883f;

  //xyz to lab
  fx = x > 0.008856f ? pow(x, 0.333333) : (7.787f * x + 0.137931f);
  fy = y > 0.008856f ? pow(y, 0.333333) : (7.787f * y + 0.137931f);
  fz = z > 0.008856f ? pow(z, 0.333333) : (7.787f * z + 0.137931f);

  pixel_lab_tmp.x = 116.0f * fy - 16.0f;
  pixel_lab_tmp.y = 500.0f * (fx - fy);
  pixel_lab_tmp.z = 200.0f * (fy - fz);
}

double RGB2CMYK(const Pixel_D &rgb,
              CMYK_Pix_D &pixel_cmyk_tmp,double rate){
  double r, g, b,c,m,y,k;
  r = rgb.x;
  g = rgb.y;
  b = rgb.z;
  k = 1-max(max(r,g),b);
  if(k==1){
    c=0;
    m=0;
    y=0;
  }
  else{
    c = (1-k-r)/(1-k);
    m = (1-k-g)/(1-k);
    y = (1-k-b)/(1-k);
  }

  pixel_cmyk_tmp.at<double>(0) = c;
  pixel_cmyk_tmp.at<double>(1) = m;
  pixel_cmyk_tmp.at<double>(2) = y;
  pixel_cmyk_tmp.at<double>(3) = k;
  return (c+m+y+rate*k)/(3+rate);
}

int gridKmeans(const vector<Pixel_D> &pixels_rgb_d,
               vector<Pixel_D> &palette_lab,
               vector<uint> &palette_count,
               uint grid_num = 16) {

  const uint K = (uint) palette_count.size();
  const size_t image_area = pixels_rgb_d.size();
  const uint grid_num_square = grid_num * grid_num;
  const uint grid_size = grid_num_square * grid_num;
  const double step_size = 1.0 / (grid_num - 1);

  /**
   * in each grid, count the number the pixels which fall into,
   * and calculate the LAB channels sum of these pixels.
   */
  vector<uint> pixels_count(grid_size);//一维下标对应bin中的像素个数
  vector<Pixel_D> pixels_lab_sum(grid_size);//lab三通道的和值

  for (size_t i = 0; i < image_area; ++i) {

    Pixel_D pixel_lab_tmp;
    uint bin, bin_r, bin_g, bin_b;

    //pixels_rgb_d [0,1]
    RGB2LAB(pixels_rgb_d[i], pixel_lab_tmp);

    bin_r = (uint) round(pixels_rgb_d[i].x / step_size);
    bin_g = (uint) round(pixels_rgb_d[i].y / step_size);
    bin_b = (uint) round(pixels_rgb_d[i].z / step_size);
    bin = bin_r * grid_num_square
        + bin_g * grid_num
        + bin_b;//取到一维下标

    pixels_count[bin] += 1;
    //LAB
    pixels_lab_sum[bin] += pixel_lab_tmp;
  }


  /**
   * rearrange the data, remove unused grid
   */
  size_t valid_grid_size = 0;
  for (size_t i = 0; i < grid_size; ++i) {
    if (pixels_count[i] > 0) {
      ++valid_grid_size;
    }
  }

  vector<uint> valid_pixels_count(valid_grid_size);//有效bin中像素个数
  vector<Pixel_D> valid_pixels_lab_mean(valid_grid_size);//lab均值
  valid_grid_size = 0;

  //remove unused grids
  for (size_t i = 0; i < grid_size; ++i) {
    uint count = pixels_count[i];
    if (count > 0) {
      valid_pixels_count[valid_grid_size] = count;
      valid_pixels_lab_mean[valid_grid_size] = pixels_lab_sum[i] /
          (double) count;
      ++valid_grid_size;
    }
  }

  {
    /**
     * choose the CENTER_NUM initial K-means point
     */
    vector<uint> valid_pixels_count_copy(valid_pixels_count);

    for (size_t i = 0; i < K; ++i) {
      // 样本数最多的颜色格子（取了个下标，为啥要auto？）
      auto idx = distance(valid_pixels_count_copy.begin(),
                          max_element(valid_pixels_count_copy.begin(),
                                      valid_pixels_count_copy.end()));
      assert(idx >= 0);

      palette_lab[i] = valid_pixels_lab_mean[idx];

      for (size_t j = 0; j < valid_grid_size; ++j) {
        double distance =
            Norm2Distance(valid_pixels_lab_mean[j], palette_lab[i]);
        distance = distance / (80 * 80);
        // 如果离当前K-mean起始点越远
        // 在下一次选择K-mean起始点中
        // 该颜色格子的样本数衰减系数会较低
        valid_pixels_count_copy[j] =
            (uint) (valid_pixels_count_copy[j] * (1 - exp(-distance)));
      }
    }
  }

  //used for count the pixels which fall into each kmeans group
  vector<uint> kmeans_count(K);

  //kmeans 20 iterations
  for (size_t iter = 0; iter < 20; ++iter) {

    fill(kmeans_count.begin(), kmeans_count.end(), 0);
    vector<Pixel_D> kmeans_sum(K);

    for (size_t i = 0; i < valid_grid_size; ++i) {

      int min_idx = -1;
      double min_dis = FLT_MAX;

      for (uint j = 0; j < K; ++j) {

        double
            distance = Norm2Distance(valid_pixels_lab_mean[i], palette_lab[j]);
        if (distance < min_dis) {
          min_dis = distance;
          min_idx = j;
        }
      }//更新到min_idx中
      kmeans_count[min_idx] += valid_pixels_count[i];
      kmeans_sum[min_idx] +=
          valid_pixels_lab_mean[i] * (double) valid_pixels_count[i];
    }

    for (size_t i = 0; i < K; ++i) {
      if (kmeans_count[i] > 0) {
        palette_lab[i] = kmeans_sum[i] / (double) kmeans_count[i];
        palette_count[i] = kmeans_count[i];
      }
    }
  }//end --- kmeans iterations
  return 0;
}

/**
 *
 * @param pixels_rgb_d image pixels data, gamut is [0, 1]
 * @param palette_rgb_d get Kmeans rgb color,
 * @param palette_count
 * @return
 */
int doKmeans(const vector<Pixel_D> &pixels_rgb_d,
             vector<Pixel_D> &palette_rgb_d,
             vector<uint> &palette_count,double rate) {

  assert(palette_rgb_d.size() == palette_count.size());//不等则程序终止

  uint K = (uint) palette_count.size();
  vector<Pixel_D> palette_lab(K);

  gridKmeans(pixels_rgb_d, palette_lab, palette_count, 128);//grid remain to search

  //sort by L component by LAB(按L选择排序（升序））
//  for (long i = 0; i < K; ++i) {
//    for (long j = i + 1; j < K; ++j) {
//      if (palette_lab[i].x > palette_lab[j].x) {
//        Pixel_D tmp = palette_lab[j];
//        palette_lab[j] = palette_lab[i];
//        palette_lab[i] = tmp;
//        uint count_tmp = palette_count[j];
//        palette_count[j] = palette_count[i];
//        palette_count[i] = count_tmp;
//      }
//    }
//  }


  for (size_t i = 0; i < K; ++i) {
    Pixel_D rgb_double;
    LAB2RGB(palette_lab[i], rgb_double);
    palette_rgb_d[i].x = max(0.0, min(1.0, rgb_double.x));
    palette_rgb_d[i].y = max(0.0, min(1.0, rgb_double.y));
    palette_rgb_d[i].z = max(0.0, min(1.0, rgb_double.z));
  }//lab转rgb？

  //排序
//  vector<CMYK_Pix_D> palette_cmyk(K);
//  for (long i = 0; i < K; ++i) {
//    CMYK_Pix_D temp_cmyk(4,1,CV_64F);
//    palette_cmyk[i] = temp_cmyk;
//  }
//
//  vector<double> depth(K);
//  for (long i = 0; i < K; ++i) {
//    depth[i] = RGB2CMYK(palette_rgb_d[i],palette_cmyk[i],rate);
//  }
//  for (long i = 0; i < K; ++i) {
//    for (long j = i + 1; j < K; ++j) {
//      if (depth[i]> depth[j]) {
//        Pixel_D tmp = palette_rgb_d[j];
//        palette_rgb_d[j] = palette_rgb_d[i];
//        palette_rgb_d[i] = tmp;
//        uint count_tmp = palette_count[j];
//        palette_count[j] = palette_count[i];
//        palette_count[i] = count_tmp;
//      }
//    }
//  }
  return 0;
}