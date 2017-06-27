//
// author: Wu Shensen
//
#include "include/all_in_one.h"
#include "include/utils.h"

using namespace std;

void cal_rbf_single_grid(const Pixel_D &rgb, vector<double> &res,
                         const uint RBF_GRID_NUM) {
  double step_1 = RBF_GRID_NUM + 1;
  double step_2 = step_1 * step_1;
  double r_bias, g_bias, b_bias;
  double main_base_grid_idx;
  double r_idx, g_idx, b_idx;

  r_idx = rgb.x * RBF_GRID_NUM;
  r_bias = r_idx - floor(r_idx);
  r_idx = floor(r_idx);
  if (r_idx == RBF_GRID_NUM) {
    r_idx = RBF_GRID_NUM - 1;
    r_bias = 1;
  }//?

  g_idx = rgb.y * RBF_GRID_NUM;
  g_bias = g_idx - floor(g_idx);
  g_idx = floor(g_idx);
  if (g_idx == RBF_GRID_NUM) {
    g_idx = RBF_GRID_NUM - 1;
    g_bias = 1;
  }

  b_idx = rgb.z * RBF_GRID_NUM;
  b_bias = b_idx - floor(b_idx);
  b_idx = floor(b_idx);
  if (b_idx == RBF_GRID_NUM) {
    b_idx = RBF_GRID_NUM - 1;
    b_bias = 1;
  }

  main_base_grid_idx = r_idx + g_idx * step_1 + b_idx * step_2;//这里为什么是13？

  res[0] = main_base_grid_idx;
  res[1] = main_base_grid_idx + step_2;
  res[2] = main_base_grid_idx + step_1;
  res[3] = main_base_grid_idx + step_1 + step_2;
  res[4] = main_base_grid_idx + 1;
  res[5] = main_base_grid_idx + step_2 + 1;
  res[6] = main_base_grid_idx + step_1 + 1;
  res[7] = main_base_grid_idx + step_1 + step_2 + 1;

  res[8] = (1 - r_bias) * (1 - g_bias) * (1 - b_bias);
  res[9] = (1 - r_bias) * (1 - g_bias) * b_bias;
  res[10] = (1 - r_bias) * g_bias * (1 - b_bias);
  res[11] = (1 - r_bias) * g_bias * b_bias;
  res[12] = r_bias * (1 - g_bias) * (1 - b_bias);
  res[13] = r_bias * (1 - g_bias) * b_bias;
  res[14] = r_bias * g_bias * (1 - b_bias);
  res[15] = r_bias * g_bias * b_bias;
}

void prepare_rbf_color(const vector<Pixel_D> &pixels_rgb_d,
                       vector<Pixel_D> &rbf_lab_colors,
                       vector<vector<double>> &rbf_weight_index,
                       vector<vector<double>> &rbf_weight_map,
                       const uint RBF_GRID_NUM) {

  const double step = 1.0 / RBF_GRID_NUM;
  uint iterations = RBF_GRID_NUM + 1;
  uint idx = 0;
  for (uint i = 0; i < iterations; ++i) {//b-channels
    for (uint j = 0; j < iterations; ++j) {//g-channels
      for (uint k = 0; k < iterations; ++k) {//r-channels
        //反序索引累加, k, j, i，只是按照rgb的顺序而已，没有骚操作
        RGB2LAB(Pixel_D(k * step, j * step, i * step), rbf_lab_colors[idx]);
        ++idx;
      }
    }
  }//分成13^3个顶点的值

  const size_t image_area = pixels_rgb_d.size();
  //magic data structure from original js code for rbf color calculation.
  vector<double> res(16);
  for (size_t i = 0; i < image_area; i++) {

    cal_rbf_single_grid(pixels_rgb_d[i], res, RBF_GRID_NUM);

    for (uint8_t j = 0; j < 8; j++) {
      rbf_weight_index[i][j] = res[j];
      rbf_weight_map[i][j] = res[j + 8];//res[0]和res[8]对应？？感觉是反的
	  
    }
  }
}

void modifyLuminance(vector<Pixel_D> current_palette_lab_copy,
                     vector<Pixel_D> &target_palette_rgb_d,
                     const uint transfer_id,
                     const Pixel_D &target_color_rgb) {

  uint CENTER_NUM_WITH_BLACK = (uint) target_palette_rgb_d.size();

  Pixel_D target_color_lab;
  RGB2LAB(target_color_rgb, target_color_lab);

  // 因保留了BLACK色，需要做出偏移
  double old_L = current_palette_lab_copy[transfer_id + 1].x;

  double new_L = target_color_lab.x;
  double delta_L_absolute = old_L - new_L;

  for (uint i = 1; i < CENTER_NUM_WITH_BLACK; ++i) {

    if (i != (transfer_id + 1)) {
      if (current_palette_lab_copy[i].x < old_L) {

        current_palette_lab_copy[i].x = new_L
            - smoothL(delta_L_absolute,
                      old_L - current_palette_lab_copy[i].x);
      } else {

        current_palette_lab_copy[i].x = new_L
            + smoothL(-delta_L_absolute,
                      current_palette_lab_copy[i].x - old_L);
      }
    }
  }

  current_palette_lab_copy[transfer_id + 1] = target_color_lab;

  for (size_t i = 0; i < CENTER_NUM_WITH_BLACK; ++i) {
    cout << current_palette_lab_copy[i].x << " - "
         << current_palette_lab_copy[i].y << " - "
         << current_palette_lab_copy[i].z << endl;
  }
  cout << endl;

  for (uint i = 0; i < CENTER_NUM_WITH_BLACK; ++i) {

    Pixel_D rgb_double;
    LAB2RGB(current_palette_lab_copy[i], rgb_double);
    target_palette_rgb_d[i].x = max(0.0, min(1.0, rgb_double.x));
    target_palette_rgb_d[i].y = max(0.0, min(1.0, rgb_double.y));
    target_palette_rgb_d[i].z = max(0.0, min(1.0, rgb_double.z));
  }
}

bool isBigChange(Pixel_D &point) {
  return (abs(point.x) > 0.5 || abs(point.y) > 0.5 || abs(point.z) > 0.5);
}

bool isOutOfBoundary(Pixel_D &point) {
  double out_threshold = 0.5;
  return (point.x < -out_threshold || point.x > 255 + out_threshold
      || point.y < -out_threshold || point.y > 255 + out_threshold
      || point.z < -out_threshold || point.z > 255 + out_threshold);
}

double findBoundary(const Pixel_D &point, const Pixel_D &direction,
                    double start, double end) {
  //Assume dir is large
  double middle;
  for (uint8_t iter = 0; iter < 15; ++iter) {
    middle = 0.5 * (start + end);
    Pixel_D testRGB;
    LAB2RGB((point + middle * direction), testRGB);
    if (isOutOfBoundary(testRGB)) {
      end = middle;
    } else {
      start = middle;
    }
  }
  return start;
}//在lab空间上找到对应的rgb空间的boundary

Pixel_D calculateSinglePoint(vector<Pixel_D> &current_palette_lab,
                             double rbf_param,
                             Pixel_D &in_color_rbf,
                             vector<Pixel_D> &diffs_rbf) {

  uint CENTER_NUM = (uint) current_palette_lab.size();

  // weightMat: e^current_palette_lab distance
  cv::Mat weightMat(CENTER_NUM, CENTER_NUM, CV_64F);//新建k*k的矩阵
  for (uint i = 0; i < CENTER_NUM; ++i) {
    for (uint j = 0; j < CENTER_NUM; ++j) {
      //首元素为黑，做偏移
      weightMat.at<double>(i, j) = exp(rbf_param * -1 *
          LAB2Distance(current_palette_lab[i], current_palette_lab[j]));//LAB2Distance返回的是距离平方？rbf_param为1/(2*(6r)^2)?
    }
  }

  cv::Mat diffMat(CENTER_NUM, 1, CV_64F);//k*1
  for (uint i = 0; i < CENTER_NUM; ++i) {
    //首元素为黑，做偏移?没做
    diffMat.at<double>(i) = exp(rbf_param * -1 *
        LAB2Distance(current_palette_lab[i], in_color_rbf));
  }

  cv::Mat paramMat = weightMat.inv(cv::DECOMP_LU) * diffMat;//weightMat T *diffMat（对称矩阵转置干啥？）

  Pixel_D delta;
  double scale = 0;
  for (uint i = 0; i < CENTER_NUM; ++i) {
    scale += max(paramMat.at<double>(i), 0.0);
  }
  for (uint i = 0; i < CENTER_NUM; ++i) {
    if (paramMat.at<double>(i) > 0) {
      delta += paramMat.at<double>(i) / scale * diffs_rbf[i];//Wi(X)=paramMat.at<double>(i) / scale???
    }
  }
  return in_color_rbf + delta;
}

void calculateGridResult(vector<Pixel_D> &current_palette_lab,
                         vector<Pixel_D> &palette_diffs,
                         const vector<Pixel_D> &rbf_lab_colors,
                         vector<Pixel_D> &rbf_RGB_grids) {
  const uint RBF_PARAM_COFF = 5;

  uint CENTER_NUM = (uint) current_palette_lab.size();

  double average_distance = 0;
  uint distance_count = 0;
  //首元素为黑，做偏移
  for (uint i = 0; i < CENTER_NUM; ++i) {
    for (uint j = i + 1; j < CENTER_NUM; ++j) {
      ++distance_count;
      average_distance += sqrt(LAB2Distance(current_palette_lab[i],
                                            current_palette_lab[j]));
    }
  }

  if (CENTER_NUM > 1) {
    average_distance /= distance_count;
  } else {
    average_distance = 1.0;
  }

  double rbf_param = RBF_PARAM_COFF / (average_distance * average_distance);//5/调色盘平均距离平方，这个参数干啥的？

  const size_t rbf_color_size = rbf_lab_colors.size();

  for (size_t i = 0; i < rbf_color_size; ++i) {

    Pixel_D in_color_rbf = rbf_lab_colors[i];//一个格子的端点

    vector<Pixel_D> diffs_rbf(CENTER_NUM);

    //首元素为黑，做偏移
    for (size_t j = 0; j < CENTER_NUM; ++j) {//对调色盘上每一个值

      Pixel_D palette_diff = palette_diffs[j];

      if (isBigChange(palette_diff)) {//存在大于0.5的就算bigchange

        Pixel_D out_color_rbf;
        LAB2RGB((in_color_rbf + palette_diff), out_color_rbf);//out是rgb的
        if (isOutOfBoundary(out_color_rbf)) { //with 0.5 judge diff 是255的？
			//X+CC'超过边界时
          Pixel_D out_palette_color_lab(current_palette_lab[j] +
              palette_diff);//out_palette_color_lab
          Pixel_D diff(in_color_rbf - current_palette_lab[j]);//格子端点和该原始调色盘值的差别CX

          double ratio1 = findBoundary(out_palette_color_lab, diff, 0, 1);//Cj',C'X0,0,1,返回与RGB空间边界的交点Xb(CX即C'X0)
          double ratio2 = findBoundary(current_palette_lab[j],
                                       palette_diff, 1, 300);//C,CC',1,300,取Cb，这里300有什么意义，返回CCb/CC'

          // validation
          LAB2RGB((out_palette_color_lab + diff), out_color_rbf);
          if (!isOutOfBoundary(out_color_rbf)) {
            cerr << "Something wrong in [calculateGridResult] function --- "
                "isOutOfBoundary(out_color_rbf)!\n";
            exit(-1);
          }
          if (ratio1 > 1) {
            cerr << "Something wrong in [calculateGridResult] function --- "
                "ratio1 > 1 -- case 1!\n";
            exit(-1);
          }//因为X+CC'>boundary,所以必然小于1
          if (ratio2 < 1) {
            cerr << "Something wrong in [calculateGridResult] function --- "
                "ratio2 < 1 -- case 1!\n";
            exit(-1);
          }//C'在边界内，故必然>=1

          diffs_rbf[j] = (palette_diff - (1 - ratio1) * diff) / ratio2;//(XX0-XbX0)*CC'/CCb=XX'(XX0-XbX0=XXb)
          diffs_rbf[j].x *= ratio2;//???
        } else {
          double ratio1 = findBoundary(in_color_rbf, palette_diff, 1, 300);
          double ratio2 = findBoundary(current_palette_lab[j],
                                       palette_diff, 1, 300);
          if (ratio2 < 1) {
            cerr << "Something wrong in [calculateGridResult] function --- "
                "ratio2 < 1 -- case 2!\n";
            exit(-1);
          }
          double lambda = min(ratio1 / ratio2, 1.0);//x比C离边界远时XX'=CC'，否则按比例压缩
          diffs_rbf[j] = palette_diff * lambda;
          diffs_rbf[j].x /= lambda;//L做了啥操作？？等看完整个之后问
        }//没超边界
      } else {
        diffs_rbf[j] = palette_diff;
      }// end if(！isBigChange(palette_diff))对应调色盘小改变就直接对所有X都直接加了
    }// end iteration in CENTER_NUM，即对一个顶点X计算K个Fi

    Pixel_D transfered_rbf_grid =
        calculateSinglePoint(current_palette_lab, rbf_param,
                             in_color_rbf, diffs_rbf);

    Pixel_D rgb_double;
    LAB2RGB(transfered_rbf_grid, rgb_double);
    rbf_RGB_grids[i].x = max(0.0, min(1.0, rgb_double.x));//保证在范围内
    rbf_RGB_grids[i].y = max(0.0, min(1.0, rgb_double.y));
    rbf_RGB_grids[i].z = max(0.0, min(1.0, rgb_double.z));
  }
}

void colorTransfer(vector<Pixel_D> &pixels_rgb_d,
                   vector<Pixel_D> &current_palette_lab,
                   const vector<Pixel_D> &target_palette_lab,
                   vector<Pixel_D> &rbf_lab_colors,//隔出来的那些。。
                   vector<Pixel_D> &rbf_RGB_grids,//空的？
                   vector<vector<double>> &rbf_weight_index,
                   vector<vector<double>> &rbf_weight_map) {
  uint CENTER_NUM_WITH_BLACK = (uint) current_palette_lab.size();

  // transfer directions in LAB color space
  // between output palette color and input palette color
  vector<Pixel_D> palette_diffs(CENTER_NUM_WITH_BLACK);
  for (size_t i = 0; i < CENTER_NUM_WITH_BLACK; ++i) {
    palette_diffs[i] = target_palette_lab[i] - current_palette_lab[i];
  }//k个基向量cc'

  calculateGridResult(current_palette_lab, palette_diffs,
                      rbf_lab_colors, rbf_RGB_grids);

  const size_t image_area = pixels_rgb_d.size();
  for (size_t i = 0; i < image_area; ++i) {
    Pixel_D tmp_rgb(0, 0, 0);
    for (size_t j = 0; j < 8; ++j) {
      tmp_rgb += (rbf_RGB_grids[rbf_weight_index[i][j]] * rbf_weight_map[i][j]);
    }
    //value gamut judgement
    pixels_rgb_d[i] = tmp_rgb;
  }//每个点按原来和周围8个点的比例关系转化
}

int recolor5(string input_image_path,
             string output_image_path,
             vector<Pixel_D> desired_palette_list) {

  const uint CENTER_NUM = (uint) desired_palette_list.size();
//  size_t CENTER_NUM_WITH_BLACK = CENTER_NUM + 1;

  const uint RBF_GRID_NUM = 12;


  // read color image by BGR
  cv::Mat input_image = cv::imread(input_image_path);

  if (input_image.empty()) {
    cerr << "Error >> " << input_image_path << " cannot be loaded!!" << endl;
    return -1;
  }

  size_t image_area = (size_t) input_image.size().area();//w*h

#ifndef NDEBUG
  cout << "[recolor5] >> input image size: " << image_area << endl;
  cout << "[recolor5] >> height: " << input_image.rows << endl;
  cout << "[recolor5] >> width: " << input_image.cols << endl;
  cout << endl;
#endif

  vector<Pixel_D> pixels_rgb_d(image_area);//Pixel_D 含x,y,z

  //read b-g-r three channels from opencv Mat
  for (size_t i = 0; i < image_area; ++i) {
    //r
    pixels_rgb_d[i].x = input_image.data[3 * i + 2] / 255.0;
    //g
    pixels_rgb_d[i].y = input_image.data[3 * i + 1] / 255.0;
    //b
    pixels_rgb_d[i].z = input_image.data[3 * i] / 255.0;
  }

  vector<Pixel_D> current_palette_rgb_d(CENTER_NUM);
  vector<uint> palette_count(CENTER_NUM);

  cout << "[recolor5] >> Kmeans doing\n";
  doKmeans(pixels_rgb_d, current_palette_rgb_d, palette_count);//原调色盘，配色的时候怎么对应回来？（对所有的有效颜色，有一个对应关系，不用再对应到坐标）

  //magic data structure
  const uint rbf_color_size = (RBF_GRID_NUM + 1u) * (RBF_GRID_NUM + 1u) *
      (RBF_GRID_NUM + 1u);//+1？
  vector<Pixel_D> rbf_lab_colors(rbf_color_size);
  vector<vector<double>> rbf_weight_index(image_area, vector<double>(8));//格子的8个顶点？
  vector<vector<double>> rbf_weight_map(image_area, vector<double>(8));//小数部分？
  vector<Pixel_D> rbf_RGB_grids(rbf_color_size);


  /**
   * 原始算法中调色盘只取5色，第6色一般为黑色
   * 而作者默认黑色不适合过于调整，顾先将黑色作为一个聚类色，但最后抛弃掉
   * 此版改写中，我们不抛弃黑色，方便颜色统计。
   * 但在后续的颜色转换中，！！！！调色盘的index要做相应的位移！！！！
   */

  // start transfer
  //如果在配色方案中有黑色的存在，在应该将黑色也作为颜色迁移的一个中心点
  //目前这一版考虑颜色方案中缺少黑色的存在，故少一色进行迁移
  for (int transfer_id = 0; transfer_id < CENTER_NUM; ++transfer_id) {//这外面一层的CENTER_NUM是什么意思？一次只改一个调色盘颜色

    Pixel_D target_color_rgb_d = desired_palette_list[transfer_id];//下标transfer_id的目标调色盘颜色值

    cout << "[recolor5] >> transfer color: "
         << target_color_rgb_d.x << " - "
         << target_color_rgb_d.y << " - "
         << target_color_rgb_d.z << endl;
    cout << "[recolor5] >> transfer idx: " << transfer_id << endl;
    cout << endl;

    prepare_rbf_color(pixels_rgb_d, rbf_lab_colors,
                      rbf_weight_index, rbf_weight_map, RBF_GRID_NUM);//整个的prepare_rbf_color

    //Modifying luminance to maintain monotonicity
    vector<Pixel_D> current_palette_lab(CENTER_NUM);

    for (size_t i = 0; i < CENTER_NUM; ++i) {
      RGB2LAB(current_palette_rgb_d[i], current_palette_lab[i]);
    }//当前调色盘转lab

    cout << "[recolor5] >> modify Luminance" << endl;

    vector<Pixel_D> target_palette_rgb_d(current_palette_rgb_d);//复制

//    modifyLuminance(current_palette_lab, target_palette_rgb_d,
//                    (uint)transfer_id, target_color_rgb_d);
//    cout << "[recolor5] >> modify Luminance Done! " << endl;

    target_palette_rgb_d[transfer_id] = target_color_rgb_d;//第targetid个变

    vector<Pixel_D> target_palette_lab(CENTER_NUM);
    for (size_t i = 0; i < CENTER_NUM; ++i) {
      RGB2LAB(target_palette_rgb_d[i], target_palette_lab[i]);
    }//全转lab

#ifndef NDEBUG
    cout << "[recolor5] >> current palette:\n";

    for (size_t i = 0; i < CENTER_NUM; ++i) {
      cout << current_palette_lab[i].x << " - "
           << current_palette_lab[i].y << " - "
           << current_palette_lab[i].z << endl;
    }
    cout << endl;
    cout << "[recolor5] >> target palette:\n";
    for (size_t i = 0; i < CENTER_NUM; ++i) {
      cout << target_palette_lab[i].x << " - "
           << target_palette_lab[i].y << " - "
           << target_palette_lab[i].z << endl;
    }
    cout << endl;
#endif

    cout << "[recolor5] >> color transferring\n\n";
    colorTransfer(pixels_rgb_d, current_palette_lab, target_palette_lab,
                  rbf_lab_colors, rbf_RGB_grids,
                  rbf_weight_index, rbf_weight_map);
    cout << "[recolor5] >> color transferring DONE!\n\n";

    current_palette_rgb_d = target_palette_rgb_d;

  }
  //merge rgb channels back to image Mat
  for (size_t i = 0; i < image_area; ++i) {
    //r
    input_image.data[3 * i + 2] = (uint8_t) round(pixels_rgb_d[i].x * 255);
    //g
    input_image.data[3 * i + 1] = (uint8_t) round(pixels_rgb_d[i].y * 255);
    //b
    input_image.data[3 * i] = (uint8_t) round(pixels_rgb_d[i].z * 255);
  }

  cv::imwrite(output_image_path, input_image);

  return 0;
}

int recolor4one(string input_image_path,
                string output_image_path,
                vector<Pixel_D> current_palette_rgb_d,
                int transfer_id,
                Pixel_D target_color_rgb_d) {
  const uint CENTER_NUM = (uint) current_palette_rgb_d.size();

  const uint RBF_GRID_NUM = 12;

  // read color image by BGR
  cv::Mat input_image = cv::imread(input_image_path);

  if (input_image.empty()) {
    cerr << "Error >> " << input_image_path << " cannot be loaded!!" << endl;
    return -1;
  }

  size_t image_area = (size_t) input_image.size().area();

  vector<Pixel_D> pixels_rgb_d(image_area);

  //read b-g-r three channels from opencv Mat
  for (size_t i = 0; i < image_area; ++i) {
    //r
    pixels_rgb_d[i].x = input_image.data[3 * i + 2] / 255.0;
    //g
    pixels_rgb_d[i].y = input_image.data[3 * i + 1] / 255.0;
    //b
    pixels_rgb_d[i].z = input_image.data[3 * i] / 255.0;
  }

  //magic data structure
  const uint rbf_color_size = (RBF_GRID_NUM + 1u) * (RBF_GRID_NUM + 1u) *
      (RBF_GRID_NUM + 1u);
  vector<Pixel_D> rbf_lab_colors(rbf_color_size);
  vector<vector<double>> rbf_weight_index(image_area, vector<double>(8));
  vector<vector<double>> rbf_weight_map(image_area, vector<double>(8));
  vector<Pixel_D> rbf_RGB_grids(rbf_color_size);

  cout << "[recolor5] >> transfer color: "
       << target_color_rgb_d.x << " - "
       << target_color_rgb_d.y << " - "
       << target_color_rgb_d.z << endl;
  cout << "[recolor5] >> transfer idx: " << transfer_id << endl;
  cout << endl;

  prepare_rbf_color(pixels_rgb_d, rbf_lab_colors,
                    rbf_weight_index, rbf_weight_map, RBF_GRID_NUM);

  //Modifying luminance to maintain monotonicity
  vector<Pixel_D> current_palette_lab(CENTER_NUM);

  for (size_t i = 0; i < CENTER_NUM; ++i) {
    RGB2LAB(current_palette_rgb_d[i], current_palette_lab[i]);
  }

  cout << "[recolor5] >> modify Luminance" << endl;

  vector<Pixel_D> target_palette_rgb_d(current_palette_rgb_d);

  target_palette_rgb_d[transfer_id] = target_color_rgb_d;

  vector<Pixel_D> target_palette_lab(CENTER_NUM);
  for (size_t i = 0; i < CENTER_NUM; ++i) {
    RGB2LAB(target_palette_rgb_d[i], target_palette_lab[i]);
  }

  cout << "[recolor5] >> color transferring\n\n";
  colorTransfer(pixels_rgb_d, current_palette_lab, target_palette_lab,
                rbf_lab_colors, rbf_RGB_grids,
                rbf_weight_index, rbf_weight_map);
  cout << "[recolor5] >> color transferring DONE!\n\n";

  //merge rgb channels back to image Mat
  for (size_t i = 0; i < image_area; ++i) {
    //r
    input_image.data[3 * i + 2] = (uint8_t) round(pixels_rgb_d[i].x * 255);
    //g
    input_image.data[3 * i + 1] = (uint8_t) round(pixels_rgb_d[i].y * 255);
    //b
    input_image.data[3 * i] = (uint8_t) round(pixels_rgb_d[i].z * 255);
  }

  cv::imwrite(output_image_path, input_image);

  return 0;
}
