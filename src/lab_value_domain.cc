//
// author: Wu Shensen
//
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

void lab_value_domain() {
  Mat img(256, 256, CV_32FC3);

  vector<float> pv(256);

  // pixel color is in [0.f, 1.f];
  for (int t = 0; t < 256; t++){
    pv[t] = t / 255.f;
  }

  Mat img_lab;
  vector<Mat> mv;

  float min_l = FLT_MAX, max_l = FLT_MIN;
  float min_a = FLT_MAX, max_a = FLT_MIN;
  float min_b = FLT_MAX, max_b = FLT_MIN;

  for (int c1 = 0; c1 < 256; c1++){
    for (int c2 = 0; c2 < 256; c2++){
      for (int c3 = 0; c3 < 256; c3++){
        img.at<Vec3f>(c2, c3) = Vec3f(pv[c1], pv[c2], pv[c3]);
      }
    }

    cvtColor(img, img_lab, CV_BGR2Lab);

    cv::split(img_lab, mv);
    double min_val, max_val;
    cv::Point min_loc, max_loc;
    cv::minMaxLoc(mv[0], &min_val, &max_val, &min_loc, &max_loc); // L
    if (min_l > min_val) min_l = min_val;
    if (max_l < max_val) max_l = max_val;

    cv::minMaxLoc(mv[1], &min_val, &max_val, &min_loc, &max_loc); // A
    if (min_a > min_val) min_a = min_val;
    if (max_a < max_val) max_a = max_val;

    cv::minMaxLoc(mv[2], &min_val, &max_val, &min_loc, &max_loc); // B
    if (min_b > min_val) min_b = min_val;
    if (max_b < max_val) max_b = max_val;

    imshow("test", img);
    waitKey(10);
    cout << "step : " << c1 << endl;
  }

  cout << "L [" << min_l << ", " << max_l << "]" << endl;
  cout << "A [" << min_a << ", " << max_a << "]" << endl;
  cout << "B [" << min_b << ", " << max_b << "]" << endl;
}