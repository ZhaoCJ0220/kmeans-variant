
#include <opencv2/opencv.hpp>

using namespace cv;

int main(int argc, char **argv) {
  // Read image
  Mat im = imread( argv[1], IMREAD_GRAYSCALE );

// Set up the detector with default parameters.
  SimpleBlobDetector detector;

// Detect blobs.
  std::vector<KeyPoint> keypoints;
  detector.detect( im, keypoints);

// Draw detected blobs as red circles.
// DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
  Mat im_with_keypoints;
  drawKeypoints( im, keypoints, im_with_keypoints, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

// Show blobs
  imshow("keypoints", im_with_keypoints);
  waitKey(0);
  return 0;
}
