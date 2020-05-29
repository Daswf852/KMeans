#include <iostream>

#include <opencv4/opencv2/opencv.hpp>

#include "kmeans.hpp"
#include "palette.hpp"

int main(int argc, char **argv) {
    std::vector<Point::APoint_t> palette = PaletteExtractor::GetPalette(6, argv[1]);

    int rectHeight = 128;
    int rectWidth = 32;
    int initialX = 0;
    int initialY = 0;

    cv::Mat drawMat(rectHeight, rectWidth * 6, CV_8UC3, cv::Scalar::all(0));


    std::string winName = "Display";
    cv::namedWindow(winName);

    for (size_t i = 0; i < palette.size(); i++) {
        const Point::APoint_t &pt = palette.at(i);
        cv::Rect rect(initialX + (rectWidth * i), initialY, rectWidth, rectHeight);
        cv::rectangle(drawMat, rect,
            cv::Scalar(
                (unsigned char)pt.coordinates[0], 
                (unsigned char)pt.coordinates[1], 
                (unsigned char)pt.coordinates[2]), cv::FILLED);
    }
    cv::imshow(winName, drawMat);
    cv::waitKey(0);

    cv::destroyAllWindows();
    return 0;
}