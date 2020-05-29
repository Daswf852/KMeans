/*
Copyright © 2020 daswf852@protonmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <cmath>
#include <vector>

#ifdef _THREADED
#include <thread>
#endif

#include <opencv4/opencv2/opencv.hpp>

#include "kmeans.hpp"

namespace PaletteExtractor {
    void GetPointsFromMat(std::vector<Point::APoint_t> &dst, const cv::Mat &image) {
        cv::Mat downscaled;

        float ratio = 0.f;
        size_t target = 1280*720;

        if (image.rows * image.cols > target) {
            ratio = std::sqrt((target)/(image.rows * image.cols));
        } else {
            ratio = 1.0f;
        }

            cv::resize(image, downscaled, cv::Size(image.cols, image.rows), ratio, ratio, cv::INTER_AREA);

        dst = std::vector<Point::APoint_t>(downscaled.rows * downscaled.cols);
        for (size_t i = 0; i < downscaled.rows; i++) {
            for (size_t j = 0; j < downscaled.cols; j++) {
                cv::Vec3b color = downscaled.at<cv::Vec3b>(i, j);
                dst.at(j+(i*downscaled.cols)) = Point::NewPoint(0, {(double)color.val[2], (double)color.val[1], (double)color.val[0]});
            }
        }
    }

    std::vector<Point::APoint_t> GetPalette(unsigned int numColors, std::vector<Point::APoint_t> &points) {
        if (numColors == 0) numColors = 1;
        KMeans *km = new KMeans(numColors, points);

#ifdef _THREADED
        unsigned numThreads = std::thread::hardware_concurrency();
#endif
        if (!numThreads) numThreads = 1;
        km->SetThreadCount(numThreads);

        km->IterateUntilVariance(1.f, true);
        std::vector<Point::APoint_t> retPoints(numColors);
        for (size_t i = 0; i < numColors; i++) {
            retPoints.at(i) = Point::NewPoint(km->GetMeans().at(i));
        }
        delete km;
        return retPoints;
    }

    std::vector<Point::APoint_t> GetPalette(unsigned int numColors, std::string path) {
        cv::Mat image = cv::imread(path);
        std::vector<Point::APoint_t> dataPoints;
        PaletteExtractor::GetPointsFromMat(dataPoints, image);
        std::vector<Point::APoint_t> retVec = PaletteExtractor::GetPalette(numColors, dataPoints);
        Point::FreePoints(dataPoints);
        return retVec;
    }
}