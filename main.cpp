#include <iostream>
#include <random>
#include <utility>
#include <unordered_map>
#include <vector>

#include "lib/kmeans.hpp"

struct color {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

int main(int, char **) {
    //Puny test:
    /*std::vector<color> colors = {
        {10,20,30},
        {255,0,0},
        {0,0,0},
        {20,40,50},
    };
    std::vector<Point> points = Point::GeneratePoints<color>(3, colors, [](const color &col) -> Point {
        return Point({(double)col.r, (double)col.g, (double)col.b});
    });

    KMeans km(2, points);
    km.RandomiseMeans();
    km.CalculatePointsClusters();

    std::vector<int> vals = {-2, -1, 0, 1, 2};
    double avg = 0.f;
    for (size_t i = 0; i < 5; i++) {
        avg = Point::RunningAverage(avg, vals.at(i), i);
    }

    std::cout<<"Initial means: "<<std::endl;
    for (auto m : km.GetMeans()) {
        std::cout<<m<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"Points: "<<std::endl;
    for (auto p : points) {
        std::cout<<p<<std::endl;
    }
    std::cout<<std::endl;

    double target = 0.1f;
    std::cout<<"Iterating until variance target "<<target<<" is met..."<<std::endl;
    int iterations = km.IterateUntilVariance(target);
    std::cout<<"Did "<<iterations<<" iterations."<<std::endl;

    std::cout<<"New means: "<<std::endl;
    for (auto m : km.GetMeans()) {
        std::cout<<m<<std::endl;
    }

    std::cout<<"Now for some benchmarks:"<<std::endl;*/

    //Chad benchmarks:
    enum {
        FullHDNoise,
        UHDNoise,
        UHDHDRNoise
    } benchType = UHDHDRNoise;

    size_t pointCount;
    double minimum;
    double maximum;
    size_t clusters;
    double bVarianceTarget;

    int setAmountOfIterations = 50;

    switch (benchType) {
        default:
        case FullHDNoise:
            pointCount = 1920*1080;
            minimum = 0.f;
            maximum = 255.f;
            clusters = 6;
            bVarianceTarget = 3.f;
            break;
        case UHDNoise:
            pointCount = 3840*2160;
            minimum = 0.f;
            maximum = 255.f;
            clusters = 6;
            bVarianceTarget = 3.f;
            break;
        case UHDHDRNoise:
            pointCount = 3840*2160;
            minimum = 0.f;
            maximum = 1023.f;
            clusters = 6;
            bVarianceTarget = 6.f;
            break;
    }

    std::cout<<"pointCount = "<<pointCount<<std::endl;
    std::cout<<"minimum = "<<minimum<<std::endl;
    std::cout<<"maximum = "<<maximum<<std::endl;
    std::cout<<"clusters = "<<clusters<<std::endl;
    std::cout<<"bVarianceTarget = "<<bVarianceTarget<<std::endl;

    std::cout<<"Generating a dataset"<<std::endl;

    std::vector<Point> benchPoints(pointCount);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minimum, maximum);
    
    struct color {
        unsigned short r,g,b,a;
    };

    std::vector<color> colorArray(pointCount);
    for (auto it = colorArray.begin(); it != colorArray.end(); it++) {
        *it = {
            (unsigned short)dist(gen),
            (unsigned short)dist(gen),
            (unsigned short)dist(gen),
            (unsigned short)dist(gen),
        };
    }

    benchPoints = Point::GeneratePoints<color>(colorArray, [](const color &b) -> Point {
        return Point({(double)b.r, (double)b.g, (double)b.b, (double)b.a});
    });

    std::cout<<"Constructing a KMeans object"<<std::endl;
    KMeans bkm(clusters, benchPoints);
    bkm.SetThreadCount(5);
    bkm.ResetMeans();
    bkm.RandomiseMeans();

    std::cout<<bkm.GetMinimum()<<std::endl<<bkm.GetMaximum()<<std::endl;

    std::cout<<"Starting iterating"<<std::endl;

    int benchSteps = setAmountOfIterations;
    if (setAmountOfIterations) {
        for (;setAmountOfIterations;setAmountOfIterations--) {
            bkm.JustIterateParallel();
        }
    } else {
        benchSteps = bkm.IterateUntilVariance(bVarianceTarget);
    }
    std::cout<<"Did "<<benchSteps<<" iterations."<<std::endl;


    return 0;
}