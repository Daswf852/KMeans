#include <iostream>

#include "lib/kmeans.hpp"

struct color {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

int main() {
    std::vector<color> colors = {
        {10,20,30},
        {255,0,0},
        {0,0,0},
        {20,40,50},
    };
    std::vector<Point> points = Point::GeneratePoints<color>(colors, [](const color &col) -> Point {
        return Point({(double)col.r, (double)col.g, (double)col.b});
    });

    KMeans km(2, points);
    km.ResetMeans();
    km.RandomiseMeans();

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
}