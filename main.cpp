#include <iostream>

#include "lib/kmeans.hpp"

int main(int, char **) {
    /*std::vector<Point::APoint_t> points;
    for (int i = 0; i < 64; i++) {
        Point::APoint_t pt = Point::NewPoint(0, 8);
        points.push_back(pt);
    }
    Point::FreePoints(points);*/

    Point::APoint_t pt = Point::NewPoint(0, 8);
    std::cout<<pt.coordinates<<std::endl;
    std::cout<<pt.dimensions<<std::endl;
    Point::FreePoint(pt);
    return 0;
}