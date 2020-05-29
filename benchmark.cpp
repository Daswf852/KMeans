#include <iostream>

#include "lib/kmeans.hpp"

enum EBenchType {
    FullHDNoise,
    UHDNoise,
    Big
};

struct RunType {
    size_t pointCount;
    double minimum;
    double maximum;
    size_t clusters;
    double bVarianceTarget;
    int setAmountOfIterations;
    bool MT;
};

int main(int argc, char **argv) {
    const std::unordered_map<std::string, EBenchType> presets = {
        {"FullHDNoise", FullHDNoise},
        {"UHDNoise", UHDNoise},
        {"Big", Big}
    };

    const std::unordered_map<EBenchType, RunType> presetInfo = {
        {FullHDNoise, {
            1920*1080, 0.f, 255.f,
            6,
            3.f, 0,
            false
        }},
        {UHDNoise, {
            3840*2160, 0.f, 255.f,
            6,
            3.f, 0,
            false
        }},
        {Big, {
            3840*2160, 0.f, 255.f,
            6,
            0.1f, 50,
            false
        }}
    };

    RunType runType;

    //pointCount min max clusters varTarget iterCount
    //preset
    if (argc == 2) {
        try {
            runType = presetInfo.at(presets.at(std::string(argv[1])));
        } catch (...) {
            std::cout<<"Enter a valid preset"<<std::endl;
            return 1;
        }
    } else if (argc == 7) {

    } else {
        runType = presetInfo.at(FullHDNoise);
    }

    std::cout<<"pointCount = "<<runType.pointCount<<std::endl;
    std::cout<<"minimum = "<<runType.minimum<<std::endl;
    std::cout<<"maximum = "<<runType.maximum<<std::endl;
    std::cout<<"clusters = "<<runType.clusters<<std::endl;
    std::cout<<"bVarianceTarget = "<<runType.bVarianceTarget<<std::endl;

    std::cout<<"Generating a dataset"<<std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(runType.minimum, runType.maximum);
    
    struct color {
        unsigned char r,g,b;
    };

    std::vector<color> colorArray(runType.pointCount);
    for (auto it = colorArray.begin(); it != colorArray.end(); it++) {
        *it = {
            (unsigned char)dist(gen),
            (unsigned char)dist(gen),
            (unsigned char)dist(gen),
        };
    }

    std::vector<Point::APoint_t> benchPoints;
    Point::GeneratePoints<color>(benchPoints, colorArray, [](const color &b) -> Point::APoint_t {
        return Point::NewPoint(0, {(double)b.r, (double)b.g, (double)b.b});
    });

    colorArray = std::vector<color>(1);

    std::cout<<"Constructing a KMeans object"<<std::endl;
    KMeans bkm(runType.clusters, benchPoints);
    bkm.SetThreadCount(8);
    bkm.ResetMeans();
    bkm.RandomiseMeans();

    std::cout<<bkm.GetMinimum()<<std::endl<<bkm.GetMaximum()<<std::endl;

    std::cout<<"Starting iterating"<<std::endl;

    int benchSteps = runType.setAmountOfIterations;
    if (runType.setAmountOfIterations) {
        for (;runType.setAmountOfIterations;runType.setAmountOfIterations--) {
            bkm.JustIterateParallel();
        }
    } else {
        benchSteps = bkm.IterateUntilVariance(runType.bVarianceTarget);
    }
    std::cout<<"Did "<<benchSteps<<" iterations."<<std::endl;
    for (auto m : bkm.GetMeans()) {
        std::cout<<m<<std::endl;
    }

    Point::FreePoints(benchPoints);


    return 0;
}