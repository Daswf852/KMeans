/*
Copyright © 2020 daswf852@protonmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#define _THREADED

#include <cassert>
#include <cfloat>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <ostream>
#include <random>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _THREADED
#include <thread>
#endif

class DifferentDimensionsException : public std::exception {
    public:
        DifferentDimensionsException() : dimsGiven(false) {}
        DifferentDimensionsException(int first, int second)
        : dimsGiven(true)
        , first(first)
        , second(second) {}

        const char *what () const throw () {
            std::stringstream wh;
            wh<<"Got 2 different dimensions";
            if (dimsGiven) wh<<": "<<first<<" and "<<second;
            return wh.str().c_str();
        }

    private:
        bool dimsGiven = false;
        int first = 0;
        int second = 0;
};

class CantDeduceDimensionsException : public std::exception {
    public:
        CantDeduceDimensionsException() {}

        const char *what () const throw () {
            return "Not enough information to get dimensions";
        }
};

class ZeroDimensionsException : public std::exception {
    public:
        ZeroDimensionsException() {}

        const char *what () const throw () {
            return "A Point_t with dimensions=0 got passed";
        }
};


namespace Point {

    struct Point_t {
        int clusterID;
        std::vector<double> coordinates;
    };

    struct APoint_t { //12+(coordinate*8) bytes per point (74MB for 1080p RGB data (obviously will do downscaling))
        uint16_t clusterID;
        uint16_t dimensions;
        double *coordinates;
    };

    APoint_t NewPoint(int cluster, std::vector<double> coordinates) {
        APoint_t retPoint;
        retPoint.clusterID = cluster;
        retPoint.dimensions = coordinates.size();
        double *pcoords = new double[coordinates.size()];
        std::copy(coordinates.begin(), coordinates.end(), pcoords);
        retPoint.coordinates = pcoords;
        return retPoint;
    }

    APoint_t NewPoint(int cluster, uint16_t dimensions) {
        APoint_t retPoint;
        retPoint.clusterID = cluster;
        retPoint.dimensions = dimensions;
        double *pcoords = new double[(size_t)dimensions];
        retPoint.coordinates = pcoords;
        return retPoint;
    }

    APoint_t NewPoint(const APoint_t &_point) {
        APoint_t retPoint;
        retPoint.clusterID = _point.clusterID;
        retPoint.dimensions = _point.dimensions;
        double *pcoords = new double[_point.dimensions];
        std::copy(_point.coordinates, _point.coordinates + _point.dimensions, pcoords);
        retPoint.coordinates = pcoords;
        return retPoint;
    }

    void FreePoint(APoint_t &point) {
        point.clusterID = -1;
        point.dimensions = 0;
        delete[] point.coordinates;
        point.coordinates = nullptr;
    }

    void CopyPoints(std::vector<APoint_t> &dst, const std::vector<APoint_t> &src) {
        dst = std::vector<APoint_t>(src.size());
        for (size_t i = 0; i < src.size(); i++) {
            dst.at(i) = NewPoint(src.at(i));
        }
    }

    void FreePoints(std::vector<APoint_t> &points) {
        for (auto it = points.begin(); it != points.end(); it++) {
            FreePoint(*it);
        }
    }

    std::ostream &operator<<(std::ostream &stream, const Point::APoint_t &point) {
        stream<<"["<<point.clusterID<<"]: (";
        bool first = true;
        for (int i = 0; i < point.dimensions; i++) {
            if (!first) stream<<", ";
            stream<<point.coordinates[i];
            first = false;
        }
        stream<<')';
        return stream;
    }

    //GetDimensions gets the dimensions of the Points in a given Point vector
    //Throws CantDeduceDimensionsException if the said vector's size is 0
    //Throws DifferentDimensionsException if any 2 Points in the vector have different dimensions
    static int GetDimensions(const std::vector<Point::APoint_t> &points) {
        if (points.size() == 0) throw CantDeduceDimensionsException();
        
        uint16_t currentDimension = points[0].dimensions;
        for (const Point::APoint_t &p : points) {
            if (currentDimension != p.dimensions) throw DifferentDimensionsException(currentDimension, p.dimensions);
        }

        return currentDimension;
    }

    //GeneratePoints generates Points from a vector of Ts
    //The description on how to generate a Point from a T is given with getPointFunction
    template <typename T, typename C>
    static void GeneratePoints(std::vector<Point::APoint_t> &dst, const C &values, std::function<Point::APoint_t(const T&)> getPointFunction) {
        if (values.size() <= 0) throw CantDeduceDimensionsException();

        dst = std::vector<Point::APoint_t>(values.size());

        for (size_t i = 0; i < values.size(); i++) {
            dst[i] = getPointFunction(values[i]);
        }

        try {
            GetDimensions(dst);
        } catch (DifferentDimensionsException ex) {
            throw ex;
        } catch (CantDeduceDimensionsException ex) {
            throw ex;
        }
    }

    //EuclideanDistance calculates the euclidean distance between 2 points
    //Throws DifferentDimensionsException if the 2 points have different dimensions
    //Far slower than ManhattanDistance
    //Decided to use this over manhattan distance after reading https://pdfs.semanticscholar.org/a630/316f9c98839098747007753a9bb6d05f752e.pdf
    ///TODO: make distance algorithm optional
    ///TODO: make user-suppliable distance functions
    ///Maybe TODO: switch to big numbers
    static double EuclideanDistance(const Point::APoint_t &p0, const Point::APoint_t &p1) {
        if (p0.dimensions != p1.dimensions) throw DifferentDimensionsException(p0.dimensions, p1.dimensions);

        double interiorSum = 0.f; //This can overflow a *lot* with far apart points
                                    //If an overflow occurs, it throws off the distance and a point's clusterID can be calculated incorrectly which leads to a std::out_of_range since it'll get assigned -1 as clusterID

        for (size_t i = 0; i < p0.dimensions; i++) {
            interiorSum += std::pow(p0.coordinates[i] - p1.coordinates[i], 2);
        }

        return std::sqrt(interiorSum);
    }

    //ManhattanDistance calculates the manhattan distance between 2 points
    //Throws DifferentDimensionsException if the 2 points have different dimensions
    //Is faster than EuclideanDistance
    //Might overflow, a lot of times over
    ///Maybe TODO: switch to big numbers
    static double ManhattanDistance(const Point::APoint_t &p0, const Point::APoint_t &p1) {
        if (p0.dimensions != p1.dimensions) throw DifferentDimensionsException(p0.dimensions, p1.dimensions);
        
        double distance = 0.f;

        for (size_t i = 0; i < p0.dimensions; i++) {
            distance += std::fabs(p0.coordinates[i] - p1.coordinates[i]);
        }

        return distance;
    }

    //RunningAverage calculates the running average of many numbers
    //prevCount contains the amount of numbers processed beforehand
    //if prevCount is 0, next is returned as-is
    static double RunningAverage(double prev, double next, unsigned int prevCount) {
        if (prevCount == 0) return next; //not-so-special case but faster than calculating
        //else if (prevCount == 1) return (prev+next)/(double)2; //profiling indicates that this made it faster too :> //a larger, controlled test indicates that it does not
        return (prev * (double)(prevCount) + next) / (double)(prevCount+1);
    }
}


//The class for calculating the cluster means of a given Point vector
//You can have as many mean counts and dimensions as you want
class KMeans {
    public:
        //Throws std::out_of_range if meanCount is not nonzero
        //Inherits the throws from Point::GetDimensions
        KMeans(unsigned int meanCount, std::vector<Point::APoint_t> &points)
            : points(points)
            , meanCount(meanCount)
            , dimensions(Point::GetDimensions(points)) {
                if (!meanCount) throw std::out_of_range("Mean count is 0");

                means = std::vector<Point::APoint_t>(meanCount);
                for (size_t i = 0; i < meanCount; i++) {
                    means[i] = Point::NewPoint(0, std::vector<double>(dimensions));
                }

                ResetMeans();
                RandomiseMeans();
                CalculatePointsClusters();
            };
        
        ~KMeans() {
            Point::FreePoints(means);
            //Freeing this->points is up to the caller
        };

        inline void ResetMeans() {
            for (std::vector<Point::APoint_t>::iterator it = means.begin(); it != means.end(); it++) {
                std::fill((*it).coordinates, (*it).coordinates + (*it).dimensions, 0.f);
            }
        }


        //IterateUntilVariance keeps iterating until a variance target is met and returns the amount of iterations done
        int IterateUntilVariance(double targetMinimumVariance, bool parallel = false) {
            double lastVariance = DBL_MAX;
            int iterations = 0;
            while (lastVariance > targetMinimumVariance) {
                lastVariance = Iterate(parallel);
                ++iterations;
            }

            return iterations;
        }

        //Iterate iterates the means and clusters and returns how much the means moved on average
        //Recommended function: IterateUntilVariance()
        double Iterate(bool parallel = false) {
            std::vector<Point::APoint_t> oldmeans;
            Point::CopyPoints(oldmeans, means);

            JustIterate(parallel);

            double averageDistance = 0.f;

            for (size_t i = 0; i < meanCount; i++) {
                averageDistance = Point::RunningAverage(averageDistance, Point::EuclideanDistance(oldmeans[i], means[i]), i);
            }

            Point::FreePoints(oldmeans);

            return averageDistance;
        }

        //JustIterate iterates the means and clusters but nothing more
        void JustIterate(bool parallel = false) {
            if (parallel) {
#ifdef _THREADED
                CalculateNewMeans();
                CalculatePointsClustersParallel();
#else
                CalculateNewMeans();
                CalculatePointsClusters();
#endif
            } else {
                CalculateNewMeans();
                CalculatePointsClusters();
            }
        }

        //GetMinimum returns the geometric minimum of the Points of a KMeans object
        Point::APoint_t GetMinimum() {
            Point::APoint_t minimum = Point::NewPoint(0, dimensions);

            for (int i = 0; i < dimensions; i++) {
                minimum.coordinates[i] = DBL_MAX;
            }

            for (const Point::APoint_t &p : points) {
                for (int i = 0; i < dimensions; i++) {
                    if (minimum.coordinates[i] > p.coordinates[i]) {
                        minimum.coordinates[i] = p.coordinates[i];
                    }
                }
            }

            return minimum;
        }

        //GetMinimum returns the geometric maximum of the Points of a KMeans object
        Point::APoint_t GetMaximum() {
            Point::APoint_t maximum = Point::NewPoint(0, dimensions);

            for (int i = 0; i < dimensions; i++) {
                maximum.coordinates[i] = DBL_MIN;
            }

            for (const Point::APoint_t &p : points) {
                for (int i = 0; i < dimensions; i++) {
                    if (maximum.coordinates[i] < p.coordinates[i]) {
                        maximum.coordinates[i] = p.coordinates[i];
                    }
                }
            }

            return maximum;
        }

        //RandomiseMeans randomises the current means
        void RandomiseMeans() {
            std::random_device rd;
            std::mt19937 gen(rd());

            Point::APoint_t minimum = GetMinimum();
            Point::APoint_t maximum = GetMaximum();

            for (size_t i = 0; i < meanCount; i++) {
                for (size_t j = 0; j < dimensions; j++) {
                    std::uniform_real_distribution<double> dist(minimum.coordinates[j], maximum.coordinates[j]);
                    means[i].coordinates[j] = dist(gen);
                    means[i].clusterID = i;
                }
            }
            
            Point::FreePoint(minimum);
            Point::FreePoint(maximum);
        }

        //GetMeans returns a const ref to the current calculated means
        const std::vector<Point::APoint_t> &GetMeans() {
            return means;
        }

        void SetThreadCount(size_t threadCount) {
            this->threadCount = threadCount;
        }

    private:
        //CalculateNewMeans calculates new means based on clusterIDs of the Points
        void CalculateNewMeans() {
            ResetMeans();
            std::vector<int> runningCounts(meanCount);
            std::fill(runningCounts.begin(), runningCounts.end(), 0);

            for (const Point::APoint_t &p : points) {
                if (p.dimensions != dimensions) throw DifferentDimensionsException(dimensions, p.dimensions);
                if (p.clusterID >= meanCount) throw std::out_of_range("p.clusterID >= meanCount: "+std::to_string(p.clusterID)+", "+std::to_string(meanCount));

                for (size_t dim = 0; dim < dimensions; dim++) {
                    means[p.clusterID].coordinates[dim] = Point::RunningAverage(means[p.clusterID].coordinates[dim], p.coordinates[dim], runningCounts[p.clusterID]);
                }
                
                ++runningCounts[p.clusterID];
            }
        }

        //CalculatePointsClusters calculates the ClusterID values of the object's Points based on the means
        void CalculatePointsClusters(bool parallel = false, size_t threadID = 0) {
            std::vector<double> distances(meanCount);
            
            size_t workload = (points.size()/threadCount);
            size_t excess = ((threadID+1 >= threadCount)?points.size()%threadCount:0);

            for (size_t i = ((parallel)?(threadID*workload):0); i < ((parallel)?workload+excess:points.size()); i++) {
                std::vector<double>::iterator distIt = distances.begin();

                for (const Point::APoint_t &m : means) {
                    *distIt = Point::EuclideanDistance(points[i], m);
                    ++distIt;
                }

                int closestClusterID = -1;
                double closestMeanDistance = DBL_MAX;
                for (size_t i = 0; i < distances.size(); i++) {
                    if (closestMeanDistance >= distances[i]) {
                        closestMeanDistance = distances[i];
                        closestClusterID = i;
                    }
                }

                points[i].clusterID = closestClusterID;
            }
        }

#ifdef _THREADED
        void CalculatePointsClustersParallel() {
            std::vector<std::thread> threads(threadCount);

            for (size_t i = 0; i < threadCount; i++) {
                threads.at(i) = std::thread(&KMeans::CalculatePointsClusters, this, true, i);
            }

            for (size_t i = 0; i < threadCount; i++) {
                threads.at(i).join();
            }
        }
#endif

        size_t threadCount = 1; //this can stay outside, we just want to filter out <thread> stuffs

        unsigned int meanCount = 1;
        std::vector<Point::APoint_t> means;

        std::vector<Point::APoint_t> &points;
        
        int dimensions = 0;
};