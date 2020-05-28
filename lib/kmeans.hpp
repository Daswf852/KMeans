/*
Copyright © 2020 daswf852@protonmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#define _THREADED

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

class Point {
    public:
        int clusterID;
        std::vector<double> coordinates;

        Point() : clusterID(0) {};

        Point(size_t dims)
            : clusterID(0)
            , coordinates(std::vector<double>(dims)) {};

        Point(std::initializer_list<double> initialiser)
            : coordinates(initialiser)
            , clusterID(0) {};

        ~Point() {};

        friend std::ostream &operator<<(std::ostream &stream, const Point &point) {
            stream<<"["<<point.clusterID<<"]: (";
            bool first = true;
            for (const double d : point.coordinates) {
                if (!first) stream<<", ";
                stream<<d;
                if (first) first = false;
            }
            stream<<')';
            return stream;
        }

        /// Some helpers:

        //GeneratePoints generates Points from a vector of Ts
        //The description on how to generate a Point from a T is given with getPointFunction
        template <typename T, typename C>
        static std::vector<Point> GeneratePoints(const C &values, std::function<Point(const T&)> getPointFunction) {
            if (values.size() <= 0) throw CantDeduceDimensionsException();

            std::vector<Point> retVec;
            for (const T &v : values) {
                Point temp = getPointFunction(v);
                /*if (temp.coordinates.size() != expectedDimensions) {
                    throw DifferentDimensionsException(expectedDimensions, temp.coordinates.size());
                }*/
                retVec.push_back(getPointFunction(v));
            };

            try {
                GetDimensions(retVec);
            } catch (DifferentDimensionsException ex) {
                throw ex;
            } catch (CantDeduceDimensionsException ex) {
                throw ex;
            }

            return retVec;
        }

        //GetDimensions gets the dimensions of the Points in a given Point vector
        //Throws CantDeduceDimensionsException if the said vector's size is 0
        //Throws DifferentDimensionsException if any 2 Points in the vector have different dimensions
        static int GetDimensions(const std::vector<Point> &points) {
            if (points.size() == 0) throw CantDeduceDimensionsException();
            
            int currentDimension = points[0].coordinates.size();
            for (const Point &p : points) {
                if (currentDimension != p.coordinates.size()) throw DifferentDimensionsException(currentDimension, p.coordinates.size());
            }

            return currentDimension;
        }

        //EuclideanDistance calculates the euclidean distance between 2 points
        //Throws DifferentDimensionsException if the 2 points have different dimensions
        //Far slower than ManhattanDistance
        //Decided to use this over manhattan distance after reading https://pdfs.semanticscholar.org/a630/316f9c98839098747007753a9bb6d05f752e.pdf
        ///TODO: make distance algorithm optional
        ///TODO: make user-suppliable distance functions
        ///Maybe TODO: switch to big numbers
        static double EuclideanDistance(const Point &p0, const Point &p1) {
            size_t p0size = p0.coordinates.size();
            size_t p1size = p1.coordinates.size();
            
            if (p0size != p1size) throw DifferentDimensionsException(p0size, p1size);

            double interiorSum = 0.f; //This can overflow a *lot* with far apart points
                                      //If an overflow occurs, it throws off the distance and a point's clusterID can be calculated incorrectly which leads to a std::out_of_range since it'll get assigned -1 as clusterID

            for (size_t i = 0; i < p0size; i++) {
                interiorSum += std::pow(p0.coordinates[i] - p1.coordinates[i], 2);
            }

            return std::sqrt(interiorSum);
        }

        //ManhattanDistance calculates the manhattan distance between 2 points
        //Throws DifferentDimensionsException if the 2 points have different dimensions
        //Is faster than EuclideanDistance
        //Might overflow, a lot of times over
        ///Maybe TODO: switch to big numbers
        static double ManhattanDistance(const Point &p0, const Point &p1) {
            size_t p0size = p0.coordinates.size();
            size_t p1size = p1.coordinates.size();
            
            if (p0size != p1size) throw DifferentDimensionsException(p0size, p1size);
            
            double distance = 0.f;

            for (size_t i = 0; i < p0size; i++) {
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
};

//The class for calculating the cluster means of a given Point vector
//You can have as many mean counts and dimensions as you want
class KMeans {
    public:
        //Throws std::out_of_range if meanCount is not nonzero
        //Inherits the throws from Point::GetDimensions
        KMeans(unsigned int meanCount, std::vector<Point> &points)
            : points(points)
            , meanCount(meanCount)
            , dimensions(Point::GetDimensions(points)) {
                if (!meanCount) throw std::out_of_range("Mean count is 0");

                means = std::vector<Point>(meanCount);
                for (size_t i = 0; i < meanCount; i++) {
                    means[i] = Point(dimensions);
                    means[i].clusterID = i;
                }

                ResetMeans();
                RandomiseMeans();
                CalculatePointsClusters();
            };
        
        ~KMeans() {};

        inline void ResetMeans() {
            for (std::vector<Point>::iterator it = means.begin(); it != means.end(); it++) {
                std::fill((*it).coordinates.begin(), (*it).coordinates.end(), 0.f);
            }
        }


        //IterateUntilVariance keeps iterating until a variance target is met and returns the amount of iterations done
        int IterateUntilVariance(double targetMinimumVariance) {
            double lastVariance = DBL_MAX;
            int iterations = 0;

            std::cout<<lastVariance<<std::endl;
            while (lastVariance > targetMinimumVariance) {
                lastVariance = Iterate();
                ++iterations;
                std::cout<<lastVariance<<std::endl;
            }

            return iterations;
        }

        //Iterate iterates the means and clusters and returns how much the means moved on average
        //Recommended function: IterateUntilVariance()
        double Iterate() {
            std::vector<Point> oldmeans = means;
            JustIterate();
            double averageDistance = 0.f;

            for (size_t i = 0; i < meanCount; i++) {
                averageDistance = Point::RunningAverage(averageDistance, Point::EuclideanDistance(oldmeans[i], means[i]), i);
            }

            return averageDistance;
        }

        //JustIterate iterates the means and clusters but nothing more
        //Useful for speed if you want to iterate for a set number of times
        void JustIterate() {
            CalculateNewMeans();
            CalculatePointsClusters();
        }

        //RandomiseMeans randomises the current means
        void RandomiseMeans() {
            std::random_device rd;
            std::mt19937 gen(rd());

            Point minimum = GetMinimum();
            Point maximum = GetMaximum();

            for (size_t i = 0; i < meanCount; i++) {
                for (size_t j = 0; j < dimensions; j++) {
                    std::uniform_real_distribution<double> dist(minimum.coordinates[j], maximum.coordinates[j]);
                    means[i].coordinates[j] = dist(gen);
                    means[i].clusterID = i;
                }
            }
        }

        //GetMinimum returns the geometric minimum of the Points of a KMeans object
        Point GetMinimum() {
            Point minimum;

            for (int i = 0; i < dimensions; i++) {
                minimum.coordinates.push_back(DBL_MAX);
            }

            for (const Point &p : points) {
                for (int i = 0; i < dimensions; i++) {
                    if (minimum.coordinates[i] > p.coordinates[i]) {
                        minimum.coordinates[i] = p.coordinates[i];
                    }
                }
            }

            return minimum;
        }

        //GetMinimum returns the geometric maximum of the Points of a KMeans object
        Point GetMaximum() {
            Point maximum;

            for (int i = 0; i < dimensions; i++) {
                maximum.coordinates.push_back(DBL_MIN);
            }

            for (const Point &p : points) {
                for (int i = 0; i < dimensions; i++) {
                    if (maximum.coordinates[i] < p.coordinates[i]) {
                        maximum.coordinates[i] = p.coordinates[i];
                    }
                }
            }

            return maximum;

        }

        //GetMeans returns a const ref to the current calculated means
        const std::vector<Point> &GetMeans() {
            return means;
        }

#ifdef _THREADED

        void SetThreadCount(size_t threadCount) {
            this->threadCount = threadCount;
        }

#endif

        void JustIterateParallel() {
            CalculateNewMeans();
            CalculatePointsClustersParallel();
        }

    private:
        //CalculateNewMeans calculates new means based on clusterIDs of the Points
        void CalculateNewMeans() {
            ResetMeans();
            std::vector<int> runningCounts(meanCount);
            std::fill(runningCounts.begin(), runningCounts.end(), 0);

            for (const Point &p : points) {
                if (p.coordinates.size() != dimensions) throw DifferentDimensionsException(dimensions, p.coordinates.size());
                if (p.clusterID >= meanCount) throw std::out_of_range("p.clusterID >= meanCount: "+std::to_string(p.clusterID)+", "+std::to_string(meanCount));

                for (size_t dim = 0; dim < dimensions; dim++) {
                    means[p.clusterID].coordinates[dim] = Point::RunningAverage(means[p.clusterID].coordinates[dim], p.coordinates[dim], runningCounts[p.clusterID]);
                }
                
                ++runningCounts[p.clusterID];
            }
        }

        //CalculatePointsClusters calculates the ClusterID values of the object's Points based on the means
        void CalculatePointsClusters() {
            std::vector<double> distances(meanCount);
            for (std::vector<Point>::iterator p = points.begin(); p != points.end(); p++) {
                std::vector<double>::iterator distIt = distances.begin();

                for (const Point &m : means) {
                    *distIt = Point::EuclideanDistance(*p, m);
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

                (*p).clusterID = closestClusterID;
            }
        }

#ifdef _THREADED

        //NYI
        void CalculateNewMeansParallel() {
        }

        //NYI
        void CalculateNewMeansParallelFunc(size_t threadID) {

        }

        void CalculatePointsClustersParallel() {
            std::vector<std::thread> threads(threadCount);

            for (size_t i = 0; i < threadCount; i++) {
                threads.at(i) = std::thread(&KMeans::CalculatePointsClustersParallelFunc, this, i);
            }

            for (size_t i = 0; i < threadCount; i++) {
                threads.at(i).join();
            }
        }

        //Now was this easy ._.
        void CalculatePointsClustersParallelFunc(size_t threadID) {
            std::vector<double> distances(meanCount);
            for (std::vector<Point>::iterator p = points.begin() + threadID; p < points.end(); p += threadCount) {
                std::vector<double>::iterator distIt = distances.begin();

                for (const Point &m : means) {
                    *distIt = Point::EuclideanDistance(*p, m);
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

                (*p).clusterID = closestClusterID;
            }
        }

        size_t threadCount = 1;

#endif

        unsigned int meanCount = 1;
        std::vector<Point> means;

        std::vector<Point> &points;
        
        int dimensions = 0;
};