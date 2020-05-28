# What?
The name is self explanatory, this is a simple, small and generic K-Means implementation for C++  
  
## How to use
First step is preparing data.  

```cpp
    struct SomeDataStruct {
        int foo, bar;
        double baz;
    };

    std::array<SomeDataStruct> OrigDataPoints = {...};

    auto TranslatorFunction = [](const SomeDataStruct &dataPoint) -> Point::Point_t {
        return Point::Point_t{3, {(double)foo, (double)bar, baz}};
    };

    std::vector<Point::Point_t> TranslatedDataPoints;

    TranslatedDataPoints = Point::GeneratePoints(OrigDataPoints, TranslatorFunction);
```

The second and the last step is iterating

```cpp
    KMeans km(2 /*<-- mean/cluster count*/, TranslatedDataPoints);
    
    km.IterateUntilVariance(4.f);

    for (const Point &p : km.GetMeans()) {
        std::cout<<p<<std::endl;
    }
```

`km.IterateUntilVariance(double)` will run until a targeted variance/deviation is reached (the deviation is calculated as the average of the displacement of all means on the last iteration).  
  
Alternatively you can run a set number of times:  

```cpp
    for (int i = 0; i < 8; i++) km.JustIterate();
```
  
## Note to self
Never trust `auto`
(removing it made the program 15 seconds faster from an original speed of 20 seconds with a generated data set)