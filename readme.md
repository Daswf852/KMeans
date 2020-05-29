# What?
The name is self explanatory, this is a simple, small and generic K-Means implementation for C++  
  
## How to use
First get opencv then build the examples:  
`mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..`  
`./bench <Big|FullHDNoise|UHDNoise>` for benchmarking  
`./palette <path to image>` for color palette extraction  
  
There are the following examples included:  
- `palette.hpp`: good practical example but uses its own data conversion things  
- `main.cpp`: color palette extraction using `palette.hpp`  
- `bench.cpp`: benchmark with big datasets  