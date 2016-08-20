# velox

Velox is an R package for performing fast extraction and manipulation operations on geographic raster data. velox is fast because all raster computations are performed in C++ (using the excellent [Rcpp API](http://www.rcpp.org/)), and all data is held in memory. velox is intended to be used together with the [raster](https://cran.r-project.org/web/packages/raster/index.html) package, to which it provides a straightforward interface.

Currently, the following operations are implemented in velox:
+ Focal value calculation (i.e. moving window filters)
+ Raster value extraction given polygons
+ Rasterization of polygons
+ Raster aggregation
+ Cropping
+ Image patch flattening (similar to Matlab's im2col)
