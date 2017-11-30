# velox 0.2.0

Second version.

## Major changes

+ New C++ back-end based on the boost geometry library.

+ Polygon, linestring and point collections can be stored directly in C++ as BoostGeometries objects.

+ Implemented bg_intersects function for intersecting BoostGeometries objects.

+ Extraction and rasterize operations now accept sf objects.

+ Extraction and rasterize operations now accept line geometries.

+ The extract method now has a `small` argument similar to that of the raster::extract function.

+ The extract method now returns 'raw' raster values if argument `fun` is set to `NULL` (default).

+ The extract method now has a `df` boolean argument. If true, extract returns (a list of) data frames.

+ The velox function now accepts RasterBrick objects. 

+ New as.RasterBrick method.

+ Added extract_points method.


## Bug fixes

+ None



# velox 0.1.0

First version.

## Major changes

+ None

## Bug fixes

+ None
