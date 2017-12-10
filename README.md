<br>

![](docs/logo.svg)

<br>

Fast raster extraction and manipulation in R
============================================

velox is an R package for performing fast extraction and manipulation
operations on geospatial raster data. velox is fast because:

-   All raster operations are performed in C++.
-   Geometric operations are implemented with the [Boost
    Geometry](http://www.boost.org/doc/libs/1_65_1/libs/geometry/doc/html/index.html)
    libraries.
-   All data is held in memory.

velox is fully interoperable with the
[raster](https://cran.r-project.org/package=raster),
[sp](https://cran.r-project.org/package=raster), and
[sf](https://cran.r-project.org/package=raster) packages.

<br>

Features
--------

velox currently offers the following features:

-   Raster value extraction given *polygons*, *lines*, or *points*
-   Focal value calculation (i.e. moving window filters)
-   Rasterization of polygons or lines
-   Raster aggregation
-   Cropping
-   Image patch flattening and reconstruction

For more information, see the [velox project
website](https://hunzikp.github.io/velox/).

<br>

Status
------

[![Travis-CI Build
Status](https://travis-ci.org/hunzikp/velox.svg?branch=master)](https://travis-ci.org/hunzikp/velox)
[![codecov](https://codecov.io/gh/hunzikp/velox/branch/master/graph/badge.svg)](https://codecov.io/gh/hunzikp/velox)
[![CRAN
Version](http://www.r-pkg.org/badges/version/velox)](https://cran.r-project.org/package=velox)
[![develVersion](https://img.shields.io/badge/devel%20version-0.2.0.9001-green.svg?style=flat)](https://github.com/hunzikp/velox)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/last-week/velox)](https://www.r-pkg.org/pkg/velox)
