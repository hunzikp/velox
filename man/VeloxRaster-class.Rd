% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/velox.R
\docType{class}
\name{VeloxRaster-class}
\alias{VeloxRaster-class}
\alias{VeloxRaster}
\title{A Reference Class for velox rasters}
\description{
A Reference Class for velox rasters
}
\section{Fields}{

\describe{
\item{\code{rasterbands}}{A list of matrices containing the raster data}

\item{\code{dim}}{Raster dimensions}

\item{\code{extent}}{Raster extent}

\item{\code{res}}{Raster resolution}

\item{\code{nbands}}{Number of raster bands}

\item{\code{crs}}{Coordinate reference system (Proj4 string)}
}}

\section{Methods}{

\describe{
\item{\code{aggregate(factor, aggtype = c("sum", "mean", "min", "max", "median"))}}{See \code{\link{VeloxRaster_aggregate}}.}

\item{\code{as.matrix(band = 1)}}{See \code{\link{VeloxRaster_as.matrix}}.}

\item{\code{as.RasterBrick(assign_data_type = FALSE)}}{See \code{\link{VeloxRaster_as.RasterBrick}}.}

\item{\code{as.RasterLayer(band = 1, assign_data_type = FALSE)}}{See \code{\link{VeloxRaster_as.RasterLayer}}.}

\item{\code{as.RasterStack(assign_data_type = FALSE)}}{See \code{\link{VeloxRaster_as.RasterStack}}.}

\item{\code{col2im(
  mat,
  wrow,
  wcol,
  band,
  rowframe = 0,
  colframe = 0,
  rowstride = 1,
  colstride = 1
)}}{See \code{\link{VeloxRaster_col2im}}.}

\item{\code{crop(x)}}{See \code{\link{VeloxRaster_crop}}.}

\item{\code{drop(bands)}}{See \code{\link{VeloxRaster_drop}}.}

\item{\code{extract(sp, fun = NULL, df = FALSE, small = FALSE, legacy = FALSE)}}{See \code{\link{VeloxRaster_extract}}.}

\item{\code{extract_points(sp)}}{See \code{\link{VeloxRaster_extract_points}}.}

\item{\code{get_data_type()}}{See \code{\link{VeloxRaster_get_data_type}}.}

\item{\code{getCoordinates()}}{See \code{\link{VeloxRaster_getCoordinates}}.}

\item{\code{im2col(
  wrow,
  wcol,
  band,
  padval = 0,
  rowframe = 0,
  colframe = 0,
  rowstride = 1,
  colstride = 1
)}}{See \code{\link{VeloxRaster_im2col}}.}

\item{\code{meanFocal(weights, bands = 1)}}{See \code{\link{VeloxRaster_meanFocal}}.}

\item{\code{medianFocal(wrow, wcol, bands = 1)}}{See \code{\link{VeloxRaster_medianFocal}}.}

\item{\code{rasterize(spdf, field, band = 1, background = NULL, small = FALSE)}}{See \code{\link{VeloxRaster_rasterize}}.}

\item{\code{sumFocal(weights, bands = 1)}}{See \code{\link{VeloxRaster_sumFocal}}.}

\item{\code{write(path, overwrite = FALSE)}}{See \code{\link{VeloxRaster_write}}.}
}}

\examples{
## Make VeloxRaster objects using the 'velox' function
mat <- matrix(1:100, 10, 10)
vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
class(vx)

}
