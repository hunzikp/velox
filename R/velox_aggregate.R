#' @title Aggregate
#'
#' @name VeloxRaster_aggregate
#'
#' @description
#' Aggregates a VeloxRaster object to a lower resolution.
#'
#' @details
#' \code{aggtype} must be one of the following: "sum", "mean", "min", "max", "median".
#'
#' @param factor A numeric vector of length 1 or 2 indicating the aggregation factor in the x and y dimensions.
#' Must be positive integers > 1.
#' @param aggtype A character string indicating the aggregation type. See Details.
#'
#' @return Void.
#'
#' @examples
#' ## Make VeloxRaster
#' mat <- matrix(1:100, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Aggregate
#' vx$aggregate(factor=c(2,2), aggtype='sum')
NULL
VeloxRaster$methods(aggregate = function(factor, aggtype = c("sum", "mean", "min", "max", "median")) {
  "See \\code{\\link{VeloxRaster_aggregate}}."

  if (any(factor <= 1 | factor%%1!=0)) {
    stop("factor must be a numeric vector of positive integers > 1.")
  }
  aggint <- switch(match.arg(aggtype), "sum"=0, "mean"=1, "min"=2, "max"=3, "median"=4)

  if (length(factor) == 1) {
    factor <- c(factor, factor)
  }

  agg.ls <- aggregate_cpp(rasterbands, dim, res, factor, aggint)
  rasterbands <<- agg.ls[[1]]
  dim <<- agg.ls[[2]]
  res <<- agg.ls[[3]]

  xmax <- extent[1] + dim[2]*res[1]
  ymin <- extent[4] - dim[1]*res[2]
  extent[2] <<- xmax
  extent[3] <<- ymin
})
