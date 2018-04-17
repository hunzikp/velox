#' @title Disaggregate
#'
#' @name VeloxRaster_disaggregate
#'
#' @description
#' Disaggregates a VeloxRaster object to a higher resolution.
#'
#' @param factor A numeric vector of length 1 or 2 indicating the disaggregation factor in the x and y dimensions.
#' Must be positive integers > 1.
#'
#' @return Void.
#'
#' @examples
#' ## Make VeloxRaster
#' mat <- matrix(1:100, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Disaggregate
#' vx$disaggregate(factor=c(2,2))
NULL
VeloxRaster$methods(disaggregate = function(factor) {
  if (any(factor < 1 | factor %% 1 != 0)) {
    stop("factor must be a numeric vector of positive integers >= 1")
  }

  if (length(factor) == 1) {
    factor <- c(factor, factor)
  } else if (length(factor) != 2) {
    stop("factor must be a numeric vector of length 1 or 2")
  }

  rasterbands <<- lapply(rasterbands, function(band) disaggregate_cpp(band, factor))
  dim <<- dim * factor
  res <<- res / factor
})
