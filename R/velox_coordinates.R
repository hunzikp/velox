
#' @title Get coordinates
#'
#' @name VeloxRaster_getCoordinates
#'
#' @description
#' Returns a matrix containing the x-y coordinates of all cell center points of a VeloxRaster.
#'
#' @return A numeric matrix.
#'
#' @examples
#' ## Make VeloxRaster
#' mat <- matrix(1:100, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Get coordinate matrix
#' cmat <- vx$getCoordinates()
#'
NULL
VeloxRaster$methods(getCoordinates = function() {
  "See \\code{\\link{VeloxRaster_getCoordinates}}."
  cmat <- getcoordinates_cpp(dim, res, extent)
  return(cmat)
})
