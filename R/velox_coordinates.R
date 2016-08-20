
#' @title Get coordinates
#'
#' @name VeloxRaster_getCoordinates
#'
#' @description
#' Returns a matrix containing the x-y coordinates of all cell center points of a VeloxRaster.
#'
#' @return A numeric matrix.
NULL
VeloxRaster$methods(getCoordinates = function() {
  "See \\code{\\link{VeloxRaster_getCoordinates}}."
  cmat <- getcoordinates(dim, res, extent)
  return(cmat)
})
