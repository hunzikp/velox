VeloxRaster$methods(overlapsExtent = function(x) {

  if (is(x, "Extent")) {
    cext <- as.vector(x)
  } else if (is(x, "numeric")) {
    cext <- x
  } else {
    cext <- as.vector(extent(x))
  }

  overlaps <- (cext[1] < extent[2] & cext[2] > extent[1] &
                 cext[3] < extent[4] & cext[4] > extent[3])

  return(overlaps)

})


#' @title Crop a VeloxRaster object
#'
#' @name VeloxRaster_crop
#'
#' @description
#' \code{crop}s a VeloxRaster object
#'
#' @details
#' Crops a VeloxRaster object to the extent of y. Note that currently the
#' extent of y must overlap with the extent of x, otherwise an error is thrown.
#'
#' @param y An object from which an \code{extent} object can be extracted. Usually a Spatial* or Raster* object.
#'
#' @return Void.
NULL
VeloxRaster$methods(crop = function(x) {
  "See \\code{\\link{VeloxRaster_crop}}."

  overlaps <- .self$overlapsExtent(x)

  if (overlaps) {

    if (is(x, "Extent")) {
      cext <- as.vector(x)
    } else if (is(x, "numeric")) {
      cext <- x
    } else {
      cext <- as.vector(extent(x))
    }

    nrow <- dim[1]
    ncol <- dim[2]

    xres <- res[1]
    yres <- res[2]

    xmindiff <- cext[1] - extent[1]
    xmaxdiff <- extent[2] - cext[2]
    ymindiff <- cext[3] - extent[3]
    ymaxdiff <- extent[4] - cext[4]

    if (xmindiff < 0) {
      mincol <- 1
      new.xmin <- extent[1]
    } else {
      mincol <- floor(xmindiff/xres)
      new.xmin <- extent[1] + (mincol-1)*xres
    }

    if (xmaxdiff < 0) {
      maxcol <- ncol
      new.xmax <- extent[2]
    } else {
      maxcol <- ncol - floor(xmaxdiff/xres)
      new.xmax <- extent[1] + (maxcol-1)*xres
    }

    if (ymindiff < 0) {
      maxrow <- nrow
      new.ymin <- extent[3]
    } else {
      maxrow <- nrow - floor(ymindiff/yres)
      new.ymin <- extent[4] - (maxrow-1)*yres
    }

    if (ymaxdiff < 0) {
      minrow <- 1
      new.ymax <- extent[4]
    } else {
      minrow <- floor(ymaxdiff/yres)
      new.ymax <- extent[4] - (minrow-1)*yres
    }

    for (k in 1:nbands) {
      rasterbands[[k]] <<- (rasterbands[[k]])[minrow:maxrow, mincol:maxcol]
    }
    new.extent <- c(new.xmin, new.xmax, new.ymin, new.ymax)
    new.dim <- c(maxrow-minrow, maxcol-mincol)
    extent <<- new.extent
    dim <<- new.dim
  } else {
    warning("Provided object does not overlap with VeloxRaster. No cropping occurred.")
  }
})