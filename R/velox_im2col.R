#' @title im2col
#'
#' @name VeloxRaster_im2col
#'
#' @description
#' Creates a matrix of flattened image patches from a VeloxRaster band.
#' Order is left-to-right, top-to-bottom.
#' Note that if \code{any(c(rowframe, colframe)>0)}, the image patches are (partially) overlapping.
#'
#' @param wrow Patch size in the y dimension.
#' @param wcol Patch size in the x dimension.
#' @param band The band to be flattened.
#' @param padval A padding value.
#' @param rowframe A non-negative integer specifying the size of the frame around
#' the image patches in the y dimension.
#' @param colframe A non-negative integer specifying the size of the frame around
#' the image patches in the x dimension.
#' @param rowstride A positive integer denoting the stride between extracted patches
#' in the y dimension. I.e. only every \code{rowstride}th patch is extracted.
#' @param colstride A positive integer denoting the stride between extracted patches
#' in the x dimension. I.e. only every \code{colstride}th patch is extracted.
#'
#' @return A numeric matrix with \code{(wrow+2*rowframe)*(wcol+2*colframe)} columns.
NULL
VeloxRaster$methods(im2col = function(wrow, wcol, band, padval=0, rowframe=0, colframe=0,
                                      rowstride=1, colstride=1) {
  "See \\code{\\link{VeloxRaster_im2col}}."

  ## Some safety checks
  if (any(!(band %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }
  if (any(c(wrow,wcol)<1)) {
    stop(paste("wrow and wcol must be positive integers."))
  }
  if (any(c(rowframe,colframe)<0)) {
    stop(paste("rowframe and colframe must be non-negative integers."))
  }
  if (any(c(rowstride, colstride)<1)) {
    stop("rowstride and colstride must be positive integers.")
  }

  out <- im2col_cpp(rasterbands[[band]], dim,
                    wrow, wcol, band, padval, rowframe, colframe,
                    rowstride, colstride)
  return(out)
})



#' @title col2im
#'
#' @name VeloxRaster_col2im
#'
#' @description
#' Assigns values to a VeloxRaster band from a matrix of flattened image patches.
#' Patch frames, as specified by \code{rowframe} and \code{rowframe}, are not assigned.
#' This function is intended to be used with \code{mat} matrices constructed with the \code{im2col} function.
#'
#' @param mat The matrix of flattened image patches.
#' @param wrow Patch size in the y dimension.
#' @param wcol Patch size in the x dimension.
#' @param band The band to be assigned.
#' @param rowframe A non-negative integer specifying the size of the frame around
#' the image patches in the y dimension.
#' @param colframe A non-negative integer specifying the size of the frame around
#' the image patches in the x dimension.
#' @param rowstride A positive integer denoting the stride between extracted patches
#' in the y dimension.
#' @param colstride A positive integer denoting the stride between extracted patches
#' in the x dimension.
#'
#' @return Void.
NULL
VeloxRaster$methods(col2im = function(mat, wrow, wcol, band, rowframe=0, colframe=0,
                                      rowstride=1, colstride=1) {
  "See \\code{\\link{VeloxRaster_col2im}}."

  ## Some safety checks
  if (any(!(band %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }
  if (prod(dim(mat)) < prod(ceiling(dim/c(rowstride,colstride)))) {
    stop(paste("mat does not have enough entries."))
  }
  if ((wrow+2*rowframe)*(wcol+2*colframe) != ncol(mat)) {
    stop(paste("(wrow+2*rowframe)*(wcol+2*colframe) != ncol(mat)."))
  }
  if (any(c(rowstride, colstride)<1)) {
    stop("rowstride and colstride must be positive integers.")
  }

  rasterbands[[band]] <<- col2im_cpp(rasterbands[[band]], dim,
                                     mat, wrow, wcol, band, rowframe, colframe,
                                     rowstride, colstride)
})
