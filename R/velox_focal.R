#' @title Median focal
#'
#' @name VeloxRaster_medianFocal
#'
#' @description
#' Applies a median filter of dimension \code{wcol x wrow} to a VeloxRaster.
#'
#' @details
#' Padding is currently not implemented.
#'
#' @param wrow y dimension of filter. Must be uneven integer.
#' @param wcol x dimension of filter. Must be uneven integer.
#' @param bands Numeric vector indicating bands where filter is applied.
#'
#' @return Void.
NULL
VeloxRaster$methods(medianFocal = function(wrow, wcol, bands=1) {
  "See \\code{\\link{VeloxRaster_medianFocal}}."
  if (any(!(bands %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }
  if (wrow < 0 | wcol < 0 | (wrow %% 2) == 0 | (wcol %% 2) == 0) {
    stop(paste("wrow and wcol must be positive uneven integers."))
  }
  for (i in bands) {
    rasterbands[[i]] <<- medianfocal_cpp(rasterband=rasterbands[[i]], wrow=wrow, wcol=wcol, band=i)
  }
})



#' @title Sum focal
#'
#' @name VeloxRaster_sumFocal
#'
#' @description
#' Applies a focal sum with weights matrix \code{weights} to a VeloxRaster.
#'
#' @details
#' Padding is currently not implemented.
#'
#' @param weights A numeric matrix of weights. Both dimensions must be uneven.
#' @param bands Numeric vector indicating bands where filter is applied.
#'
#' @return Void.
NULL
VeloxRaster$methods(sumFocal = function(weights, bands=1) {
  "See \\code{\\link{VeloxRaster_sumFocal}}."
  if (any(!(bands %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }
  wrow <- dim(weights)[1]
  wcol <- dim(weights)[2]
  if (wrow < 0 | wcol < 0 | (wrow %% 2) == 0 | (wcol %% 2) == 0) {
    stop(paste("The dimensions of the weights matrix must be uneven."))
  }
  for (i in bands) {
    rasterbands[[i]] <<- sumfocal_cpp(rasterband=rasterbands[[i]], weights=weights, wrow=nrow(weights), wcol=ncol(weights), band=i)
  }
})



#' @title Mean focal
#'
#' @name VeloxRaster_meanFocal
#'
#' @description
#' Applies a mean filter with weights matrix \code{weights} to a VeloxRaster.
#'
#' @details
#' Padding is currently not implemented.
#'
#' @param weights A numeric matrix of weights. Both dimensions must be uneven.
#' @param bands Numeric vector indicating bands where filter is applied.
#'
#' @return Void.
NULL
VeloxRaster$methods(meanFocal = function(weights, bands=1) {
  "See \\code{\\link{VeloxRaster_meanFocal}}."
  if (any(!(bands %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }
  wrow <- dim(weights)[1]
  wcol <- dim(weights)[2]
  if (wrow < 0 | wcol < 0 | (wrow %% 2) == 0 | (wcol %% 2) == 0) {
    stop(paste("The dimensions of the weights matrix must be uneven."))
  }
  for (i in bands) {
    rasterbands[[i]] <<- meanfocal_cpp(rasterbands[[i]], weights=weights, wrow=nrow(weights), wcol=ncol(weights), band=i)
  }
})
