

#' A Reference Class for velox rasters
#'
#' @field rasterbands A list of matrices containing the raster data
#' @field dim Raster dimensions
#' @field extent Raster extent
#' @field res Raster resolution
#' @field nbands Number of raster bands
#' @field crs Coordinate reference system (Proj4 string)
#'
#' @examples
#' ## Make VeloxRaster objects using the 'velox' function
#' mat <- matrix(1:100, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' class(vx)
#'
#' @export
VeloxRaster <- setRefClass("VeloxRaster",
                           fields = c("rasterbands",
                                      "dim",
                                      "extent",
                                      "res",
                                      "nbands",
                                      "crs")
)


#' @title Create a VeloxRaster object
#'
#' @description
#' \code{velox} creates a VeloxRaster object.
#'
#' @details
#' Creates a VeloxRaster object. Note that VeloxRaster objects are Reference Class objects and thus mutable.
#' Hence, the usual R copy on modify semantics do not apply.
#'
#' Note that if \code{x} is a list of VeloxRasters, the \code{extent} and \code{crs} attributes are copied
#' from the first list element.
#'
#' @param x A RasterLayer, RasterStack, matrix, list of matrices, list of VeloxRaster objects,
#'  or character string pointing to a GDAL-readable file.
#' @param extent An \code{extent} object or a numeric vector of length 4. Required if \code{x} is a matrix or list
#' of matrices, ignored otherwise.
#' @param res The x and y resolution of the raster as a numeric vector of length 2. Required if \code{x} is a matrix or list
#' of matrices, ignored otherwise.
#' @param crs Optional. A character string describing a projection and datum in the PROJ.4 format.
#' Ignored if \code{x} is a Raster* object.
#'
#' @return A VeloxRaster object.
#'
#' @examples
#' ## Create VeloxRaster from list of matrices
#' mat1 <- matrix(1:100, 10, 10)
#' mat2 <- matrix(100:1, 10, 10)
#' mat.ls <- list(mat1, mat2)
#' vx <- velox(mat.ls, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#'
#' @import Rcpp
#' @import methods
#' @import raster
#' @import rgdal
#' @useDynLib velox
#'
#' @export
velox <- function(x, extent=NULL, res=NULL, crs=NULL) {

  if (is(x, "RasterLayer")) {

    extent = as.vector(extent(x))
    dim = c(nrow(x), ncol(x))
    res = res(x)
    if (is.na(crs(x))) {
      crs = ""
    } else {
      crs = crs(x, asText=TRUE)
    }

    rasterbands <- list(as.matrix(x))

    obj <- VeloxRaster$new(rasterbands=rasterbands, dim=dim, extent=extent, res=res, nbands=1, crs=crs)

    return(obj)
  }

  if (is(x, "RasterStack")) {

    extent = as.vector(extent(x))
    origin = c(extent[1], extent[4])
    dim = c(nrow(x), ncol(x))
    res = res(x)
    if (is.na(crs(x))) {
      crs = ""
    } else {
      crs = crs(x, asText=TRUE)
    }
    nbands = raster::nlayers(x)
    rasterbands <- vector("list", nbands)
    for (k in 1:length(rasterbands)) {
      rasterbands[[k]] <- as.matrix(x[[k]])
    }

    obj <- VeloxRaster$new(rasterbands=rasterbands, dim=dim, extent=extent, res=res, nbands=nbands, crs=crs)

    return(obj)
  }

  if (is(x, "matrix")) {

    if (is.null(extent) | is.null(res)) {
      stop("extent and res arguments needed.")
    }

    origin = c(extent[1], extent[4])
    dim = c(nrow(x), ncol(x))
    if (is.null(crs)) {
      crs <- ""
    } else {
      if (is.na(crs)) {
        crs <- ""
      }
    }
    rasterbands <- list(x)

    obj <- VeloxRaster$new(rasterbands=rasterbands, dim=dim, extent=extent, res=res, nbands=1, crs=crs)

    return(obj)
  }

  if (is(x, "list")) {
    if (is(x[[1]], "matrix")) {
      ## List of matrices

      if (is.null(extent) | is.null(res)) {
        stop("extent and res arguments needed.")
      }
      if (length(unique(lapply(x,dim))) != 1) {
        stop("All list elements must have the same dimension.")
      }

      origin = c(extent[1], extent[4])
      dim = c(nrow(x[[1]]), ncol(x[[1]]))
      if (is.null(crs)) {
        crs <- ""
      } else {
        if (is.na(crs)) {
          crs <- ""
        }
      }

      rasterbands <- x

      obj <- VeloxRaster$new(rasterbands=rasterbands, dim=dim, extent=extent, res=res, nbands=length(rasterbands), crs=crs)

      return(obj)

    } else if (is(x[[1]], "VeloxRaster")) {
      ## List of VeloxRasters

      ndim <- nrow(unique(do.call("rbind", lapply(x, function(x) x$dim))))
      if (ndim > 1) {
        stop("All list elements must have the same dimension.")
      }

      dim <- x[[1]]$dim
      extent <- x[[1]]$extent
      origin <- c(extent[1], extent[4])
      res <- x[[1]]$res
      crs <- x[[1]]$crs

      rasterbands <- list()
      counter <- 0
      for (k in 1:length(x)) {
        this.nbands <- x[[k]]$nbands
        for (l in 1:this.nbands) {
          counter <- counter + 1
          rasterbands[[counter]] <- x[[k]]$rasterbands[[l]]
        }
      }
      nbands <- counter

      obj <- VeloxRaster$new(rasterbands=rasterbands, dim=dim, extent=extent, res=res, nbands=nbands, crs=crs)

      return(obj)

    } else {
      stop("If x is a list, its elements must be of class 'matrix' or 'VeloxRaster'.")
    }
  }

  if (is(x, "character")) {

    if (!file.exists(x)) {
      stop(paste("File", x, "does not exist."))
    }

    info <- suppressWarnings(GDALinfo(x))
    crs <- GDALSpatialRef(x)
    res <- info[6:7]
    nbands <- info[3]
    dim <- info[1:2]
    origin <- info[4:5]
    extent <- rep(NA, 4)
    extent[1] <- origin[1]
    extent[2] <- origin[1] + dim[2]*res[1]
    extent[3] <- origin[2]
    extent[4] <- origin[2] + dim[1]*res[2]

    rasterbands <- vector("list", nbands)
    gds <- new("GDALReadOnlyDataset", x)
    for (i in 1:nbands) {
      rasterbands[[i]] <- t(getRasterData(gds, band = i, offset = c(0, 0), as.is = TRUE, list_out=FALSE))
    }
    GDAL.close(gds)

    obj <- VeloxRaster$new(rasterbands=rasterbands, dim=dim, extent=extent, res=res, nbands=nbands, crs=crs)

    return(obj)
  }
}



#' Write a VeloxRaster to disk as a GeoTiff file
#'
#' @name VeloxRaster_write
#'
#' @param path Output filename as character string.
#' @param overwrite Boolean indicating whether target file should be overwritten.
#'
#' @return Void.
#'
#' @import rgdal
NULL
VeloxRaster$methods(write = function(path, overwrite=FALSE) {
  "See \\code{\\link{VeloxRaster_write}}."

  ## Some safe programming
  dir.path <- dirname(path)
  if (!dir.exists(dir.path)) {
    stop(paste("Directory", dir.path, "does not exists."))
  }
  if (file.exists(path) & !overwrite) {
    stop("File already exists. use overwrite=FALSE to overwrite.")
  }
  if (overwrite & file.exists(path)) {
    file.remove(path)
  }

  ## Determine data type
  chk.vec <- checktype_cpp(.self$rasterbands)
  isint <- as.logical(chk.vec[1])
  isneg <- as.logical(chk.vec[2])
  maxval <- chk.vec[3]
  if (isint) {
    for (i in 1:.self$nbands) {
      storage.mode(.self$rasterbands[[i]]) <- "integer"
    }
    if (!isneg) {
      int16 <- maxval<65534
      if (int16) {
        byte <- maxval<=255
        if (byte) {
          type <- "Byte"
        } else {
          type <- "UInt16"
        }
      } else {
        type <- "UInt32"
      }
    } else {
      int16 <- maxval<32767
      if (int16) {
        type <- "Int16"
      } else {
        type <- "Int32"
      }
    }
  } else {
    type <- "Float32"
  }

  ## Make driver object
  dr <- new("GDALDriver", "GTiff")

  ## Create GDS and write raster data
  gtd <- new("GDALTransientDataset", driver=dr, rows=.self$dim[1], cols=.self$dim[2], bands=.self$nbands, type=type)
  for (i in 1:.self$nbands) {
    putRasterData(gtd, t(.self$rasterbands[[i]]), band = i, offset = c(0, 0))
  }

  ## Set Metadata
  gt = c(.self$extent[1], .self$res[1], 0, .self$extent[4], 0, -.self$res[2])
  GDALcall(gtd, "SetGeoTransform", gt)
  GDALcall(gtd, "SetProject", .self$crs)

  ## Save and close
  saveDataset(gtd, filename=path)
  GDAL.close(gtd)
})


#' Delete a raster band from a VeloxRaster
#'
#' @name VeloxRaster_drop
#'
#' @param bands Numeric vector containing IDs of bands to be dropped.
#'
#' @return Void.
#'
#' @examples
#' ## Make VeloxRaster with 2 bands
#' mat1 <- matrix(1:100, 10, 10)
#' mat2 <- matrix(100:1, 10, 10)
#' vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
#'       crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Delete band 2
#' vx$drop(bands=2)
#'
NULL
VeloxRaster$methods(drop = function(bands) {
  "See \\code{\\link{VeloxRaster_drop}}."
  if (!(all(bands %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }

  oldbands <- 1:nbands
  keep <- oldbands[!(oldbands%in%bands)]
  if (length(keep)==0) {
    stop("Cannot drop all bands of a raster.")
  }

  nbands <<- length(keep)
  rasterbands <<- rasterbands[keep]
})
