

#' A Reference Class for velox rasters
#'
#' @field vRaster External pointer to C++ object of class 'vRaster'
#' @field dim Raster dimensions
#' @field extent Raster extent
#' @field res Raster resolution
#' @field nbands Number of raster bands
#' @field crs Coordinate reference system (Proj4 string)
VeloxRaster <- setRefClass("VeloxRaster",
                           fields = c("vRaster",
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
#' Thus, the usual R copy on modify semantics do not apply.
#'
#' @param x A RasterLayer, RasterStack, matrix, list of matrices, or text string pointing to a GDAL-readable file.
#' @param extent An \code{extent} object or a numeric vector of length 4. Required if \code{x} is a matrix or list
#' of matrices, ignored otherwise.
#' @param res The x and y resolution of the raster as a numeric vector of length 2.
#' @param crs A character string describing a projection and datum in the PROJ.4 format.
#'
#' @return A VeloxRaster object.
velox <- function(x, extent=NULL, res=NULL, crs=NULL) {

  vras <- new(vRaster)

  if (is(x, "RasterLayer")) {

    extent = raster::as.vector(extent(x))
    origin = c(extent[1], extent[4])
    dim = c(nrow(x), ncol(x))
    res = res(x)
    if (is.na(crs(x))) {
      crs = ""
    } else {
      crs = crs(x, asText=TRUE)
    }

    vras$setRasterBands(list(raster::as.matrix(x)), origin, res, 1, crs)

    obj <- VeloxRaster$new(vRaster=vras, dim=dim, extent=extent, res=res, nbands=1, crs=crs)

    return(obj)
  }

  if (is(x, "RasterStack")) {

    extent = raster::as.vector(extent(x))
    origin = c(extent[1], extent[4])
    dim = c(nrow(x), ncol(x))
    res = res(x)
    if (is.na(crs(x))) {
      crs = ""
    } else {
      crs = crs(x, asText=TRUE)
    }
    nbands = raster::nlayers(x)
    x.ls <- vector("list", nbands)
    for (k in 1:length(x.ls)) {
      x.ls[[k]] <- raster::as.matrix(x[[k]])
    }

    vras$setRasterBands(x.ls, origin, res, nbands, crs)

    obj <- VeloxRaster$new(vRaster=vras, dim=dim, extent=extent, res=res, nbands=nbands, crs=crs)

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
    vras$setRasterBands(list(x), origin, res, 1, crs)

    obj <- VeloxRaster$new(vRaster=vras, dim=dim, extent=extent, res=res, nbands=1, crs=crs)

    return(obj)
  }

  if (is(x, "list")) {

    if (is.null(extent) | is.null(res)) {
      stop("extent and res arguments needed.")
    }
    if (!is(x[[1]], "matrix")) {
      stop("If x is a list, its elements must be of class 'matrix'.")
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
    vras$setRasterBands(x, origin, res, length(x), crs)

    obj <- VeloxRaster$new(vRaster=vras, dim=dim, extent=extent, res=res, nbands=length(x), crs=crs)

    return(obj)
  }

  if (is(x, "character")) {
    vras$read(x)
    dim <- vras$getDim()
    res <- vras$getResolution()
    extent <- vras$getExtent()
    crs <- vras$getCRS()
    nbands <- vras$getNBands()

    obj <- VeloxRaster$new(vRaster=vras, dim=dim, extent=extent, res=res, nbands=nbands, crs=crs)

    return(obj)
  }
}


#' @title Cast a VeloxRaster band as a RasterLayer object
#'
#' @description
#' \code{as.RasterLayer} creates a RasterLayer object from a VeloxRaster band.
#'
#' @param obj A VeloxRaster object.
#' @param band Integer indicating the RasterLayer band to be transformed.
#'
#' @return A RasterLayer object.
as.RasterLayer <- function(obj, band=1) {
  UseMethod("as.RasterLayer", obj)
}
as.RasterLayer.VeloxRaster <- function(obj, band=1) {

  if (band > obj$nbands) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }

  ext <- obj$extent
  crs <- obj$crs
  if (crs=="") {
    crs <- NA
  }
  ras <- raster::raster(obj$vRaster$getRasterBand(band), xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4], crs=crs)
  return(ras)
}

#' @title Cast a VeloxRaster as a RasterStack object
#'
#' @description
#' \code{as.RasterStack} creates a RasterStack object from a VeloxRaster.
#'
#' @param obj A VeloxRaster object.
#'
#' @return A RasterStack object.
as.RasterStack <- function(obj) {
  UseMethod("as.RasterStack", obj)
}
as.RasterStack.VeloxRaster <- function(obj) {
  ext <- obj$extent
  crs <- obj$crs
  if (crs=="") {
    crs <- NA
  }
  stk.ls <- vector("list", obj$nbands)
  for (k in 1:obj$nbands) {
    stk.ls[[k]] <- raster::raster(obj$vRaster$getRasterBand(k), xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4], crs=crs)
  }
  stk <- raster::stack(stk.ls)
  return(stk)
}


#' @title Write VeloxRaster to a file
#'
#' @description
#' \code{write} writes a VeloxRaster object as a GeoTIFF file.
#'
#' @details
#' Use the '.tif' file ending.
#'
#' @param obj The VeloxRaster to write.
#' @param path Output filename as character string.
#' @param overwrite Boolean indicating whether target file should be overwritten.
#'
#' @return Void.
write <- function(obj, path, overwrite=FALSE) {
  UseMethod("write", obj)
}
write.VeloxRaster <- function(obj, path, overwrite=FALSE) {
  if (overwrite & file.exists(path)) {
    file.remove(path)
  }
  if (!overwrite & file.exists(path)) {
    stop("File already exists.")
  }
  (obj$vRaster)$write(path)
}


#' @title Crop a VeloxRaster object
#'
#' @description
#' \code{crop}s a VeloxRaster object
#'
#' @details
#' Crops a VeloxRaster object to the extent of y. Note that currently the
#' extent of y must overlap with the extent of x, otherwise an error is thrown.
#'
#' @param obj A VeloxRaster object.
#' @param y An object from which an \code{extent} object can be extracted. Usually a Spatial* or Raster* object.
#' @param copy Boolean. If true, \code{crop} returns a new VeloxRaster object containing the cropped raster.
#' Otherwise the input VeloxRaster is cropped.
#'
#' @return If \code{copy=TRUE}: A VeloxRaster object.
crop <- function(obj, y, copy=FALSE) {
  UseMethod("crop", obj)
}
crop.VeloxRaster <- function(obj, y, copy=FALSE) {
  y.ext <- raster::as.vector(extent(y))
  overlaps <- (obj$vRaster)$overlaps(y.ext)
  if (!overlaps) {
    stop("No cells within extent of object y.")
  }
  if (!copy) {
    (obj$vRaster)$crop(y.ext)
    obj$dim <- (obj$vRaster)$getDim()
    obj$extent <- (obj$vRaster)$getExtent()
  } else {
    cropbox <- (obj$vRaster)$getCropBox(y.ext)
    cropext <- (obj$vRaster)$getCropExtent(cropbox)
    cls <- (obj$vRaster)$getCropBands(cropbox)
    res <- obj$res
    crs <- obj$crs
    nbands <- obj$nbands
    out <- velox(cls, cropext, res, crs)
    return(out)
  }
}


#' @title Median focal
#'
#' @description
#' Applies a median filter of dimension \code{wcol x wrow} to a VeloxRaster.
#'
#' @details
#' Padding is currently not implemented.
#'
#' @param obj A VeloxRaster object.
#' @param wrow y dimension of filter. Must be uneven integer.
#' @param wcol x dimension of filter. Must be uneven integer.
#' @param bands Numeric vector indicating bands where filter is applied.
#'
#' @return Void.
medianFocal <- function(obj, wrow, wcol, bands=1) {
  UseMethod("medianFocal", obj)
}
medianFocal.VeloxRaster <- function(obj, wrow, wcol, bands=1) {
  if (any(!(bands %in% 1:obj$nbands))) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }
  if (wrow < 0 | wcol < 0 | (wrow %% 2) == 0 | (wcol %% 2) == 0) {
    stop(paste("wrow and wcol must be positive uneven integers."))
  }
  for (i in bands) {
    (obj$vRaster)$medianFocal(wrow, wcol, i)
  }
}

#' @title Sum focal
#'
#' @description
#' Applies a focal sum with weights matrix \code{weights} to a VeloxRaster.
#'
#' @details
#' Padding is currently not implemented.
#'
#' @param obj A VeloxRaster object.
#' @param weights A numeric matrix of weights. Both dimensions must be uneven.
#' @param bands Numeric vector indicating bands where filter is applied.
#'
#' @return Void.
sumFocal <- function(obj, weights, bands=1) {
  UseMethod("sumFocal", obj)
}
sumFocal.VeloxRaster <- function(obj, weights, bands=1) {
  if (any(!(bands %in% 1:obj$nbands))) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }
  wrow <- dim(weigths)[1]
  wcol <- dim(weights)[2]
  if (wrow < 0 | wcol < 0 | (wrow %% 2) == 0 | (wcol %% 2) == 0) {
    stop(paste("The dimensions of the weights matrix must be uneven."))
  }
  for (i in bands) {
    (obj$vRaster)$sumFocal(weights, nrow(weights), ncol(weights), i)
  }
}

#' @title Mean focal
#'
#' @description
#' Applies a mean filter with weights matrix \code{weights} to a VeloxRaster.
#'
#' @details
#' Padding is currently not implemented.
#'
#' @param obj A VeloxRaster object.
#' @param weights A numeric matrix of weights. Both dimensions must be uneven.
#' @param bands Numeric vector indicating bands where filter is applied.
#'
#' @return Void.
meanFocal <- function(obj, weights, bands=1) {
  UseMethod("meanFocal", obj)
}
meanFocal.VeloxRaster <- function(obj, weights, bands=1) {
  if (any(!(bands %in% 1:obj$nbands))) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }
  wrow <- dim(weigths)[1]
  wcol <- dim(weights)[2]
  if (wrow < 0 | wcol < 0 | (wrow %% 2) == 0 | (wcol %% 2) == 0) {
    stop(paste("The dimensions of the weights matrix must be uneven."))
  }
  for (i in bands) {
    (obj$vRaster)$meanFocal(weights, nrow(weights), ncol(weights), i)
  }
}


disect <- function(sp) {
  ## sp:    SpatialPolygon* object

  ## Disect spatial polygon object
  ring.ls <- list()
  hole.ls <- list()

  polygons.ls <- sp@polygons
  for (polygons in polygons.ls) {
    Polygon.ls <- polygons@Polygons
    holeOrder <- as.numeric(unlist(strsplit(comment(polygons), " ")))
    ring.idx <- which(holeOrder==0)
    for (i in ring.idx) {

      ## Add rings
      Polygon <- Polygon.ls[[i]]
      ring.ls[[length(ring.ls)+1]] <- Polygon@coords

      ## Add nested hole list
      if (any(holeOrder==i)) {
        hole.idx <- which(holeOrder==i)
        this.hole.ls <- vector("list", length(hole.idx))
        for (j in 1:length(hole.idx)) {
          Hole <- Polygon.ls[[hole.idx[j]]]
          this.hole.ls[[j]] <- Hole@coords
        }
        hole.ls[[length(hole.ls)+1]] <- this.hole.ls
      } else {
        hole.ls[[length(hole.ls)+1]] <- integer(0)
      }
    }
  }
  return(list(ring.ls, hole.ls))
}

#' @title Extract
#'
#' @description
#' Extracts the values of all cells whose centerpoint is in \code{SpatialPolygons*} object
#' \code{sp} and applies R function \code{fun}.
#'
#' @details
#' \code{fun} must be an R function accepting a numeric vector as its sole input.
#'
#' @param obj A VeloxRaster object.
#' @param sp A SpatialPolygons* object.
#' @param fun An R function. See Details.
#'
#' @return A numeric matrix. One row per element in \code{sp}, one column per band in \code{obj}.
extract <- function(obj, sp, fun) {
  UseMethod("extract", obj)
}
extract.VeloxRaster <- function(obj, sp, fun) {
  ## obj:    VeloxRaster object
  ## sp:      SpatialPolygons* object

  out <- matrix(NA, length(sp), obj$nbands)
  for (p in 1:length(sp)) {

    ## Disect SpatialPolygons object
    sp.ls <- disect(rgeos::createSPComment(sp[p,]))
    ring.ls <- sp.ls[[1]]
    hole.ls <- sp.ls[[2]]

    ## Intersections
    hitmat.ls <- vector("list", length(ring.ls))
    for (i in 1:length(ring.ls)) {  ## For each outer ring

      ## Get this ring
      ring <- ring.ls[[i]]
      ring <- ring[-nrow(ring),]

      ## Get intersection with ring
      hitmat.ls[[i]] <- (obj$vRaster)$hittest(ring[,1], ring[,2], nrow(ring))

      ## For each hole: remove points that intersect with hole

      thishole.ls <- hole.ls[[i]]
      if (length(thishole.ls) > 0) {
        for (j in 1:length(thishole.ls)) {
          hole <- thishole.ls[[j]]
          hole <- hole[-nrow(hole),]
          hitmat.ls[[i]] <- (obj$vRaster)$unhit(hitmat.ls[[i]], hole[,1], hole[,2], nrow(hole))
        }
      }
    }
    hitmat <- do.call(rbind, hitmat.ls)
    valmat <- hitmat[,3:ncol(hitmat),drop=FALSE]
    for (k in 1:ncol(valmat)) {
      out[p,k] <- fun(valmat[,k])
    }
  }
  return(out)
}


#' @title Copy
#'
#' @description
#' Creates a deep copy of a VeloxRaster.
#'
#' @param obj A VeloxRaster object.
#'
#' @return A VeloxRaster object.
copy <- function(obj) {
  UseMethod("copy", obj)
}
copy.VeloxRaster <- function(obj) {
  ## obj:    VeloxRaster object
  out <- velox((obj$vRaster)$getRasterBands(), obj$extent, obj$res, obj$crs)
  return(out)
}


#' @title Aggregate
#'
#' @description
#' Aggregates a VeloxRaster object to a lower resolution.
#'
#' @details
#' \code{aggtype} must be one of the following: "sum", "mean", "min", "max", "median".
#'
#' @param obj A VeloxRaster object.
#' @param factor A numeric vector of length 1 or 2 indicating the aggregation factor in the x and y dimensions.
#' @param aggtype A character string indicating the aggregation type. See Details.
#'
#' @return Void.
aggregate <- function(obj, factor, aggtype) {
  UseMethod("aggregate", obj)
}
aggregate.VeloxRaster <- function(obj, factor, aggtype="sum") {
  ## obj:    VeloxRaster object
  ## factor: Aggregation factors, length 1 or 2
  ## aggtype:Aggregation type as string

  if (aggtype=="sum") {
    aggint = 0
  } else if (aggtype == "mean") {
    aggint = 1
  } else if (aggtype == "min") {
    aggint = 2
  } else if (aggtype == "max") {
    aggint = 3
  } else if (aggtype == "median") {
    aggint = 4
  } else {
    stop("Only 'sum', 'mean', 'min', 'max', 'median' are allowed in argument aggtype.")
  }

  if (length(factor) == 1) {
    factor <- c(factor, factor)
  }

  (obj$vRaster)$aggregate(factor[1:2], aggint)
  res <- (obj$vRaster)$getResolution()
  dim <- (obj$vRaster)$getDim()
  extent <- (obj$vRaster)$getExtent()
  obj$dim <- dim
  obj$res <- res
  obj$extent <- extent
}


#' @title Get coordinates
#'
#' @description
#' Returns a matrix containing the x-y coordinates of all cell center points of a VeloxRaster.
#'
#' @param obj A VeloxRaster object.
#'
#' @return A numeric matrix.
getCoordinates <- function(obj) {
  UseMethod("getCoordinates", obj)
}
getCoordinates.VeloxRaster <- function(obj) {
  cmat <- (obj$vRaster)$getCoordinates()
  return(cmat)
}
