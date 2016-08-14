

########### REFERENCE CLASS DEFINITION
VeloxRaster <- setRefClass("VeloxRaster",
                           fields = c("vRaster",
                                         "dim",
                                         "extent",
                                         "res",
                                         "nbands",
                                         "crs")
)


########### CONSTRUCTOR
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

###### TRANSLATION
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


###### WRITE
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


###### GEOMETRIC OPERATIONS
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


##### FOCAL METHODS 
medianFocal <- function(obj, wrow, wcol, bands=1) {
  UseMethod("medianFocal", obj)
}
medianFocal.VeloxRaster <- function(obj, wrow, wcol, bands=1) {
  if (any(!(bands %in% 1:obj$nbands))) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }
  for (i in bands) {
    (obj$vRaster)$medianFocal(wrow, wcol, i)
  }
}

sumFocal <- function(obj, weights, bands=1) {
  UseMethod("sumFocal", obj)
}
sumFocal.VeloxRaster <- function(obj, weights, bands=1) {
  if (any(!(bands %in% 1:obj$nbands))) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }
  for (i in bands) {
    (obj$vRaster)$sumFocal(weights, nrow(weights), ncol(weights), i)
  }
}

meanFocal <- function(obj, weights, bands=1) {
  UseMethod("meanFocal", obj)
}
meanFocal.VeloxRaster <- function(obj, weights, bands=1) {
  if (any(!(bands %in% 1:obj$nbands))) {
    stop(paste("VeloxRaster only has", obj$nbands, "bands."))
  }
  for (i in bands) {
    (obj$vRaster)$meanFocal(weights, nrow(weights), ncol(weights), i)
  }
}


###### EXTRACT
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


###### COPY
copy <- function(obj) {
  UseMethod("copy", obj)
}
copy.VeloxRaster <- function(obj) {
  ## obj:    VeloxRaster object
  out <- velox((obj$vRaster)$getRasterBands(), obj$extent, obj$res, obj$crs)
  return(out)
}


###### AGGREGATE
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


###### COORDINATES
getCoordinates <- function(obj) {
  UseMethod("getCoordinates", obj)
}
getCoordinates.VeloxRaster <- function(obj) {
  cmat <- (obj$vRaster)$getCoordinates()
  return(cmat)
}
