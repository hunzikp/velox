#' @title Rasterize Polygons
#'
#' @name VeloxRaster_rasterize
#'
#' @description
#' Rasterizes a SpatialPolygonsDataFrame, i.e. assigns the values in the \code{field} column of the
#' SPDF to the raster cells intersecting with the respective SpatialPolygon.
#'
#' @details
#' Note that rasterization is performed sequentially. Hence, cells being contained by multiple polygons
#' are assigned the value of the last polygon in the \code{spdf} object.
#'
#'
#' @param spdf A SpatialPolygonsDataFrame object.
#' @param field A character string corresponding to the name of a numeric column in \code{spdf}.
#' @param band A positive integer denoting the ID of the band where the rasterized values are written.
#' @param background Optional. A numeric value assigned to all background cells.
#'
#' @return Void.
#'
#' @examples
#' ## Make VeloxRaster
#' mat <- matrix(0, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Make SpatialPolygonsDataFrame
#' library(sp)
#' library(rgeos)
#' coord <- cbind(0.5, 0.5)
#' spoint <- SpatialPoints(coords=coord)
#' spols <- gBuffer(spgeom=spoint, width=0.25)
#' spdf <- SpatialPolygonsDataFrame(Sr=spols, data=data.frame(value=1), match.ID=FALSE)
#' ## Rasterize, set background to -1
#' vx$rasterize(spdf=spdf, field="value", background=-1)
#'
#'@import rgeos
#'@import sp
NULL
VeloxRaster$methods(rasterize = function(spdf, field, band=1, background=NULL) {
  "See \\code{\\link{VeloxRaster_rasterize}}."

  ## Some safety checks
  if (any(!(band %in% 1:nbands))) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }
  if (!(field%in%names(spdf))) {
    stop(paste(field, " is not a column name of spdf."))
  }

  ## Fill target band with background
  if (!is.null(background)) {
    nrow <- dim[1]
    ncol <- dim[2]
    new.mat <- matrix(background, nrow, ncol)
    rasterbands[[band]] <<- new.mat
  }

  ## Iterate through polygons and color
  sp <- SpatialPolygons(spdf@polygons)
  values <- spdf@data[,field]
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

      ## Get extent of this ring
      ext <- c(min(ring[,1]), max(ring[,1]), min(ring[,2]), max(ring[,2]))

      ## Check whether ring overlaps with raster
      ring.overlaps <- .self$overlapsExtent(ext)
      if (ring.overlaps) {
        ## Crop raster
        crop.vx <- .self$copy()
        crop.vx$crop(ext)
        ## Get intersection with ring
        hitmat.ls[[i]] <- hittest_cpp(crop.vx$rasterbands, crop.vx$dim, crop.vx$extent, crop.vx$res, ring[,1], ring[,2], nrow(ring))
      } else {
        hitmat.ls[[i]] <- matrix(NA, 0, 2+nbands)
      }

      ## For each hole: remove points that intersect with hole
      thishole.ls <- hole.ls[[i]]
      if (length(thishole.ls) > 0 & ring.overlaps) {
        for (j in 1:length(thishole.ls)) {
          hole <- thishole.ls[[j]]
          hole <- hole[-nrow(hole),]
          hitmat.ls[[i]] <- unhit_cpp(hitmat.ls[[i]], hole[,1], hole[,2], nrow(hole))
        }
      }
    }

    ## Coloring
    hitmat <- do.call(rbind, hitmat.ls)
    if (nrow(hitmat) > 0) {
      value <- values[p]
      coordvalmat <- cbind(hitmat[,1:2,drop=FALSE], value)
      rasterbands[[band]] <<- color_cpp(rasterbands[[band]], coordvalmat, extent, res)
    }
  }
})
