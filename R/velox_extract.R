#' @title Extract Values Given Polygons
#'
#' @name VeloxRaster_extract
#'
#' @description
#' Extracts the values of all cells whose centerpoint is in \code{SpatialPolygons*} object
#' \code{sp} and optionally applies R function \code{fun}.
#'
#' @details
#' If passed, \code{fun} must be an R function accepting a numeric vector as its first (and only mandatory) argument, and returning a scalar.
#' If \code{fun} is \code{NULL}, \code{extract} returns a list of matrices, each matrix containing the raster values intersecting with the respective polygon.
#'
#' @param sp A SpatialPolygons* object.
#' @param fun An R function. See Details.
#'
#' @return
#' If \code{fun} is passed: A numeric matrix with one row per element in \code{sp}, one column per band in the VeloxRaster.
#'
#' Otherwise: A list of numeric matrices, with one list element per element in \code{sp}.
#' Each matrix consists of one column per band in the VeloxRaster, one row per raster cell intersecting with the polygon.
#'
#' @examples
#' ## Make VeloxRaster with two bands
#' set.seed(0)
#' mat1 <- matrix(rnorm(100), 10, 10)
#' mat2 <- matrix(rnorm(100), 10, 10)
#' vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
#'             crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Make SpatialPolygons
#' library(sp)
#' library(rgeos)
#' coord <- cbind(0.5, 0.5)
#' spoint <- SpatialPoints(coords=coord)
#' spols <- gBuffer(spgeom=spoint, width=0.5)
#' ## Extract
#' vx$extract(sp=spols, fun=mean)
#'
#' @import rgeos
#' @import sp
NULL
VeloxRaster$methods(extract = function(sp, fun = NULL) {
  "See \\code{\\link{VeloxRaster_extract}}."

  overlaps <- .self$overlapsExtent(sp)
  if (!overlaps) {
    out <- matrix(NA, length(sp), nbands)
    return(out)
  }

  if (!is.null(fun)) {
    out <- matrix(NA, length(sp), nbands)
  } else {
    out <- vector("list", length(sp))
  }


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

    if (!is.null(fun)) {
      ## Pass extracted values through 'fun'
      hitmat <- do.call(rbind, hitmat.ls)
      valmat <- hitmat[,3:ncol(hitmat),drop=FALSE]
      if (nrow(valmat) > 0) {
        for (k in 1:ncol(valmat)) {
          out[p,k] <- fun(valmat[,k])
        }
      }
    } else {
      ## Return 'raw' raster values for each polygon
      hitmat <- do.call(rbind, hitmat.ls)
      valmat <- hitmat[,3:ncol(hitmat),drop=FALSE]
      if (nrow(valmat) > 0) {
        out[[p]] <- valmat
      }
    }

  }
  return(out)
})

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

#' @title Extract Values Given Points
#'
#' @name VeloxRaster_extract_points
#'
#' @description
#' Given a set of points, returns all raster values of the cells with which they intersect.
#'
#' @param sp A SpatialPoints* object.
#'
#' @return A numeric matrix. One row per element in \code{sp}, one column per band in the VeloxRaster.
#'
#' @examples
#' ## Make VeloxRaster with two bands
#' set.seed(0)
#' mat1 <- matrix(rnorm(100), 10, 10)
#' mat2 <- matrix(rnorm(100), 10, 10)
#' vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
#'             crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Make SpatialPoints
#' library(sp)
#' library(rgeos)
#' coord <- cbind(runif(10), runif(10))
#' spoint <- SpatialPoints(coords=coord)
#' ## Extract
#' vx$extract_points(sp=spoint)
#'
#' @import rgeos
#' @import sp
NULL
VeloxRaster$methods(extract_points = function(sp) {
  "See \\code{\\link{VeloxRaster_extract_points}}."

  if (!inherits(sp, 'SpatialPoints')) {
    stop('sp must be of class SpatialPoints*.')
  }

  coords <- sp@coords
  out.mat <- pointextract_cpp(.self$rasterbands, .self$dim, .self$extent, .self$res, coords)

  return(out.mat)
})
