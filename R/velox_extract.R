#' @title Extract Values Given Polygons
#'
#' @name VeloxRaster_extract
#'
#' @description
#' Extracts the values of all cells intersecting with a spatial object (line or polygon)
#' \code{sp} and optionally applies R function \code{fun}.
#'
#' @details
#' If passed, \code{fun} must be an R function accepting a numeric vector as its first (and only mandatory) argument, and returning a scalar.
#' If \code{fun} is \code{NULL}, \code{extract} returns a list of matrices, each matrix containing the raster values intersecting with the respective polygon (but see argument \code{df}).
#' If sp contains polygons, then cell-polygon intersections are calculated based on cell centroids (but see argument \code{small}).
#' If sp contains lines, then regular cell-line intersections are calculated.
#'
#' @param sp A sf* POLYGON or MULTIPOLYGON object, a sf* LINE or MULTILINE object, a SpatialPolygons* object, or a SpatialLines* object.
#' @param fun An R function. See Details.
#' @param df Boolean. If TRUE, the return value will be a data frame (or list of data frames, see Details), otherwise a matrix (or list of matrices, see Details).
#' If TRUE, a column \code{ID_sp} will be added to each data frame containing the ID of the sp object.
#' @param small Boolean. If TRUE and sp contains polygons, then raster values for small (or oddly shaped) polygons that do not intersect with any cell centroid
#' are established by intersecting the small polygon with the entire (boxed) cells.
#' @param legacy Boolean. Whether to use legacy C++ code (pre velox 0.1.0-9007).
#'
#' @return
#' If \code{fun} is passed: A numeric matrix or data frame (see argument \code{df}) with one row per element in \code{sp}, one column per band in the VeloxRaster.
#'
#' Otherwise: A list of numeric matrices or data frames (see argument \code{df}), with one list element per element in \code{sp}.
#' Each matrix/data frame consists of one column per band in the VeloxRaster, one row per raster cell intersecting with the geometry.
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
VeloxRaster$methods(extract = function(sp, fun = NULL, df=FALSE, small = FALSE, legacy = FALSE) {
  "See \\code{\\link{VeloxRaster_extract}}."

  if (!legacy) {

    # Ensure we have sfc object
    if (inherits(sp, "sf")) {
      sp <- st_geometry(sp)
    }

    # Convert to sfc / ensure argument class is correct
    isLine <- FALSE
    if (inherits(sp, "SpatialPolygons") || inherits(sp, "SpatialPolygonsDataFrame")) {
      geomc <- st_geometry(st_as_sf(sp))
      sp_IDs <- sapply(sp@polygons, function(x) slot(x, "ID"))
    } else if (inherits(sp, "SpatialLines") || inherits(sp, "SpatialLinesDataFrame")) {
      geomc <- st_geometry(st_as_sf(sp))
      sp_IDs <- sapply(sp@lines, function(x) slot(x, "ID"))
      isLine <- TRUE
    } else if (inherits(sp, "sfc_MULTIPOLYGON") || inherits(sp, "sfc_POLYGON")) {
      geomc <- sp
      sp_IDs <- 1:length(geomc)  # No polygon IDs in sf...?
    } else if (inherits(sp, "sfc_MULTILINESTRING") || inherits(sp, "sfc_LINESTRING")) {
      geomc <- sp
      sp_IDs <- 1:length(geomc)  # No polygon IDs in sf...?
      isLine <- TRUE
    } else {
      stop("Argument sp is of wrong class.")
    }

    # Prepare out
    if (!is.null(fun)) {
      out <- matrix(NA, NROW(geomc), nbands)
      rownames(out) <- sp_IDs
    } else {
      out <- vector("list", NROW(geomc))
      names(out) <- sp_IDs
    }

    # Ensure we have an overlap
    overlaps <- .self$overlapsExtent(geomc)
    if (!overlaps) {
      return(out)
    }

    # Create boost grid, boost geometries, intersect
    if (isLine) {
      boostGrid <- boost(.self, box = TRUE)  # If geomc is line, only box intersects make sense
    } else {
      boostGrid <- boost(.self)  # If geomc is line, only box intersects make sense
    }
    geomc.boost <- boost(geomc)
    intrs.ls <- bg_intersects(geomc.boost, boostGrid)

    # Iterate over polygons, bands, apply fun if not null
    missing.idx <- c()
    for (i in 1:length(intrs.ls)) {
      idx <- intrs.ls[[i]]
      if (length(idx) > 0) {
        if (is.null(fun)) {
          valmat <- matrix(NA, length(idx), nbands)
          for (band in 1:nbands) {
            valmat[,band] <- rasterbands[[band]][idx]
          }
          out[[i]] <- valmat
        } else {
          for (band in 1:nbands) {
            out[i,band] <- fun(rasterbands[[band]][idx])
          }
        }
      } else {
        missing.idx <- c(missing.idx, i)
      }
    }

    # If small activated: Fill missings using box intersect
    if (small & length(missing.idx) > 0 & !isLine) {

      # Create box grid, boost geometries, intersect
      boostBoxGrid <- boost(.self, box = TRUE)
      missing.boost <- geomc.boost[missing.idx]
      intrs.ls <- bg_intersects(missing.boost, boostBoxGrid)

      # Iterate over polygons, bands, apply fun if not null
      for (i in 1:length(intrs.ls)) {
        geom.idx <- missing.idx[i]
        idx <- intrs.ls[[i]]
        if (length(idx) > 0) {
          if (is.null(fun)) {
            valmat <- matrix(NA, length(idx), nbands)
            for (band in 1:nbands) {
              valmat[,band] <- rasterbands[[band]][idx]
            }
            out[[geom.idx]] <- valmat
          } else {
            for (band in 1:nbands) {
              out[geom.idx,band] <- fun(rasterbands[[band]][idx])
            }
          }
        }
      }
    }

  } else {
    ## Legacy code

    # Convert to sp & check argument type
    if (inherits(sp, "sfc") || inherits(sp, "sf")) {
      sp <- as(sp, "Spatial")
    }
    if (!(inherits(sp, "SpatialPolygons") || inherits(sp, "SpatialPolygonsDataFrame"))) {
      stop("In legacy mode, argument sp must be SpatialPolygons* object.")
    }

    # Prepare out
    sp_IDs <- sapply(sp@polygons, function(x) slot(x, "ID"))
    if (!is.null(fun)) {
      out <- matrix(NA, length(sp), nbands)
      rownames(out) <- sp_IDs
    } else {
      out <- vector("list", length(sp))
      names(out) <- sp_IDs
    }

    # Early return if no overlap
    overlaps <- .self$overlapsExtent(sp)
    if (!overlaps) {
      return(out)
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
  }

  # Return as df if requested
  if(df){
    if(!is.null(fun)) {
      out <- data.frame(ID_sp=sp_IDs, out)
    } else {
      nrow_each <- sapply(out, function(x) {
        nr <- nrow(x);
        if (is.null(nr)) {
          return(0)
        } else {
          return(nr)
        }
      })
      sp_IDs_rep <- unlist(Map(rep, sp_IDs, nrow_each))
      out <- data.frame(ID_sp=sp_IDs_rep, do.call("rbind", out))
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
#' @param sp A SpatialPoints* object or a sf* POINT object.
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

  if (inherits(sp, 'sf')) {
    sp <- st_geometry(sp)
  }

  if (inherits(sp, 'SpatialPoints')) {
    coords <- sp@coords
  } else if (inherits(sp, 'sfc_POINT')) {
    coords <- st_coordinates(sp)
  } else {
    stop('sp must be of class SpatialPoints* or sf* POINT.')
  }

  out.mat <- pointextract_cpp(.self$rasterbands, .self$dim, .self$extent, .self$res, coords)

  return(out.mat)
})
