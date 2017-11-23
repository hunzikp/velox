Rcpp::loadModule("BOOSTGEOM", TRUE)

## Class defs

setOldClass("crs")
setOldClass("Rcpp_MultiPolygonCollection")
setOldClass("Rcpp_MultiLineCollection")
setOldClass("Rcpp_MultiPointCollection")

#' @export MultiPolygonCollection
#' @export MultiLineCollection
#' @export MultiPointCollection
#' @export BoostFactory
#' @export PointGrid
#' @export BoxGrid
#' @exportClass crs
NULL

#' @export
setClass("BoostGeometries",
         representation(geomcollection =  "ANY",
                        crs = "crs",
                        precision = "numeric"))
#' @export
setClass("BoostMultiPolygons", contains="BoostGeometries")

#' @export
setClass("BoostMultiLines", contains="BoostGeometries")

#' @export
setClass("BoostMultiPoints", contains="BoostGeometries")

#' @export
setClass("BoostGrid",
         representation(geomcollection =  "ANY",
                        crs = "crs",
                        precision = "numeric"))

#' @export
setClass("BoostBoxGrid", contains="BoostGrid")

#' @export
setClass("BoostPointGrid", contains="BoostGrid")


## Boost functions

#' @title Cast a sfc object as a BoostGeometries object
#'
#' @name boost
#'
#' @description
#' \code{boost} creates a BoostGeometries object from a sfc object.
#'
#' @param x An sfc object.
#'
#' @return A BoostGeometries object.
#'
#' @examples
#' ## Make sfc_POINT
#' sfc <- sf::st_sfc(st_point(c(0,1)))
#' ## Cast to BoostPoints
#' boostPoints <- boost(sfc)
#'
#' @export
boost <- function (x, ...) {
  UseMethod("boost", x)
}

#' @method boost sfc_MULTIPOLYGON
#' @export
boost.sfc_MULTIPOLYGON <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  mpc <- boostFactory$makeMultiPolygonCollection(x)
  bmp <- new("BoostMultiPolygons",
             geomcollection = mpc,
             crs = attributes(x)$crs,
             precision = attributes(x)$precision)
  return(bmp)
}

#' @method boost sfc_POLYGON
#' @export
boost.sfc_POLYGON <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  mpc <- boostFactory$makeMultiPolygonCollection(x)
  bmp <- new("BoostMultiPolygons",
             geomcollection = mpc,
             crs = attributes(x)$crs,
             precision = attributes(x)$precision)
  return(bmp)
}

#' @method boost sfc_MULTILINESTRING
#' @export
boost.sfc_MULTILINESTRING <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  mpc <- boostFactory$makeMultiLineCollection(x)
  bmp <- new("BoostMultiLines",
             geomcollection = mpc,
             crs = attributes(x)$crs,
             precision = attributes(x)$precision)
  return(bmp)
}

#' @method boost sfc_LINESTRING
#' @export
boost.sfc_LINESTRING <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  mpc <- boostFactory$makeMultiLineCollection(x)
  bmp <- new("BoostMultiLines",
             geomcollection = mpc,
             crs = attributes(x)$crs,
             precision = attributes(x)$precision)
  return(bmp)
}

#' @method boost sfc_MULTIPOINT
#' @export
boost.sfc_MULTIPOINT <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  mpc <- boostFactory$makeMultiPointCollection(x)
  bmp <- new("BoostMultiPoints",
             geomcollection = mpc,
             crs = attributes(x)$crs,
             precision = attributes(x)$precision)
  return(bmp)
}

#' @method boost sfc_POINT
#' @export
boost.sfc_POINT <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  mpc <- boostFactory$makeMultiPointCollection(x)
  bmp <- new("BoostMultiPoints",
             geomcollection = mpc,
             crs = attributes(x)$crs,
             precision = attributes(x)$precision)
  return(bmp)
}

#' @method boost VeloxRaster
#' @export
boost.VeloxRaster <- function(x, box = FALSE) {

  origin <- x$extent[c(1,4)]
  dim <- x$dim
  res <- x$res

  crsString <- ifelse(x$crs == "", NA, x$crs)
  crs <- sf::st_crs(crsString)

  boostFactory <- BoostFactory$new()

  if (box) {
    bg <- boostFactory$makeBoxGrid(origin = origin, dim = dim, res = res)
    bbg <- new("BoostBoxGrid",
              geomcollection = bg,
              crs = crs,
              precision = 0)
    return(bbg)
  } else {
    pg <- boostFactory$makePointGrid(origin = origin, dim = dim, res = res)
    bpp <- new("BoostPointGrid",
              geomcollection = pg,
              crs = crs,
              precision = 0)
    return(bpp)
  }
}

## Unboost functions

#' @title Cast a BoostGeometries object as a sfc object
#'
#' @name unboost
#'
#' @description
#' \code{unboost} creates a sfc object from a BoostGeometries object. Note that all sfc objects
#' created by unboost are of type MULTI.
#'
#' @param x A BoostGeometries object.
#'
#' @return A sfc object.
#'
#' @import sf
#'
#' @examples
#' ## Make sfc_MULTIPOINT
#' sfc <- sf::st_sfc(st_multipoint(cbind(0,1)))
#' ## Cast to BoostPoints
#' boostPoints <- boost(sfc)
#' ## Unboost
#' sfc2 <- unboost(boostPoints)
#' print(identical(sfc, sfc2))
#'
#' @export
unboost <- function (x, ...) {
  UseMethod("unboost", x)
}

#' @method unboost BoostMultiPolygons
#' @export
unboost.BoostMultiPolygons <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  geom.ls <- boostFactory$makeMultiPolygonList(x@geomcollection)
  sfc <- st_sfc(lapply(geom.ls, st_multipolygon))
  st_crs(sfc) <- x@crs
  return(sfc)
}

#' @method unboost BoostMultiLines
#' @export
unboost.BoostMultiLines <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  geom.ls <- boostFactory$makeMultiLineList(x@geomcollection)
  sfc <- st_sfc(lapply(geom.ls, st_multilinestring))
  st_crs(sfc) <- x@crs
  return(sfc)
}

#' @method unboost BoostMultiPoints
#' @export
unboost.BoostMultiPoints <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  geom.ls <- boostFactory$makeMultiPointList(x@geomcollection)
  sfc <- st_sfc(lapply(geom.ls, st_multipoint))
  st_crs(sfc) <- x@crs
  return(sfc)
}


#' @title Test whether two BoostGeometries objects intersect
#'
#' @name bg_intersects
#'
#' @description
#' Tests whether two BoostGeometries objects intersect (element-wise).
#'
#' @param obj1 A BoostGeometries object.
#' @param obj2 A BoostGeometries object;
#'
#' @return A list, with list element i an integer vector with the indices j for which intersects(x[i],y[j]) is TRUE.
#'
#' @import sf
#'
#' @examples
#' pts = boost(sf::st_sfc(st_point(c(.5,.5)), st_point(c(1.5, 1.5)), st_point(c(2.5, 2.5))))
#' pol = boost(sf::st_sfc(st_polygon(list(rbind(c(0,0), c(2,0), c(2,2), c(0,2), c(0,0))))))
#' bg_intersects(pol, pts)
#'
#' @export
setGeneric("bg_intersects", function(obj1, obj2) {
  standardGeneric("bg_intersects")
})

## Grid intersections

#' @export
setMethod("bg_intersects", signature(obj1 = "BoostBoxGrid", obj2 = "BoostMultiPolygons"),
          function(obj1, obj2) {
            outList <- obj1@geomcollection$intersectsMultiPolygon(obj2@geomcollection)
            return(outList)
})

#' @export
setMethod("bg_intersects", signature(obj1 = "BoostBoxGrid", obj2 = "BoostMultiLines"),
          function(obj1, obj2) {
            outList <- obj1@geomcollection$intersectsMultiLine(obj2@geomcollection)
            return(outList)
          })

#' @export
setMethod("bg_intersects", signature(obj1 = "BoostPointGrid", obj2 = "BoostMultiPolygons"),
          function(obj1, obj2) {
            outList <- obj1@geomcollection$intersectsMultiPolygon(obj2@geomcollection)
            return(outList)
          })


## Polygon intersections

#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiPolygons"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiPolygon(obj1@geomcollection)
            return(outList)
          })

## Line intersections

#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiLines"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiLine(obj1@geomcollection)
            return(outList)
          })


## Point intersections

#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiPoints"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiPoint(obj1@geomcollection)
            return(outList)
          })

## Geometry subsections

#' @export
setMethod("[", signature(x = "BoostMultiPolygons"),
          function(x, i) {
            new_geomcollection <- x@geomcollection$subset(i)
            bmp <- new("BoostMultiPolygons",
                       geomcollection = new_geomcollection,
                       crs = attributes(x)$crs,
                       precision = attributes(x)$precision)
            return(bmp)
          })

#' @export
setMethod("[", signature(x = "BoostMultiLines"),
          function(x, i) {
            new_geomcollection <- x@geomcollection$subset(i)
            bmp <- new("BoostMultiLines",
                       geomcollection = new_geomcollection,
                       crs = attributes(x)$crs,
                       precision = attributes(x)$precision)
            return(bmp)
          })

#' @export
setMethod("[", signature(x = "BoostMultiPoints"),
          function(x, i) {
            new_geomcollection <- x@geomcollection$subset(i)
            bmp <- new("BoostMultiPoints",
                       geomcollection = new_geomcollection,
                       crs = attributes(x)$crs,
                       precision = attributes(x)$precision)
            return(bmp)
          })


## Geometry length

#' @export
setMethod("length", signature(x = "BoostGeometries"),
          function(x) {
            return(x@geomcollection$length())
          })

## Geometry plot

#' @export
setMethod("plot", signature(x = "BoostGeometries"),
          function(x, ...) {
            sfc <- unboost(x)
            plot(sfc, ...)
          })


