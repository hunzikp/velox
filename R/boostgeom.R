## Class defs

setOldClass("Rcpp_MultiPolygonCollection")
setOldClass("Rcpp_MultiLineCollection")
setOldClass("Rcpp_MultiPointCollection")

#' @title Rcpp pointer to MultiPolygonCollection
#'
#' @name MultiPolygonCollection
#'
#' @description
#' Rcpp pointer to MultiPolygonCollection.
#' @export MultiPolygonCollection
NULL

#' @title Rcpp pointer to MultiLineCollection
#'
#' @name MultiLineCollection
#'
#' @description
#' Rcpp pointer to MultiLineCollection.
#' @export MultiLineCollection
NULL

#' @title Rcpp pointer to MultiPointCollection
#'
#' @name MultiPointCollection
#'
#' @description
#' Rcpp pointer to MultiPointCollection.
#' @export MultiPointCollection
NULL

#' @title Rcpp pointer to BoostFactory
#'
#' @name BoostFactory
#'
#' @description
#' Rcpp pointer to BoostFactory.
#' @export BoostFactory
NULL

#' @title Rcpp pointer to PointGrid
#'
#' @name PointGrid
#'
#' @description
#' Rcpp pointer to PointGrid.
#' @export PointGrid
NULL

#' @title Rcpp pointer to BoxGrid
#'
#' @name BoxGrid
#'
#' @description
#' Rcpp pointer to BoxGrid.
#' @export BoxGrid
NULL



#' A S4 class for storing Boost objects in C++
#'
#' @description This is a virtual class for storing Rcpp pointers to C++ GeometryCollection and GridCollection objects.
#'
#' @slot geomcollection Rcpp pointer.
#' @slot crs An object of class \code{sf::crs}, storing the coordinate reference system info.
#' @slot precision A numeric scalar.
#' @export
setClass("BoostObject", contains="VIRTUAL",
         representation(geomcollection =  "ANY",
                        crs = "ANY",
                        precision = "numeric"))

#' An S4 virtual class for storing Boost geometry collections in C++
#'
#' @description This is a virtual class for storing Rcpp pointers to C++ GeometryCollection objects.
#'
#' @export
setClass("BoostGeometries", contains= c("BoostObject", 'VIRTUAL'))


#' An S4 class for storing Boost multipolygon collections in C++
#'
#' @description This is a class for storing Rcpp pointers to C++ MultiPolygonCollection objects.
#'
#' @export
setClass("BoostMultiPolygons", contains="BoostGeometries")

#' An S4 class for storing Boost multiline collections in C++
#'
#' @description This is a class for storing Rcpp pointers to C++ MultiLineCollection objects.
#'
#' @export
setClass("BoostMultiLines", contains="BoostGeometries")

#' An S4 class for storing Boost multipoint collections in C++
#'
#' @description This is a class for storing Rcpp pointers to C++ MultiPointCollection objects.
#'
#' @export
setClass("BoostMultiPoints", contains="BoostGeometries")

#' An S4 virtual class for storing Boost grids in C++
#'
#' @description This is a virtual class for storing Rcpp pointers to C++ grid objects.
#'
#' @export
setClass("BoostGrid", contains=c("VIRTUAL", "BoostObject"))

#' An S4 class for storing Boost box grids in C++
#'
#' @description This is a class for storing Rcpp pointers to C++ BoxGrid objects.
#'
#' @export
setClass("BoostBoxGrid", contains="BoostGrid")

#' An S4 class for storing Boost point grids in C++
#'
#' @description This is a class for storing Rcpp pointers to C++ PointGrid objects.
#'
#' @export
setClass("BoostPointGrid", contains="BoostGrid")


## Boost functions

#' @title Cast a sfc object as a BoostObject
#'
#' @name boost
#'
#' @description
#' \code{boost} creates a BoostObject from a sfc or VeloxRaster object.
#'
#' @param x An sfc object.
#' @param box Boolean. If \code{TRUE} and \code{x} is a \code{VeloxRaster} object, returns a BoxGrid instead of a PointGrid.
#' @param ... Currently not used.
#'
#' @return A BoostObject object.
#'
#' @examples
#' ## Make sfc_POINT
#' sfc <- sf::st_sfc(sf::st_point(c(0,1)))
#' ## Cast to BoostPoints
#' boostPoints <- boost(sfc)
#'
#' @export
boost <- function (x, ...) {
  UseMethod("boost", x)
}

#' @name boost
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

#' @name boost
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

#' @name boost
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

#' @name boost
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

#' @name boost
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

#' @name boost
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

#' @name boost
#' @method boost VeloxRaster
#' @export
boost.VeloxRaster <- function(x, box = FALSE, ...) {

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
#' @param ... Currently not used.
#'
#' @return A sfc object.
#'
#' @import sf
#'
#' @examples
#' ## Make sfc_MULTIPOINT
#' sfc <- sf::st_sfc(sf::st_multipoint(cbind(0,1)))
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

#' @name unboost
#' @method unboost BoostMultiPolygons
#' @export
unboost.BoostMultiPolygons <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  geom.ls <- boostFactory$makeMultiPolygonList(x@geomcollection)
  sfc <- st_sfc(lapply(geom.ls, st_multipolygon))
  st_crs(sfc) <- x@crs
  return(sfc)
}

#' @name unboost
#' @method unboost BoostMultiLines
#' @export
unboost.BoostMultiLines <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  geom.ls <- boostFactory$makeMultiLineList(x@geomcollection)
  sfc <- st_sfc(lapply(geom.ls, st_multilinestring))
  st_crs(sfc) <- x@crs
  return(sfc)
}

#' @name unboost
#' @method unboost BoostMultiPoints
#' @export
unboost.BoostMultiPoints <- function(x, ...) {
  boostFactory <- BoostFactory$new()
  geom.ls <- boostFactory$makeMultiPointList(x@geomcollection)
  sfc <- st_sfc(lapply(geom.ls, st_multipoint))
  st_crs(sfc) <- x@crs
  return(sfc)
}


#' @title Test whether two BoostObjects intersect
#'
#' @name bg_intersects
#' @rdname bg_intersects.generic
#'
#' @description
#' Tests whether two BoostObjects intersect (element-wise).
#'
#' @param obj1 A BoostObject.
#' @param obj2 A BoostObject.
#'
#' @return A list, with list element i an integer vector with the indices j for which intersects(x[i],y[j]) is TRUE.
#'
#' @import sf
#'
#' @examples
#' pts = boost(sf::st_sfc(sf::st_point(c(.5,.5)),
#' sf::st_point(c(1.5, 1.5)),
#' sf::st_point(c(2.5, 2.5))))
#' pol = boost(sf::st_sfc(sf::st_polygon(
#' list(rbind(c(0,0), c(2,0), c(2,2), c(0,2), c(0,0))))))
#' bg_intersects(pol, pts)
#'
#' @export
setGeneric("bg_intersects", function(obj1, obj2) {
  standardGeneric("bg_intersects")
})

## Grid intersections

#' @title Test whether two BoostObjects intersect
#'
#' @rdname bg_intersects.generic
#'
setMethod("bg_intersects", signature(obj1 = "BoostMultiPolygons", obj2 = "BoostBoxGrid"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiPolygon(obj1@geomcollection)
            return(outList)
})

#' @title Test whether two BoostObjects intersect
#'
#' @rdname bg_intersects.generic
#'
#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiLines", obj2 = "BoostBoxGrid"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiLine(obj1@geomcollection)
            return(outList)
          })

#' @title Test whether two BoostObjects intersect
#'
#' @rdname bg_intersects.generic
#'
#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiPolygons", obj2 = "BoostPointGrid"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiPolygon(obj1@geomcollection)
            return(outList)
          })


## Polygon intersections

#' @title Test whether two BoostObjects intersect
#'
#' @rdname bg_intersects.generic
#'
#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiPolygons", obj2 = "BoostGeometries"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiPolygon(obj1@geomcollection)
            return(outList)
          })

## Line intersections

#' @title Test whether two BoostObjects intersect
#'
#' @rdname bg_intersects.generic
#'
#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiLines", obj2 = "BoostGeometries"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiLine(obj1@geomcollection)
            return(outList)
          })


## Point intersections

#' @title Test whether two BoostObjects intersect
#'
#' @rdname bg_intersects.generic
#'
#' @export
setMethod("bg_intersects", signature(obj1 = "BoostMultiPoints", obj2 = "BoostGeometries"),
          function(obj1, obj2) {
            outList <- obj2@geomcollection$intersectsMultiPoint(obj1@geomcollection)
            return(outList)
          })


## Geometry subsections

#' @title Subset a BoostGeometries object
#'
#' @rdname subset.BoostGeometries
#'
#' @description
#' Extract a subset of geometries from a BoostGeometries object.
#'
#' @param x A BoostGeometries object.
#' @param i An integer vector index.
#'
#' @return A BoostGeometries object.
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

#' @title Subset a BoostGeometries object
#'
#' @rdname subset.BoostGeometries
#'
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

#' @title Subset a BoostGeometries object
#'
#' @rdname subset.BoostGeometries
#'
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

#' @title BoostGeometries Length
#'
#' @rdname length.BoostGeometries
#'
#' @description
#' Returns the length (number of Geometries) of a BoostGeometries object.
#'
#' @param x A BoostGeometries object.
#'
#' @return An integer scalar.
#' @export
setMethod("length", signature(x = "BoostGeometries"),
          function(x) {
            return(x@geomcollection$length())
          })


## Geometry plot

#' @title Plot BoostGeometries
#'
#' @rdname plot.BoostGeometries
#'
#' @description
#' Plot a BoostGeometries object using the \code{sf} plotting function.
#'
#' @param x A BoostGeometries object.
#' @param y Not used.
#' @param ... Passed to \code{sf::plot}.
#'
#' @return Void.
#' @export
setMethod("plot", signature(x = "BoostGeometries"),
          function(x, ...) {
            sfc <- unboost(x)
            plot(sfc, ...)
          })





