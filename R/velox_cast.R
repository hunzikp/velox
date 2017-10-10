

#' @title Cast a VeloxRaster band as a matrix
#'
#' @name VeloxRaster_as.matrix
#'
#' @description
#' \code{as.matrix} creates a matrix from a VeloxRaster band.
#'
#' @param band Integer indicating the VeloxRaster band to be transformed.
#'
#' @return A matrix.
#'
#' @examples
#' ## Make VeloxRaster
#' mat <- matrix(1:100, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Cast to matrix
#' vx.mat <- vx$as.matrix(band=1)
#' identical(mat, vx.mat)
#'
NULL
VeloxRaster$methods(as.matrix = function(band=1) {
  "See \\code{\\link{VeloxRaster_as.matrix}}."

  if (band > nbands) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }

  mat = rasterbands[[band]]
  return(mat)
})



#' @title Cast a VeloxRaster band as a RasterLayer object
#'
#' @name VeloxRaster_as.RasterLayer
#'
#' @description
#' \code{as.RasterLayer} creates a RasterLayer object from a VeloxRaster band.
#'
#' @param band Integer indicating the VeloxRaster band to be transformed.
#' @param assign_data_type Boolean indicating whether the dataType attribute of the returned RasterLayer should be set.
#' If TRUE, the dataType attribute is set to the smallest possible data type.
#'
#' @return A RasterLayer object.
#'
#' @examples
#' ## Make VeloxRaster
#' mat <- matrix(1:100, 10, 10)
#' vx <- velox(mat, extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Cast to RasterLayer
#' library(raster)
#' rl <- vx$as.RasterLayer(band=1)
#'
NULL
VeloxRaster$methods(as.RasterLayer = function(band=1, assign_data_type=FALSE) {
  "See \\code{\\link{VeloxRaster_as.RasterLayer}}."
  if (band > nbands) {
    stop(paste("VeloxRaster only has", nbands, "bands."))
  }

  ras <- raster(rasterbands[[band]], xmn=extent[1], xmx=extent[2], ymn=extent[3], ymx=extent[4], crs=crs)

  if (assign_data_type) {
    type <- .self$get_data_type()
    if (type == "Byte") {
      rtype <- 'INT1U'
    } else if (type == "UInt16") {
      rtype <- 'INT2U'
    } else if (type == 'UInt32') {
      rtype <- 'FLT4S'  # Avoid as per ?dataType
    } else if (type == 'Int16') {
      rtype <- 'INT2S'
    } else if (type == 'Int32') {
      rtype <- 'INT4S'
    } else if (type == 'Float32') {
      rtype <- 'FLT4S'
    } else if (type == 'Float64') {
      rtype <- 'FLT8S'
    }
    dataType(ras) <- rtype
  }

  return(ras)
})



#' @title Cast a VeloxRaster as a RasterStack object
#'
#' @name VeloxRaster_as.RasterStack
#'
#' @description
#' \code{as.RasterStack} creates a RasterStack object from a VeloxRaster.
#'
#' @param assign_data_type Boolean indicating whether the dataType attribute of the returned RasterStack should be set.
#' If TRUE, the dataType attribute is set to the smallest possible data type.
#'
#' @return A RasterStack object.
#'
#' @examples
#' ## Make VeloxRaster with two bands
#' mat1 <- matrix(1:100, 10, 10)
#' mat2 <- matrix(100:1, 10, 10)
#' vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
#'       crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Cast to RasterStack
#' library(raster)
#' rs <- vx$as.RasterStack()
#'
NULL
VeloxRaster$methods(as.RasterStack = function(assign_data_type=FALSE) {
  "See \\code{\\link{VeloxRaster_as.RasterStack}}."
  stk.ls <- vector("list", nbands)
  for (k in 1:nbands) {
    stk.ls[[k]] <- raster(rasterbands[[k]], xmn=extent[1], xmx=extent[2], ymn=extent[3], ymx=extent[4], crs=crs)
  }
  stk <- stack(stk.ls)

  if (assign_data_type) {
    type <- .self$get_data_type()
    if (type == "Byte") {
      rtype <- 'INT1U'
    } else if (type == "UInt16") {
      rtype <- 'INT2U'
    } else if (type == 'UInt32') {
      rtype <- 'FLT4S'  # Avoid as per ?dataType
    } else if (type == 'Int16') {
      rtype <- 'INT2S'
    } else if (type == 'Int32') {
      rtype <- 'INT4S'
    } else if (type == 'Float32') {
      rtype <- 'FLT4S'
    } else if (type == 'Float64') {
      rtype <- 'FLT8S'
    }
    dataType(stk) <- rtype
  }

  return(stk)
})


#' @title Cast a VeloxRaster as a RasterBrick object
#'
#' @name VeloxRaster_as.RasterBrick
#'
#' @description
#' \code{as.RasterBrick} creates a RasterBrick object from a VeloxRaster.
#'
#' @param assign_data_type Boolean indicating whether the dataType attribute of the returned RasterBrick should be set.
#' If TRUE, the dataType attribute is set to the smallest possible data type.
#'
#' @return A RasterBrick object.
#'
#' @examples
#' ## Make VeloxRaster with two bands
#' mat1 <- matrix(1:100, 10, 10)
#' mat2 <- matrix(100:1, 10, 10)
#' vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
#'       crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Cast to RasterBrick
#' library(raster)
#' rs <- vx$as.RasterBrick()
#'
NULL
VeloxRaster$methods(as.RasterBrick = function(assign_data_type=FALSE) {
  "See \\code{\\link{VeloxRaster_as.RasterBrick}}."
  brk.ls <- vector("list", nbands)
  for (k in 1:nbands) {
    brk.ls[[k]] <- raster(rasterbands[[k]], xmn=extent[1], xmx=extent[2], ymn=extent[3], ymx=extent[4], crs=crs)
  }
  brk <- brick(brk.ls)

  if (assign_data_type) {
    type <- .self$get_data_type()
    if (type == "Byte") {
      rtype <- 'INT1U'
    } else if (type == "UInt16") {
      rtype <- 'INT2U'
    } else if (type == 'UInt32') {
      rtype <- 'FLT4S'  # Avoid as per ?dataType
    } else if (type == 'Int16') {
      rtype <- 'INT2S'
    } else if (type == 'Int32') {
      rtype <- 'INT4S'
    } else if (type == 'Float32') {
      rtype <- 'FLT4S'
    } else if (type == 'Float64') {
      rtype <- 'FLT8S'
    }
    dataType(brk) <- rtype
  }

  return(brk)
})
