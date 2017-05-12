library(sp)
library(rgeos)
library(velox)
library(raster)
context("Extract Points")

test_that("point extraction works", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1), crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPoints
  coord <- cbind(runif(100, 0, 2), runif(100)) # some points outside raster extent
  spoint <- SpatialPoints(coords=coord)

  ## Extract with velox
  vx_mat <- vx$extract_points(sp=spoint)

  ## Extract with ras
  ras <- vx$as.RasterStack()
  rs_mat <- extract(x = ras, y = spoint)

  ## Compare
  colnames(vx_mat) <- colnames(rs_mat)
  expect_equal(vx_mat, rs_mat)
})
