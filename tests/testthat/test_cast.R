library(velox)
library(raster)
context("cast")

test_that("cast from RasterLayer works", {

  rs <- raster(matrix(rep(0, 9), 3, 3))
  vx <- velox(rs)
  expect_equal(vx$as.RasterLayer(), rs)

})

test_that("cast from RasterStack works", {

  rs1 <- raster(matrix(rep(0, 9), 3, 3))
  rs2 <- raster(matrix(rep(0, 9), 3, 3))
  stk <- stack(list(rs1, rs2))
  vx <- velox(stk)
  expect_equal(vx$as.RasterStack(), stk)

})

test_that("cast from RasterBrick works", {

  rs1 <- raster(matrix(rep(0, 9), 3, 3), crs="+proj=longlat +datum=WGS84 +no_defs")
  rs2 <- raster(matrix(rep(0, 9), 3, 3), crs="+proj=longlat +datum=WGS84 +no_defs")
  brk <- brick(list(rs1, rs2))
  vx <- velox(brk)
  expect_equal(vx$as.RasterBrick(), brk)
  expect_equal(vx$crs, as.character(crs(brk)))
})

test_that("cast from matrix works", {

  mat <- matrix(rep(0, 9), 3, 3)
  vx <- velox(mat, extent = c(0,1,0,1), res = 1)
  expect_equal(vx$as.matrix(), mat)

})

test_that("cast from list of matrices works", {

  mat1 <- matrix(rep(0, 9), 3, 3)
  mat2 <- matrix(rep(0, 9), 3, 3)
  stk <- stack(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent = extent(stk), res = res(stk))
  expect_equal(vx$as.RasterStack(), stk)

})

test_that("cast from list of velox rasters works", {

  mat1 <- matrix(rep(0, 9), 3, 3)
  ras1 <- raster(mat1)
  vx1 <- velox(ras1)
  mat2 <- matrix(rep(0, 9), 3, 3)
  ras2 <- raster(mat2)
  vx2 <- velox(ras2)
  stk <- stack(list(ras1, ras2))
  vxf <- velox(list(vx1, vx2))
  expect_equal(vxf$as.RasterStack(), stk)

})

test_that("cast from unkown type fails", {

  expect_error(velox(c(1,2,3)))

})
