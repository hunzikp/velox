library(sp)
library(rgeos)
library(velox)
library(raster)
context("Extract")

test_that("summary extract works", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)

  ## Extract Mean
  vx.emat <- vx$extract(sp=spols, fun = mean)
  rs.emat <- extract(x = brk, y = spols, fun = mean)
  colnames(vx.emat) <- colnames(rs.emat)
  rownames(vx.emat) <- rownames(rs.emat)

  ## Comparison
  expect_equal(vx.emat, rs.emat)
})

test_that("raw extract works", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)

  ## Extract raw values as lists
  vx.elist <- vx$extract(sp=spols, fun = NULL)
  rs.elist <- extract(x = brk, y = spols, fun = NULL)

  ## Comparison
  expect_true(all(vx.elist[[1]][,1] %in% rs.elist[[1]][,1]))
  expect_true(all(vx.elist[[1]][,2] %in% rs.elist[[1]][,2]))
  expect_true(all(vx.elist[[2]][,1] %in% rs.elist[[2]][,1]))
  expect_true(all(vx.elist[[2]][,2] %in% rs.elist[[2]][,2]))

})

test_that("both extract work with df=TRUE", {
  
  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")
  
  ## Make SpatialPolygons
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)
  
  ## Extract raw values as data-frame
  vx.raw.df <- vx$extract(sp=spols, fun = NULL, df=TRUE)
  rs.raw.df <- extract(x = brk, y = spols, fun = NULL, df=TRUE)
  dimnames(rs.raw.df) <- dimnames(vx.raw.df)

  ## Extract aggregated values as data-frame
  vx.foo.df <- vx$extract(sp=spols, fun = mean, df=TRUE)
  rs.foo.df <- extract(x = brk, y = spols, fun = mean, df=TRUE)
  dimnames(rs.foo.df) <- dimnames(vx.foo.df)
  
  ## Comparison
  expect_true(all(vx.raw.df == rs.raw.df))
  expect_true(all(vx.foo.df == rs.foo.df))
})

test_that("raw extract works if polygon too small to intersect", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.0000001, byid = TRUE)

  ## Extract raw values as lists
  vx.elist <- vx$extract(sp=spols, fun = NULL)

  ## Comparison
  expect_equal(vx.elist[[1]], NULL)
  expect_equal(vx.elist[[2]], NULL)
})
