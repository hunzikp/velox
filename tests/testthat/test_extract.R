library(sp)
library(sf)
library(rgeos)
library(velox)
library(raster)
context("Extract")

test_that("summary extract works with lines and polygons", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons / sfc_MPG
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)
  pols.sf <- st_as_sf(spols)

  ## Extract Mean using sp
  vx.emat <- vx$extract(sp=spols, fun = mean)
  rs.emat <- extract(x = brk, y = spols, fun = mean)
  colnames(vx.emat) <- colnames(rs.emat)
  rownames(vx.emat) <- rownames(rs.emat)

  ## Comparison using sp
  expect_equal(vx.emat, rs.emat)

  ## Extract Mean using sf
  vx.emat <- vx$extract(sp=pols.sf, fun = mean)
  rs.emat <- extract(x = brk, y = spols, fun = mean)
  colnames(vx.emat) <- colnames(rs.emat)
  rownames(vx.emat) <- rownames(rs.emat)

  ## Comparison using sf
  expect_equal(vx.emat, rs.emat)

  ## Make SpatialLines / sfc_MLS
  lines.sf <- st_sfc(st_linestring(matrix(c(0.1,0,
                                     1,1), nrow = 2, ncol = 2, byrow=TRUE)))
  lines.sl <- as(lines.sf, 'Spatial')

  ## Extract mean using sp
  vx.emat <- vx$extract(sp=lines.sl, fun = mean)
  rs.emat <- extract(x = brk, y = lines.sl, fun = mean)
  colnames(vx.emat) <- colnames(rs.emat)
  rownames(vx.emat) <- rownames(rs.emat)

  ## Comparison using sp
  expect_equal(vx.emat, rs.emat)

  ## Extract mean using sf
  vx.emat <- vx$extract(sp=lines.sf, fun = mean)
  rs.emat <- extract(x = brk, y = lines.sl, fun = mean)
  colnames(vx.emat) <- colnames(rs.emat)
  rownames(vx.emat) <- rownames(rs.emat)

  ## Comparison using sf
  expect_equal(vx.emat, rs.emat)

})

test_that("raw extract works with lines and polygons", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons / sfc_MPG
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)
  pols.sf <- st_as_sf(spols)

  ## Extract raw values as lists using sp
  vx.elist <- vx$extract(sp=spols, fun = NULL)
  rs.elist <- extract(x = brk, y = spols, fun = NULL)

  ## Comparison using sp
  expect_true(all(vx.elist[[1]][,1] %in% rs.elist[[1]][,1]))
  expect_true(all(vx.elist[[1]][,2] %in% rs.elist[[1]][,2]))
  expect_true(all(vx.elist[[2]][,1] %in% rs.elist[[2]][,1]))
  expect_true(all(vx.elist[[2]][,2] %in% rs.elist[[2]][,2]))

  ## Extract raw values as lists using sf
  vx.elist <- vx$extract(sp=pols.sf, fun = NULL)
  rs.elist <- extract(x = brk, y = spols, fun = NULL)

  ## Comparison
  expect_true(all(vx.elist[[1]][,1] %in% rs.elist[[1]][,1]))
  expect_true(all(vx.elist[[1]][,2] %in% rs.elist[[1]][,2]))
  expect_true(all(vx.elist[[2]][,1] %in% rs.elist[[2]][,1]))
  expect_true(all(vx.elist[[2]][,2] %in% rs.elist[[2]][,2]))

  ## Make SpatialLines / sfc_MLS
  lines.sf <- st_sfc(st_linestring(matrix(c(0.1,0,
                                            1,1), nrow = 2, ncol = 2, byrow=TRUE)))
  lines.sl <- as(lines.sf, 'Spatial')

  ## Extract raw values as lists using sp
  vx.elist <- vx$extract(sp = lines.sl, fun = NULL)
  rs.elist <- extract(x = brk, y = lines.sl, fun = NULL)

  ## Comparison using sp
  expect_true(all(vx.elist[[1]][,1] %in% rs.elist[[1]][,1]))
  expect_true(all(vx.elist[[1]][,2] %in% rs.elist[[1]][,2]))

  ## Extract raw values as lists using sf
  vx.elist <- vx$extract(sp=lines.sf, fun = NULL)
  rs.elist <- extract(x = brk, y = lines.sl, fun = NULL)

  ## Comparison
  expect_true(all(vx.elist[[1]][,1] %in% rs.elist[[1]][,1]))
  expect_true(all(vx.elist[[1]][,2] %in% rs.elist[[1]][,2]))

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
  coord <- matrix(c(0,0,1,1, 0.5, 0.55), nrow=3, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=c(0.5, 0.5, 0.04), byid = TRUE)

  ## Extract raw values as data-frame
  vx.raw.df <- vx$extract(sp=spols, fun = NULL, df=TRUE)
  rs.raw.df <- extract(x = brk, y = spols, fun = NULL, df=TRUE, small=FALSE)
  dimnames(rs.raw.df) <- dimnames(vx.raw.df)

  ## Extract aggregated values as data-frame
  vx.foo.df <- vx$extract(sp=spols, fun = mean, df=TRUE)
  rs.foo.df <- extract(x = brk, y = spols, fun = mean, df=TRUE, small=FALSE)
  dimnames(rs.foo.df) <- dimnames(vx.foo.df)

  ## Comparison
  expect_true(all(vx.raw.df$X1 %in% rs.raw.df$X1) & all(vx.raw.df$X2 %in% rs.raw.df$X2))
  expect_true(all(vx.foo.df == rs.foo.df| is.na(vx.foo.df)==is.na(rs.foo.df)))
})

test_that("small option works with small polygons", {

  ## Make VeloxRaster with two bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  mat2 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1), raster(mat2))
  vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons
  coord <- matrix(c(0.1,0.1,0.9,0.9), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.0000001, byid = TRUE)

  ## Extract raw values as lists
  vx.elist <- vx$extract(sp = spols, fun = NULL, small = TRUE)
  rs.elist <- extract(x = brk, y = spols, fun = NULL, small = TRUE)

  ## Comparison
  expect_true(all(vx.elist[[1]][,1] %in% rs.elist[[1]][,1]))
  expect_true(all(vx.elist[[1]][,2] %in% rs.elist[[1]][,2]))
  expect_true(all(vx.elist[[2]][,1] %in% rs.elist[[2]][,1]))
  expect_true(all(vx.elist[[2]][,2] %in% rs.elist[[2]][,2]))
})
