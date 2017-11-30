library(sp)
library(sf)
library(rgeos)
library(velox)
library(raster)
context("Rasterize")

test_that("rasterize works with polygons", {

  ## Make VeloxRaster with one bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1))
  vx <- velox(list(mat1), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialPolygons / sfc_MPG
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)
  spdf <- SpatialPolygonsDataFrame(spols, data.frame(id=1:2), FALSE)
  pols.sf <- st_as_sf(spdf)

  ## Rasterize using sp
  vx$rasterize(sp=spdf, field = 'id', band = 1, background = 0)
  rs.emat <- rasterize(x = spdf, y = brk, field = 'id', background = 0)
  vx.mat <- vx$as.matrix()
  rs.mat <- as.matrix(rs.emat)

  ## Comparison using sp
  expect_equal(vx.mat, rs.mat)

  ## Rasterize using sf
  vx$rasterize(sp=pols.sf, field = 'id', band = 1, background = 0)
  rs.emat <- rasterize(x = spdf, y = brk, field = 'id', background = 0)
  vx.mat <- vx$as.matrix()
  rs.mat <- as.matrix(rs.emat)

  ## Comparison using sf
  expect_equal(vx.mat, rs.mat)

})

test_that("rasterize works with lines", {

  ## Make VeloxRaster with one bands
  set.seed(0)
  mat1 <- matrix(rnorm(100), 10, 10)
  brk <- brick(raster(mat1))
  vx <- velox(list(mat1), extent=c(0,1,0,1), res=c(0.1,0.1),
              crs="+proj=longlat +datum=WGS84 +no_defs")

  ## Make SpatialLines / sfc_MLS
  lines.sf <- st_sfc(st_linestring(matrix(c(0.1,0,
                                            1,1), nrow = 2, ncol = 2, byrow=TRUE)))
  lines.sl <- as(lines.sf, 'Spatial')
  sldf <- SpatialLinesDataFrame(lines.sl, data.frame(id = 1), FALSE)
  lines.sf <- st_as_sf(sldf)

  ## Rasterize using sp
  vx$rasterize(sp=sldf, field = 'id', band = 1, background = 0)
  rs.emat <- rasterize(x = sldf, y = brk, field = 'id', background = 0)
  vx.mat <- vx$as.matrix()
  rs.mat <- as.matrix(rs.emat)

  ## Comparison using sp
  expect_equal(vx.mat, rs.mat)

  ## Rasterize using sf
  vx$rasterize(sp=lines.sf, field = 'id', band = 1, background = 0)
  rs.emat <- rasterize(x = sldf, y = brk, field = 'id', background = 0)
  vx.mat <- vx$as.matrix()
  rs.mat <- as.matrix(rs.emat)

  ## Comparison using sf
  expect_equal(vx.mat, rs.mat)

})
