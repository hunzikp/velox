library(sp)
library(sf)
library(rgeos)
library(velox)
library(raster)
context("BoostGeometries")

test_that("bg_intersects works", {

  ## Make polygons
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coords=coord)
  spols <- gBuffer(spgeom=spoint, width=0.5, byid = TRUE)
  pols.sf <- st_geometry(st_as_sf(spols))
  pols.boost <- boost(pols.sf)

  ## Make Line
  coord <- matrix(c(0,0,1,1), 2, 2, byrow = TRUE)
  line.sf <- st_sfc(st_linestring(coord))
  line.boost <- boost(line.sf)

  ## Make points
  coord <- matrix(c(0.25, 0.25, 0.75, 0.75), 2, 2, byrow = TRUE)
  spoint <- SpatialPoints(coord)
  points.sf <- st_geometry(st_as_sf(spoint))
  points.boost <- boost(points.sf)

  ## Intersect polygons
  # w polygons
  bp.ls <- bg_intersects(pols.boost, pols.boost)
  sf.ls <- st_intersects(pols.sf, pols.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)
  # w lines
  bp.ls <- bg_intersects(pols.boost, line.boost)
  sf.ls <- st_intersects(pols.sf, line.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)
  # w points
  bp.ls <- bg_intersects(pols.boost, points.boost)
  sf.ls <- st_intersects(pols.sf, points.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)

  ## Intersect lines
  # w polygons
  bp.ls <- bg_intersects(line.boost, pols.boost)
  sf.ls <- st_intersects(line.sf, pols.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)
  # w lines
  bp.ls <- bg_intersects(line.boost, line.boost)
  sf.ls <- st_intersects(line.sf, line.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)
  # w points
  bp.ls <- bg_intersects(line.boost, points.boost)
  sf.ls <- st_intersects(line.sf, points.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)

  ## Intersect points
  # w polygons
  bp.ls <- bg_intersects(points.boost, pols.boost)
  sf.ls <- st_intersects(points.sf, pols.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)
  # w lines
  bp.ls <- bg_intersects(points.boost, line.boost)
  sf.ls <- st_intersects(points.sf, line.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)
  # w points
  bp.ls <- bg_intersects(points.boost, points.boost)
  sf.ls <- st_intersects(points.sf, points.sf)
  sf.ls <- lapply(sf.ls, function(x) x)
  expect_equal(bp.ls, sf.ls)

})
