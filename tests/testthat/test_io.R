library(velox)
library(raster)
context("IO")

test_that("data type recognition works", {

  make_vx <- function(val) {
    velox(raster(matrix(rep(val, 9), 3, 3)))
  }

  byte_vx <- make_vx(1)
  expect_equal(byte_vx$get_data_type(), "Byte")

  uint16_vx <- make_vx(255 + 1)
  expect_equal(uint16_vx$get_data_type(), "UInt16")

  uint32_vx <- make_vx(65534 + 1)
  expect_equal(uint32_vx$get_data_type(), "UInt32")

  float32_vx <- make_vx(4294967296 + 1)
  expect_equal(float32_vx$get_data_type(), "Float32")

  float64_vx <- make_vx(3.4e+39)
  expect_equal(float64_vx$get_data_type(), "Float64")

  int16_vx <- make_vx(-10)
  expect_equal(int16_vx$get_data_type(), "Int16")

  int32_vx <- make_vx(-60000)
  expect_equal(int32_vx$get_data_type(), "Int32")

  float32_vx <- make_vx(-2147483647 - 1)
  expect_equal(float32_vx$get_data_type(), "Float32")

  float64_vx <- make_vx(-3.4e+39)
  expect_equal(float64_vx$get_data_type(), "Float64")
})


test_that("data type assignment works when casting as rasterlayer", {

  make_vx <- function(val) {
    velox(raster(matrix(rep(val, 9), 3, 3)))
  }

  byte_vx <- make_vx(1)
  byte_ras <- byte_vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(byte_ras), "INT1U")

  vx <- make_vx(255 + 1)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "INT2U")

  vx <- make_vx(65534 + 1)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "FLT4S")

  vx <- make_vx(4294967296 + 1)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "FLT4S")

  vx <- make_vx(3.4e+39)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "FLT8S")

  vx <- make_vx(-10)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "INT2S")

  vx <- make_vx(-60000)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "INT4S")

  vx <- make_vx(-2147483647 - 1)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "FLT4S")

  vx <- make_vx(-3.4e+39)
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  expect_equal(dataType(rs), "FLT8S")

})

test_that("saving veloxraster and reading vx rasters works", {

  # Prepare
  make_vx <- function(val) {
    velox(raster(matrix(rep(val, 9), 3, 3)))
  }
  fp <- file.path(tempdir(), 'test.tif')

  # Check byte
  vx <- make_vx(1)
  vx$write(fp, TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  writeRaster(x = rs, filename = fp, datatype=dataType(rs), overwrite=TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))

  # Check uint16
  vx <- make_vx(256)
  vx$write(fp, TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  writeRaster(x = rs, filename = fp, datatype=dataType(rs), overwrite=TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))

  # Check uint32
  vx <- make_vx(65534 + 1)
  vx$write(fp, TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  writeRaster(x = rs, filename = fp, datatype=dataType(rs), overwrite=TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal('double', storage.mode(vx2$rasterbands[[1]]))

  # Check float32
  vx <- make_vx(4294967296 + 1)
  vx$write(fp, TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  writeRaster(x = rs, filename = fp, datatype=dataType(rs), overwrite=TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))

  # Check float64
  vx <- make_vx(-3.4e+39)
  vx$write(fp, TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))
  rs <- vx$as.RasterLayer(assign_data_type = TRUE)
  writeRaster(x = rs, filename = fp, datatype=dataType(rs), overwrite=TRUE)
  vx2 <- velox(fp)
  expect_equal(vx$rasterbands[[1]], vx2$rasterbands[[1]])
  expect_equal(storage.mode(vx$rasterbands[[1]]), storage.mode(vx2$rasterbands[[1]]))

  # Clean up
  unlink(fp)
})


