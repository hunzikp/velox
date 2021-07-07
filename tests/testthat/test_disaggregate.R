library(velox)
library(testthat)

context("Disaggregate")

test_that("Disaggregate fails with bad parameters", {
  v <- velox(matrix(1:80, nrow=10, ncol=8), extent=c(0, 1, 0, 1), res=c(0.1, 0.125))

  # Non-integer factor
  expect_error(v$disaggregate(c(1, 2.2)))

  # No factor
  expect_error(v$disaggregate(c()))

  # Too many factors
  expect_error(v$disaggregate(c(1, 2, 3)))
})

test_that("Disaggregate doesn't change extent", {
  v <- velox(matrix(1:80, nrow=10, ncol=8), extent=c(0, 1, 0, 1), res=c(0.1, 0.125))

  v$disaggregate(c(13,17))

  expect_equal(v$extent, c(0, 1, 0, 1))
})

test_that("Disaggregate updates resolution", {
  v <- velox(matrix(1:80, nrow=10, ncol=8), extent=c(0, 1, 0, 1), res=c(0.1, 0.125))

  expect_equal(v$extent[c(1, 3)] + v$dim * v$res,
               v$extent[c(2, 4)])

  v$disaggregate(c(13, 17))

  expect_equal(v$extent[c(1, 3)] + v$dim * v$res,
               v$extent[c(2, 4)])
})

test_that("Disaggregate output dimensions are correct", {
  v <- velox(matrix(1:80, nrow=10, ncol=8), extent=c(0, 1, 0, 1), res=c(0.1, 0.125))

  expect_equal(v$dim, c(10, 8))

  v$disaggregate(c(13, 17))

  expect_equal(v$dim, c(130, 136))
})

test_that("Disaggregate output values are correct", {
  v <- velox(matrix(1:4, nrow=2, ncol=2), extent=c(0, 1, 0, 1), res=c(0.5, 0.5))

  v$disaggregate(c(2,3))

  expect_equal(v$rasterbands[[1]],
               rbind(c(1, 1, 1, 3, 3, 3),
                     c(1, 1, 1, 3, 3, 3),
                     c(2, 2, 2, 4, 4, 4),
                     c(2, 2, 2, 4, 4, 4))
               )
})

test_that("If only one factor is supplied, it is applied to both x and y dimensions", {
  v <- velox(matrix(1:12, nrow=3, ncol=4), extent=c(0, 3, 0, 4), res=c(1, 1))

  v$disaggregate(7)

  expect_equal(v$dim, c(21, 28))
})
