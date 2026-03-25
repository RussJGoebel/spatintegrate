# tests/testthat/test-qmc_utils.R
#
# Tests for Layer 1 QMC sampling primitives:
#   generate_qmc_unit_square()
#   map_unit_square_to_triangle()

library(testthat)
library(sf)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Point-in-triangle test using barycentric coordinates
point_in_triangle <- function(pts, tri) {
  v0 <- tri[3, ] - tri[1, ]
  v1 <- tri[2, ] - tri[1, ]
  apply(pts, 1, function(p) {
    v2    <- p - tri[1, ]
    dot00 <- sum(v0 * v0); dot01 <- sum(v0 * v1)
    dot02 <- sum(v0 * v2); dot11 <- sum(v1 * v1); dot12 <- sum(v1 * v2)
    inv   <- 1 / (dot00 * dot11 - dot01 * dot01)
    u     <- (dot11 * dot02 - dot01 * dot12) * inv
    v     <- (dot00 * dot12 - dot01 * dot02) * inv
    (u >= -1e-9) && (v >= -1e-9) && (u + v <= 1 + 1e-9)
  })
}


# ══════════════════════════════════════════════════════════════════════════════
# generate_qmc_unit_square()
# ══════════════════════════════════════════════════════════════════════════════

test_that("generate_qmc_unit_square: output shape is correct", {
  for (n in c(1, 2, 10, 64, 256)) {
    pts <- generate_qmc_unit_square(n)
    expect_true(is.matrix(pts),  info = paste("n =", n))
    expect_equal(nrow(pts), n,   info = paste("n =", n))
    expect_equal(ncol(pts), 2L,  info = paste("n =", n))
    expect_true(is.numeric(pts), info = paste("n =", n))
  }
})

test_that("generate_qmc_unit_square: all values in [0, 1]", {
  pts <- generate_qmc_unit_square(1024)
  expect_true(all(pts >= 0))
  expect_true(all(pts <= 1))
})

test_that("generate_qmc_unit_square: n = 1 works", {
  pts <- generate_qmc_unit_square(1)
  expect_equal(dim(pts), c(1L, 2L))
})

test_that("generate_qmc_unit_square: non-power-of-2 n works", {
  pts <- generate_qmc_unit_square(100)
  expect_equal(nrow(pts), 100L)
  expect_true(all(pts >= 0) && all(pts <= 1))
})

test_that("generate_qmc_unit_square: sequence is deterministic", {
  expect_identical(generate_qmc_unit_square(64), generate_qmc_unit_square(64))
})

test_that("generate_qmc_unit_square: low discrepancy vs random", {
  # Divide [0,1]^2 into 4x4 = 16 cells; QMC stddev of cell counts < random
  n        <- 256
  pts_qmc  <- generate_qmc_unit_square(n)
  set.seed(42)
  pts_rand <- matrix(runif(n * 2), ncol = 2)

  cell <- function(pts) {
    xi <- pmin(floor(pts[, 1] * 4), 3)
    yi <- pmin(floor(pts[, 2] * 4), 3)
    xi * 4 + yi
  }
  expect_lt(
    sd(tabulate(cell(pts_qmc)  + 1, nbins = 16)),
    sd(tabulate(cell(pts_rand) + 1, nbins = 16))
  )
})

test_that("generate_qmc_unit_square: bad inputs error informatively", {
  expect_error(generate_qmc_unit_square(0),         "`n` must be a positive integer")
  expect_error(generate_qmc_unit_square(-1),        "`n` must be a positive integer")
  expect_error(generate_qmc_unit_square(1.5),       "`n` must be a positive integer")
  expect_error(generate_qmc_unit_square("10"),      "`n` must be a positive integer")
  expect_error(generate_qmc_unit_square(c(10, 20)), "`n` must be a positive integer")
  expect_error(generate_qmc_unit_square(NA_real_),  "`n` must be a positive integer")
})


# ══════════════════════════════════════════════════════════════════════════════
# map_unit_square_to_triangle()
# ══════════════════════════════════════════════════════════════════════════════

test_that("map_unit_square_to_triangle: output shape matches input", {
  tri <- matrix(c(0,0, 1,0, 0,1), ncol = 2, byrow = TRUE)
  for (n in c(1, 16, 256)) {
    pts    <- generate_qmc_unit_square(n)
    mapped <- map_unit_square_to_triangle(pts, tri)
    expect_equal(dim(mapped), c(n, 2L), info = paste("n =", n))
  }
})

test_that("map_unit_square_to_triangle: all points inside standard triangle", {
  tri    <- matrix(c(0,0, 1,0, 0,1), ncol = 2, byrow = TRUE)
  pts    <- generate_qmc_unit_square(512)
  mapped <- map_unit_square_to_triangle(pts, tri)
  expect_true(all(point_in_triangle(mapped, tri)))
})

test_that("map_unit_square_to_triangle: all points inside arbitrary triangle", {
  tri    <- matrix(c(1,2, 4,1, 3,5), ncol = 2, byrow = TRUE)
  pts    <- generate_qmc_unit_square(512)
  mapped <- map_unit_square_to_triangle(pts, tri)
  expect_true(all(point_in_triangle(mapped, tri)))
})

test_that("map_unit_square_to_triangle: centroid convergence", {
  tri           <- matrix(c(1,2, 4,1, 3,5), ncol = 2, byrow = TRUE)
  true_centroid <- colMeans(tri)
  mapped        <- map_unit_square_to_triangle(generate_qmc_unit_square(1024), tri)
  expect_lt(sqrt(sum((colMeans(mapped) - true_centroid)^2)), 0.05)
})

test_that("map_unit_square_to_triangle: single point works", {
  tri    <- matrix(c(0,0, 1,0, 0,1), ncol = 2, byrow = TRUE)
  mapped <- map_unit_square_to_triangle(matrix(c(0.2, 0.3), nrow = 1), tri)
  expect_equal(dim(mapped), c(1L, 2L))
  expect_true(all(point_in_triangle(mapped, tri)))
})

test_that("map_unit_square_to_triangle: same qmc reused across triangles gives correct containment", {
  # Core use pattern: one QMC sequence, mapped into multiple different triangles
  qmc <- generate_qmc_unit_square(64)

  t1 <- matrix(c(0,0, 1,0, 0,1),   ncol = 2, byrow = TRUE)
  t2 <- matrix(c(2,0, 4,0, 2,2),   ncol = 2, byrow = TRUE)
  t3 <- matrix(c(0,3, 1,3, 0.5,4), ncol = 2, byrow = TRUE)

  expect_true(all(point_in_triangle(map_unit_square_to_triangle(qmc, t1), t1)))
  expect_true(all(point_in_triangle(map_unit_square_to_triangle(qmc, t2), t2)))
  expect_true(all(point_in_triangle(map_unit_square_to_triangle(qmc, t3), t3)))
})

test_that("map_unit_square_to_triangle: collinear vertices do not crash", {
  # Zero-area triangle — mapping runs without error even if result is degenerate
  tri <- matrix(c(0,0, 1,0, 2,0), ncol = 2, byrow = TRUE)
  expect_no_error(
    map_unit_square_to_triangle(generate_qmc_unit_square(16), tri)
  )
})

test_that("map_unit_square_to_triangle: bad inputs error informatively", {
  tri <- matrix(c(0,0, 1,0, 0,1), ncol = 2, byrow = TRUE)
  pts <- generate_qmc_unit_square(16)

  expect_error(
    map_unit_square_to_triangle(as.vector(pts), tri),
    "`qmc_points` must be an n x 2 numeric matrix"
  )
  expect_error(
    map_unit_square_to_triangle(matrix(1:6, ncol = 3), tri),
    "`qmc_points` must be an n x 2 numeric matrix"
  )
  expect_error(
    map_unit_square_to_triangle(pts, matrix(c(0,0,1,0), ncol = 2)),
    "`triangle_coords` must be a 3 x 2 numeric matrix"
  )
  expect_error(
    map_unit_square_to_triangle(pts, matrix(1:8, ncol = 2)),
    "`triangle_coords` must be a 3 x 2 numeric matrix"
  )
  expect_error(
    map_unit_square_to_triangle(matrix(as.character(pts), ncol = 2), tri),
    "must be numeric"
  )
})
