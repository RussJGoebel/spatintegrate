# tests/testthat/test-geom_triangulate.R
#
# Tests for Layer 2 triangulation utilities:
#   triangulate_sf()
#   get_triangle_coords()

library(testthat)
library(sf)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Projected CRS for all test geometries (UTM zone 19N — Boston area)
TEST_CRS <- 32619

make_sfc <- function(coords_matrix, crs = TEST_CRS) {
  closed <- rbind(coords_matrix, coords_matrix[1, ])
  sf::st_sfc(sf::st_polygon(list(closed)), crs = crs)
}

# Simple convex shapes
unit_square <- function() {
  make_sfc(matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE))
}

right_triangle <- function() {
  make_sfc(matrix(c(0,0, 2,0, 0,2), ncol = 2, byrow = TRUE))
}

# Non-convex L-shape
l_shape <- function() {
  make_sfc(matrix(c(0,0, 2,0, 2,1, 1,1, 1,2, 0,2), ncol = 2, byrow = TRUE))
}

# Multipolygon (two separate squares)
two_squares <- function() {
  s1 <- sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE)))
  s2 <- sf::st_polygon(list(matrix(c(3,0,4,0,4,1,3,1,3,0), ncol=2, byrow=TRUE)))
  sf::st_sfc(sf::st_multipolygon(list(s1, s2)), crs = TEST_CRS)
}

# Check all triangle centroids are inside the original polygon
all_centroids_inside <- function(tris, original_sfc) {
  centroids <- sf::st_centroid(tris)
  covered   <- sf::st_covered_by(centroids, sf::st_union(original_sfc), sparse = FALSE)
  all(covered)
}

# Check that a geometry is a valid triangle (3 unique exterior vertices)
is_triangle <- function(sfg) {
  coords <- sf::st_coordinates(sfg)
  if ("L2" %in% colnames(coords)) {
    ext <- coords[coords[, "L2"] == 1L, c("X", "Y"), drop = FALSE]
  } else {
    ext <- coords[, c("X", "Y"), drop = FALSE]
  }
  # Remove closing duplicate
  if (nrow(ext) >= 2 && all(ext[1,] == ext[nrow(ext),])) {
    ext <- ext[-nrow(ext), , drop = FALSE]
  }
  nrow(ext) == 3L
}


# ══════════════════════════════════════════════════════════════════════════════
# triangulate_sf()
# ══════════════════════════════════════════════════════════════════════════════

# ── Output structure ──────────────────────────────────────────────────────────

test_that("triangulate_sf: returns an sfc", {
  tris <- triangulate_sf(unit_square())
  expect_s3_class(tris, "sfc")
})

test_that("triangulate_sf: all returned geometries are triangles", {
  tris <- triangulate_sf(l_shape())
  expect_true(all(sapply(tris, is_triangle)))
})

test_that("triangulate_sf: preserves input CRS", {
  sq   <- unit_square()
  tris <- triangulate_sf(sq)
  expect_equal(sf::st_crs(tris), sf::st_crs(sq))
})

# ── Triangle counts ───────────────────────────────────────────────────────────

test_that("triangulate_sf: square produces 2 triangles", {
  # A convex quadrilateral always triangulates to n-2 = 2 triangles
  tris <- triangulate_sf(unit_square())
  expect_equal(length(tris), 2L)
})

test_that("triangulate_sf: a triangle input produces 1 triangle", {
  tris <- triangulate_sf(right_triangle())
  expect_equal(length(tris), 1L)
})

test_that("triangulate_sf: L-shape produces correct interior triangle count", {
  # L-shape has 6 vertices → Delaunay gives up to 4 triangles interior
  tris <- triangulate_sf(l_shape())
  expect_gte(length(tris), 3L)  # at least 3 interior triangles
  expect_lte(length(tris), 6L)  # sanity upper bound
})

# ── Interior guarantee ────────────────────────────────────────────────────────

test_that("triangulate_sf: all triangle centroids inside original polygon", {
  poly <- l_shape()
  tris <- triangulate_sf(poly)
  expect_true(all_centroids_inside(tris, poly))
})

test_that("triangulate_sf: area of triangles sums to area of original polygon", {
  poly      <- l_shape()
  tris      <- triangulate_sf(poly)
  area_poly <- as.numeric(sf::st_area(poly))
  area_tris <- sum(as.numeric(sf::st_area(tris)))
  expect_equal(area_tris, area_poly, tolerance = 1e-6)
})

test_that("triangulate_sf: no exterior triangles for non-convex polygon", {
  c_shape <- make_sfc(matrix(
    c(0,0, 3,0, 3,1, 1,1, 1,2, 3,2, 3,3, 0,3),
    ncol = 2, byrow = TRUE
  ))
  tris <- triangulate_sf(c_shape)
  expect_true(all_centroids_inside(tris, c_shape))
  area_orig <- as.numeric(sf::st_area(c_shape))
  area_tris <- sum(as.numeric(sf::st_area(tris)))
  expect_equal(area_tris, area_orig, tolerance = 1e-6)
})

# ── Input coercion ────────────────────────────────────────────────────────────

test_that("triangulate_sf: accepts sf object", {
  poly_sf <- sf::st_sf(geometry = unit_square())
  tris    <- triangulate_sf(poly_sf)
  expect_s3_class(tris, "sfc")
  expect_equal(length(tris), 2L)
})

test_that("triangulate_sf: accepts sfg object", {
  sfg  <- sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE)))
  sfc  <- sf::st_sfc(sfg, crs = TEST_CRS)
  tris <- triangulate_sf(sfc[[1L]] |> sf::st_sfc(crs = TEST_CRS))
  expect_s3_class(tris, "sfc")
})

test_that("triangulate_sf: accepts sfc object", {
  tris <- triangulate_sf(unit_square())
  expect_s3_class(tris, "sfc")
})

test_that("triangulate_sf: handles MULTIPOLYGON", {
  mp   <- two_squares()
  tris <- triangulate_sf(mp)
  expect_s3_class(tris, "sfc")
  expect_equal(length(tris), 4L)
  expect_true(all(sapply(tris, is_triangle)))
})

test_that("triangulate_sf: min_area = 0 keeps all non-degenerate triangles", {
  tris_default <- triangulate_sf(unit_square(), min_area = 0)
  expect_equal(length(tris_default), 2L)
})

test_that("triangulate_sf: min_area drops triangles below threshold", {
  # Unit square triangulates to 2 triangles each of area 0.5
  # Set min_area = 1 to drop both
  expect_warning(
    tris <- triangulate_sf(unit_square(), min_area = 1.0),
    "All triangles have area"
  )
  expect_equal(length(tris), 0L)
})

# ── CRS errors ───────────────────────────────────────────────────────────────

test_that("triangulate_sf: errors on geographic CRS", {
  poly_geo <- sf::st_sfc(
    sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE))),
    crs = 4326
  )
  expect_error(triangulate_sf(poly_geo), "geographic CRS")
})

test_that("triangulate_sf: errors on missing CRS", {
  poly_nocrs <- sf::st_sfc(
    sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE)))
  )
  expect_error(triangulate_sf(poly_nocrs), "no CRS")
})

# ── Bad inputs ────────────────────────────────────────────────────────────────

test_that("triangulate_sf: errors on non-sf input", {
  expect_error(triangulate_sf(matrix(1:6, ncol=2)), "sf, sfc, or sfg")
})

test_that("triangulate_sf: errors on non-polygon geometry type", {
  pts <- sf::st_sfc(sf::st_point(c(0, 0)), crs = TEST_CRS)
  expect_error(triangulate_sf(pts), "POLYGON or MULTIPOLYGON")
})

test_that("triangulate_sf: errors on bad min_area", {
  expect_error(triangulate_sf(unit_square(), min_area = -1),
               "non-negative numeric scalar")
  expect_error(triangulate_sf(unit_square(), min_area = NA),
               "non-negative numeric scalar")
})


# ══════════════════════════════════════════════════════════════════════════════
# get_triangle_coords()
# ══════════════════════════════════════════════════════════════════════════════

# Helper: make a single triangle sfg in projected CRS
single_triangle_sfg <- function() {
  coords <- matrix(c(0,0, 3,0, 1,2, 0,0), ncol = 2, byrow = TRUE)
  sf::st_polygon(list(coords))
}

single_triangle_sfc <- function() {
  sf::st_sfc(single_triangle_sfg(), crs = TEST_CRS)
}

# ── Output structure ──────────────────────────────────────────────────────────

test_that("get_triangle_coords: returns a 3x2 numeric matrix", {
  coords <- get_triangle_coords(single_triangle_sfg())
  expect_true(is.matrix(coords))
  expect_equal(dim(coords), c(3L, 2L))
  expect_true(is.numeric(coords))
})

test_that("get_triangle_coords: result has no row or column names", {
  coords <- get_triangle_coords(single_triangle_sfg())
  expect_null(rownames(coords))
  expect_null(colnames(coords))
})

test_that("get_triangle_coords: no closing duplicate in output", {
  coords <- get_triangle_coords(single_triangle_sfg())
  # First and last row should not be identical
  expect_false(isTRUE(all(coords[1, ] == coords[3, ])))
})

test_that("get_triangle_coords: coordinates match input vertices", {
  tri    <- matrix(c(1,2, 4,1, 3,5), ncol = 2, byrow = TRUE)
  sfg    <- sf::st_polygon(list(rbind(tri, tri[1,])))
  coords <- get_triangle_coords(sfg)
  expect_equal(coords, tri, tolerance = 1e-10)
})

# ── Input coercion ────────────────────────────────────────────────────────────

test_that("get_triangle_coords: accepts sfg", {
  coords <- get_triangle_coords(single_triangle_sfg())
  expect_equal(dim(coords), c(3L, 2L))
})

test_that("get_triangle_coords: accepts length-1 sfc", {
  coords <- get_triangle_coords(single_triangle_sfc())
  expect_equal(dim(coords), c(3L, 2L))
})

# ── Integration with triangulate_sf ──────────────────────────────────────────

test_that("get_triangle_coords: works on every triangle from triangulate_sf", {
  tris <- triangulate_sf(l_shape())
  for (i in seq_along(tris)) {
    coords <- get_triangle_coords(tris[[i]])
    expect_equal(dim(coords), c(3L, 2L),
                 info = paste("triangle", i))
    expect_true(is.numeric(coords),
                info = paste("triangle", i))
  }
})

test_that("get_triangle_coords + map_unit_square_to_triangle: points land inside triangle", {
  # End-to-end Layer 1+2 integration test.
  # Uses st_covered_by rather than st_within — the barycentric mapping can
  # place points exactly on triangle edges (boundary), which st_within
  # (strict interior only) would incorrectly reject.
  tris <- triangulate_sf(l_shape())
  qmc  <- generate_qmc_unit_square(32)

  for (i in seq_along(tris)) {
    coords <- get_triangle_coords(tris[[i]])
    mapped <- map_unit_square_to_triangle(qmc, coords)

    pts_sfc <- sf::st_sfc(
      lapply(seq_len(nrow(mapped)), function(k) {
        sf::st_point(mapped[k, ])
      }),
      crs = TEST_CRS
    )
    # st_covered_by: TRUE if point is inside OR on the boundary
    covered <- sf::st_covered_by(pts_sfc, tris[i], sparse = FALSE)[, 1]
    expect_true(all(covered), info = paste("triangle", i))
  }
})

# ── Bad inputs ────────────────────────────────────────────────────────────────

test_that("get_triangle_coords: errors on non-sfg input", {
  expect_error(get_triangle_coords(matrix(1:6, ncol=2)), "sfg")
  expect_error(get_triangle_coords("not a geometry"),    "sfg")
})

test_that("get_triangle_coords: errors on length > 1 sfc", {
  two_tris <- triangulate_sf(unit_square())  # 2 triangles
  expect_error(get_triangle_coords(two_tris), "single geometry")
})

test_that("get_triangle_coords: errors on non-triangle polygon", {
  # A square (4 vertices) should error
  sq_sfg <- sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE)))
  expect_error(get_triangle_coords(sq_sfg), "3 unique exterior vertices")
})
