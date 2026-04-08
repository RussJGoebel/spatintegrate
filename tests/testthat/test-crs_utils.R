# tests/testthat/test-geom_crs_utils.R
#
# Tests for CRS utilities:
#   ensure_projected()
#   assert_projected()
#   ensure_same_crs()

library(testthat)
library(sf)

# ── Helpers ───────────────────────────────────────────────────────────────────

# A single point in WGS84 (Boston)
geo_point <- function() {
  sf::st_sfc(sf::st_point(c(-71.06, 42.36)), crs = 4326)
}

# A single point in UTM 19N (already projected)
proj_point <- function() {
  sf::st_sfc(sf::st_point(c(330000, 4690000)), crs = 32619)
}

# A small polygon in WGS84
geo_polygon <- function() {
  sf::st_sfc(sf::st_polygon(list(matrix(
    c(-71.1, 42.3,  -71.0, 42.3,  -71.0, 42.4,  -71.1, 42.4,  -71.1, 42.3),
    ncol = 2, byrow = TRUE
  ))), crs = 4326)
}

# A small polygon already in UTM 19N
proj_polygon <- function() {
  sf::st_sfc(sf::st_polygon(list(matrix(
    c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE
  ))), crs = 32619)
}


# ══════════════════════════════════════════════════════════════════════════════
# ensure_projected()
# ══════════════════════════════════════════════════════════════════════════════

test_that("ensure_projected: already-projected input returned unchanged", {
  p    <- proj_point()
  out  <- ensure_projected(p, verbose = FALSE)
  expect_identical(out, p)
})

test_that("ensure_projected: geographic input is transformed to projected", {
  out <- ensure_projected(geo_point(), verbose = FALSE)
  expect_false(isTRUE(sf::st_is_longlat(out)))
})

test_that("ensure_projected: output CRS is a UTM zone for standard lon/lat", {
  out  <- ensure_projected(geo_point(), verbose = FALSE)
  epsg <- sf::st_crs(out)$epsg
  # UTM EPSGs are 32601-32660 (N) or 32701-32760 (S)
  expect_true(
    (epsg >= 32601L && epsg <= 32660L) ||
      (epsg >= 32701L && epsg <= 32760L)
  )
})

test_that("ensure_projected: southern hemisphere gets correct UTM band", {
  # Sydney, Australia — should get a 327xx EPSG
  sydney <- sf::st_sfc(sf::st_point(c(151.2, -33.9)), crs = 4326)
  out    <- ensure_projected(sydney, verbose = FALSE)
  epsg   <- sf::st_crs(out)$epsg
  expect_true(epsg >= 32700L && epsg <= 32760L)
})

test_that("ensure_projected: works on sf object not just sfc", {
  poly_sf <- sf::st_sf(geometry = geo_polygon())
  out     <- ensure_projected(poly_sf, verbose = FALSE)
  expect_false(isTRUE(sf::st_is_longlat(out)))
  expect_s3_class(out, "sf")
})

test_that("ensure_projected: verbose = TRUE emits a message", {
  expect_message(ensure_projected(geo_point(), verbose = TRUE))
})

test_that("ensure_projected: verbose = FALSE is silent", {
  expect_silent(ensure_projected(geo_point(), verbose = FALSE))
})

test_that("ensure_projected: already-projected input is silent regardless of verbose", {
  expect_silent(ensure_projected(proj_point(), verbose = TRUE))
})

test_that("ensure_projected: bad input errors informatively", {
  expect_error(ensure_projected("not sf"),    "sf or sfc object")
  expect_error(ensure_projected(list()),      "sf or sfc object")
  expect_error(ensure_projected(matrix(1:4, 2)), "sf or sfc object")
})

test_that("ensure_projected: no-CRS input errors informatively", {
  nocrs <- sf::st_sfc(sf::st_point(c(0, 0)))
  expect_error(ensure_projected(nocrs), "no CRS")
})


# ══════════════════════════════════════════════════════════════════════════════
# assert_projected()
# ══════════════════════════════════════════════════════════════════════════════

test_that("assert_projected: projected sfc returns TRUE invisibly", {
  result <- assert_projected(proj_point())
  expect_true(result)
})

test_that("assert_projected: projected sf returns TRUE invisibly", {
  result <- assert_projected(sf::st_sf(geometry = proj_polygon()))
  expect_true(result)
})

test_that("assert_projected: sfg skips CRS check and returns TRUE", {
  sfg <- sf::st_point(c(0, 0))
  expect_true(assert_projected(sfg))
})

test_that("assert_projected: geographic sfc errors", {
  expect_error(assert_projected(geo_point()), "geographic CRS")
})

test_that("assert_projected: geographic sf errors", {
  expect_error(
    assert_projected(sf::st_sf(geometry = geo_polygon())),
    "geographic CRS"
  )
})

test_that("assert_projected: no-CRS input errors", {
  nocrs <- sf::st_sfc(sf::st_point(c(0, 0)))
  expect_error(assert_projected(nocrs), "no CRS")
})

test_that("assert_projected: arg_name appears in error message", {
  expect_error(
    assert_projected(geo_point(), arg_name = "soundings"),
    "soundings"
  )
})

test_that("assert_projected: non-sf input errors informatively", {
  expect_error(assert_projected("string"),    "sf, sfc, or sfg")
  expect_error(assert_projected(1:10),        "sf, sfc, or sfg")
  expect_error(assert_projected(list(a = 1)), "sf, sfc, or sfg")
})


# ══════════════════════════════════════════════════════════════════════════════
# ensure_same_crs()
# ══════════════════════════════════════════════════════════════════════════════

test_that("ensure_same_crs: identical CRS returns TRUE invisibly", {
  crs <- sf::st_crs(32619)
  result <- ensure_same_crs(crs, crs)
  expect_true(result)
})

test_that("ensure_same_crs: same EPSG different construction is equal", {
  crs_a <- sf::st_crs(32619)
  crs_b <- sf::st_crs("EPSG:32619")
  result <- ensure_same_crs(crs_a, crs_b)
  expect_true(result)
})

test_that("ensure_same_crs: different CRS errors", {
  expect_error(
    ensure_same_crs(sf::st_crs(32619), sf::st_crs(32618)),
    "CRS mismatch"
  )
})

test_that("ensure_same_crs: error message includes context", {
  expect_error(
    ensure_same_crs(sf::st_crs(32619), sf::st_crs(4326),
                    context = "triangulation"),
    "triangulation"
  )
})

test_that("ensure_same_crs: error message includes both CRS labels", {
  err <- tryCatch(
    ensure_same_crs(sf::st_crs(32619), sf::st_crs(32618)),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "32619")
  expect_match(err, "32618")
})

test_that("ensure_same_crs: non-crs inputs error informatively", {
  crs <- sf::st_crs(32619)
  expect_error(ensure_same_crs("32619", crs),  "sf::st_crs\\(\\)")
  expect_error(ensure_same_crs(crs, 32619),    "sf::st_crs\\(\\)")
  expect_error(ensure_same_crs(32619, 32619),  "sf::st_crs\\(\\)")
})

test_that("ensure_same_crs: NA CRS objects error", {
  crs_na <- sf::st_crs(NA)
  crs_ok <- sf::st_crs(32619)
  # NA CRS != valid CRS — should error on mismatch
  expect_error(ensure_same_crs(crs_na, crs_ok))
})
