# geom_crs_utils.R
# Layer 2: CRS utilities for spatintegrate
#
# Three functions:
#
#   ensure_projected(x)              — auto-project geographic CRS to UTM
#   ensure_same_crs(crs_a, crs_b)   — error if two CRS objects differ
#   assert_projected(x, arg_name)   — error if x is not in a projected CRS
#
# All spatintegrate geometry functions require projected CRS inputs.
# ensure_projected() is the user-facing entry point for converting inputs.
# assert_projected() is the internal guard used at the top of each function.


# ------------------------------------------------------------------------------

#' Ensure an sf Object is in a Projected CRS
#'
#' Converts an \code{sf} or \code{sfc} object to a projected coordinate
#' reference system if it is currently in a geographic CRS (lon/lat). If the
#' input is already projected, it is returned unchanged.
#'
#' UTM zone is selected automatically based on the centroid of the geometry.
#' If centroid computation fails, falls back to \code{fallback_epsg}.
#'
#' @param x An \code{sf} or \code{sfc} object.
#' @param fallback_epsg Integer EPSG code used if UTM zone selection fails.
#'   Default \code{6933} (Equal Earth).
#' @param verbose Logical. If \code{TRUE}, prints a message when a
#'   transformation is applied. Default \code{TRUE}.
#'
#' @return An \code{sf} or \code{sfc} object in a projected CRS. If the input
#'   was already projected, the input is returned unchanged (no copy).
#'
#' @seealso [assert_projected()], [ensure_same_crs()]
#'
#' @examples
#' x      <- sf::st_sfc(sf::st_point(c(-71, 42)), crs = 4326)
#' x_proj <- ensure_projected(x, verbose = FALSE)
#' sf::st_is_longlat(x_proj)  # FALSE
#'
#' @export
ensure_projected <- function(x, fallback_epsg = 6933L, verbose = TRUE) {
  if (!inherits(x, c("sf", "sfc"))) {
    stop("`x` must be an sf or sfc object.")
  }

  crs <- sf::st_crs(x)
  if (is.na(crs)) {
    stop("Object has no CRS. Assign a CRS before calling ensure_projected().")
  }

  # Already projected — return unchanged
  if (!isTRUE(sf::st_is_longlat(x))) {
    return(x)
  }

  # Try to determine UTM zone from centroid
  centroid <- tryCatch(
    {
      cxy <- sf::st_coordinates(sf::st_centroid(sf::st_union(x)))[1, ]
      list(lon = cxy[1L], lat = cxy[2L])
    },
    error = function(e) NULL
  )

  if (is.null(centroid)) {
    if (verbose) {
      message(
        "Could not compute centroid for UTM selection; ",
        "falling back to EPSG:", fallback_epsg
      )
    }
    return(sf::st_transform(x, fallback_epsg))
  }

  # Compute UTM EPSG from lon/lat
  utm_zone <- floor((centroid$lon + 180) / 6) + 1L
  epsg     <- if (centroid$lat >= 0) 32600L + utm_zone else 32700L + utm_zone

  if (verbose) {
    message(sprintf(
      "Projecting from geographic CRS to UTM zone %d (EPSG:%d).",
      utm_zone, epsg
    ))
  }

  x_proj <- tryCatch(
    sf::st_transform(x, epsg),
    error = function(e) {
      if (verbose) {
        message(
          "UTM transform failed; falling back to EPSG:", fallback_epsg
        )
      }
      sf::st_transform(x, fallback_epsg)
    }
  )

  x_proj
}


# ------------------------------------------------------------------------------

#' Assert That an sf Object is in a Projected CRS
#'
#' Checks that \code{x} is in a projected (non-geographic) CRS and throws an
#' informative error if not. This is the internal CRS guard used at the top of
#' all \pkg{spatintegrate} geometry functions.
#'
#' @param x An \code{sf}, \code{sfc}, or \code{sfg} object.
#' @param arg_name Character scalar. Name of the argument as it appears in the
#'   calling function, used in the error message. Default \code{"x"}.
#'
#' @return Invisibly returns \code{TRUE} if the CRS is projected. Otherwise
#'   throws an error.
#'
#' @seealso [ensure_projected()], [ensure_same_crs()]
#'
#' @examples
#' proj <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 32619)
#' assert_projected(proj)  # silent
#'
#' \dontrun{
#' geo <- sf::st_sfc(sf::st_point(c(-71, 42)), crs = 4326)
#' assert_projected(geo)  # error
#' }
#'
#' @export
assert_projected <- function(x, arg_name = "x") {
  if (!inherits(x, c("sf", "sfc", "sfg"))) {
    stop(sprintf("`%s` must be an sf, sfc, or sfg object.", arg_name))
  }

  # sfg has no CRS — skip check (coordinates are assumed correct by caller)
  if (inherits(x, "sfg")) {
    return(invisible(TRUE))
  }

  crs <- sf::st_crs(x)

  if (is.na(crs)) {
    stop(sprintf(
      "`%s` has no CRS. Assign a projected CRS before calling this function.",
      arg_name
    ))
  }

  if (isTRUE(sf::st_is_longlat(x))) {
    stop(sprintf(
      "`%s` is in a geographic CRS (lon/lat). ",
      arg_name
    ),
    "spatintegrate requires a projected CRS (units of meters). ",
    "Use ensure_projected() to convert first."
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------

#' Ensure Two CRS Objects Are Identical
#'
#' Checks whether two CRS definitions are identical and throws an informative
#' error if not. Inputs must be \code{sf::st_crs()} objects, not raw sf
#' objects — extract the CRS first with \code{sf::st_crs(x)}.
#'
#' @param crs_a,crs_b CRS objects created by \code{sf::st_crs()}.
#' @param context Character scalar describing the operation being checked,
#'   used in the error message. Default \code{"operation"}.
#'
#' @return Invisibly returns \code{TRUE} if the CRSs match. Throws an error
#'   otherwise.
#'
#' @seealso [ensure_projected()], [assert_projected()]
#'
#' @examples
#' crs1 <- sf::st_crs(32619)
#' crs2 <- sf::st_crs(32619)
#' ensure_same_crs(crs1, crs2)  # silent
#'
#' \dontrun{
#' crs3 <- sf::st_crs(32618)
#' ensure_same_crs(crs1, crs3, context = "intersection")  # error
#' }
#'
#' @export
ensure_same_crs <- function(crs_a, crs_b, context = "operation") {
  if (!inherits(crs_a, "crs") || !inherits(crs_b, "crs")) {
    stop(
      "Both arguments must be CRS objects from sf::st_crs(). ",
      "Call sf::st_crs(x) on your sf/sfc object first."
    )
  }

  if (!isTRUE(crs_a == crs_b)) {
    # Safe label extraction — fall back to wkt snippet if input is unavailable
    label_a <- if (!is.null(crs_a$input)) crs_a$input else "<unknown>"
    label_b <- if (!is.null(crs_b$input)) crs_b$input else "<unknown>"

    stop(sprintf(
      "CRS mismatch during %s:\n  A: %s\n  B: %s\n%s",
      context,
      label_a,
      label_b,
      "Transform inputs to a common CRS before proceeding."
    ))
  }

  invisible(TRUE)
}
