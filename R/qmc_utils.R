# qmc_utils.R
# Layer 1: QMC sampling primitives for spatintegrate
#
# Pure math — no sf dependency. Two functions only:
#
#   generate_qmc_unit_square()    — Sobol points in [0,1]^2
#   map_unit_square_to_triangle() — barycentric mapping into a triangle
#
# These are called inside integrate_one_sounding() (Layer 4), which handles
# the triangle iteration, area weighting, and basis evaluation. Layer 1 has
# no knowledge of soundings, basis functions, or triangulations.


# ------------------------------------------------------------------------------

#' Generate Sobol QMC Points in the Unit Square
#'
#' Returns an \eqn{n \times 2} matrix of quasi-Monte Carlo points in
#' \eqn{[0,1]^2} using a Sobol sequence. Points are suitable for mapping into
#' triangle regions via [map_unit_square_to_triangle()].
#'
#' A single call to this function is typically made once per sounding and the
#' result reused across all triangles in that sounding's triangulation — the
#' same Sobol sequence mapped independently into each triangle is correct
#' because the barycentric mapping is independent per triangle.
#'
#' @param n Positive integer. Number of QMC points to generate.
#'
#' @return An \eqn{n \times 2} numeric matrix. All values are in \eqn{[0, 1]}.
#'
#' @examples
#' pts <- generate_qmc_unit_square(128)
#' stopifnot(all(pts >= 0), all(pts <= 1))
#'
#' @export
generate_qmc_unit_square <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L || n != as.integer(n) || is.na(n)) {
    stop("`n` must be a positive integer.")
  }
  qrng::sobol(n = as.integer(n), d = 2L)
}


# ------------------------------------------------------------------------------

#' Map Unit Square Points to a Triangle
#'
#' Maps an \eqn{n \times 2} matrix of points from \eqn{[0,1]^2} into a
#' triangle defined by three vertex coordinates, using a uniform barycentric
#' transformation. Points outside the standard simplex (\eqn{u + v > 1}) are
#' reflected across the diagonal so that the mapping is area-preserving and
#' bijective.
#'
#' This function is called once per triangle inside
#' \code{integrate_one_sounding()}, with the same QMC points (from
#' [generate_qmc_unit_square()]) reused for each triangle.
#'
#' @param qmc_points An \eqn{n \times 2} numeric matrix of points in
#'   \eqn{[0,1]^2}.
#' @param triangle_coords A \eqn{3 \times 2} numeric matrix of triangle vertex
#'   coordinates \eqn{(x, y)}, one vertex per row.
#'
#' @return An \eqn{n \times 2} numeric matrix of mapped coordinates. All
#'   returned points lie inside (or on the boundary of) the triangle.
#'
#' @examples
#' tri    <- matrix(c(0,0, 1,0, 0,1), ncol = 2, byrow = TRUE)
#' pts    <- generate_qmc_unit_square(64)
#' mapped <- map_unit_square_to_triangle(pts, tri)
#' stopifnot(all(mapped >= 0))
#'
#' @export
map_unit_square_to_triangle <- function(qmc_points, triangle_coords) {
  if (!is.matrix(qmc_points) || ncol(qmc_points) != 2L) {
    stop("`qmc_points` must be an n x 2 numeric matrix.")
  }
  if (!is.matrix(triangle_coords) || !identical(dim(triangle_coords), c(3L, 2L))) {
    stop("`triangle_coords` must be a 3 x 2 numeric matrix.")
  }
  if (!is.numeric(qmc_points) || !is.numeric(triangle_coords)) {
    stop("`qmc_points` and `triangle_coords` must be numeric.")
  }

  u <- qmc_points[, 1L]
  v <- qmc_points[, 2L]

  # Reflect points outside the standard simplex across the diagonal u + v = 1
  reflect    <- (u + v > 1)
  u[reflect] <- 1 - u[reflect]
  v[reflect] <- 1 - v[reflect]

  A <- triangle_coords[1L, ]
  B <- triangle_coords[2L, ]
  C <- triangle_coords[3L, ]

  mapped_x <- (1 - u - v) * A[1L] + u * B[1L] + v * C[1L]
  mapped_y <- (1 - u - v) * A[2L] + u * B[2L] + v * C[2L]

  cbind(mapped_x, mapped_y)
}
