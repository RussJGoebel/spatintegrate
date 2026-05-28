# fourier_integrate.R
#
# Exact integration of the Neumann cosine basis over polygons.
#
# Coordinate convention
# ---------------------
# The basis is defined on the domain-local frame [0, Lx] x [0, Ly] where:
#
#   x = column direction,  Lx = domain_bbox[3] - domain_bbox[1]
#   y = row direction,     Ly = domain_bbox[4] - domain_bbox[2]
#
# For a raster grid with m1 rows and m2 columns and pixel size delta:
#
#   Lx = m2 * delta,   Ly = m1 * delta
#
# Pixel index j (row-major, j = 1..m1*m2):
#   row r = (j-1) %/% m2 + 1,   col c = (j-1) %% m2 + 1
#   pixel centre: x = (c - 0.5) * delta,   y = (r - 0.5) * delta
#
# Mode ordering
# -------------
# fourier_freq_grid(J1, J2, domain_bbox) enumerates modes as
#   expand.grid(j1 = 1:J1, j2 = 1:J2),  j1 varies fastest.
# j1 indexes x-frequencies (Lx), j2 indexes y-frequencies (Ly).
#
# For a full raster basis call with J1 = m2, J2 = m1:
#   column k of A from fourier_integrate_basis() corresponds to
#   diagonal entry k of Q from fourier_sar_Q() in fourier_prior.R.
#
# Normalisation
# -------------
# Basis functions are orthonormal under the continuous L2 inner product:
#
#   phi_{j1,j2}(x,y) = c_{j1} * c_{j2} *
#                      cos((j1-1)*pi*x/Lx) * cos((j2-1)*pi*y/Ly)
#
# where:
#   c_{jk} = 1/sqrt(Lk)    if jk = 1  (constant mode)
#   c_{jk} = sqrt(2/Lk)    if jk > 1  (non-constant mode)
#
# This ensures integral phi_{j1,j2}(s)^2 ds = 1 over [0,Lx]x[0,Ly].
#
# Integration strategy
# --------------------
# Each polygon is triangulated. For each triangle T and each basis function,
# the exact integral uses:
#
#   int_T cos(omega_x*sx) * cos(omega_y*sy) ds
#     = 0.5 * Re(int_T e^{i(omega_x*sx - omega_y*sy)} ds)
#     + 0.5 * Re(int_T e^{i(omega_x*sx + omega_y*sy)} ds)
#
# Each complex triangle integral is computed via .fourier_triangle_integral().


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Numerically stable sinc_e(x) = (e^{ix} - 1) / (ix),  sinc_e(0) = 1.
.sinc_e <- function(x) {
  EPS    <- 1e-6
  result <- complex(length(x))
  large  <- abs(x) > EPS
  if (any(large)) {
    xl            <- x[large]
    result[large] <- (exp(1i * xl) - 1) / (1i * xl)
  }
  if (any(!large)) {
    xs             <- x[!large]
    ixs            <- 1i * xs
    result[!large] <- 1 + ixs/2 + ixs^2/6 + ixs^3/24 + ixs^4/120
  }
  result
}


# Exact integral of e^{i omega . s} over triangle T.
#
# triangle_coords : 3 x 2 matrix, vertices in domain-local coords (origin = 0)
# omega_mat       : r x 2 matrix, angular frequency vectors (rad / CRS unit)
# area            : scalar, triangle area in CRS units^2
#
# Returns complex vector of length r.
.fourier_triangle_integral <- function(triangle_coords, omega_mat, area) {

  s1  <- triangle_coords[1L, ]
  d21 <- triangle_coords[2L, ] - s1
  d31 <- triangle_coords[3L, ] - s1

  theta1 <- as.vector(omega_mat %*% s1)
  a      <- as.vector(omega_mat %*% d21)
  b      <- as.vector(omega_mat %*% d31)

  EPS     <- 1e-6
  b_small <- abs(b) < EPS
  a_small <- abs(a) < EPS
  result  <- complex(length(a))

  idx <- !b_small
  if (any(idx)) {
    ag          <- a[idx]; bg <- b[idx]
    result[idx] <- (exp(1i * bg) * .sinc_e(ag - bg) - .sinc_e(ag)) / (1i * bg)
  }
  idx <- b_small & !a_small
  if (any(idx)) {
    as_         <- a[idx]
    result[idx] <- .sinc_e(as_) - (exp(1i * as_) - .sinc_e(as_)) / (1i * as_)
  }
  result[b_small & a_small] <- 0.5 + 0i

  2 * area * exp(1i * theta1) * result
}


# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------

#' Build a Neumann Cosine Frequency Grid
#'
#' Constructs the mode index table, angular frequencies, and normalisation
#' constants for the 2D Neumann cosine basis on [0,Lx] x [0,Ly].
#'
#' Coordinate convention: j1 indexes the x-direction (columns,
#' Lx = bbox[3]-bbox[1]); j2 indexes the y-direction (rows,
#' Ly = bbox[4]-bbox[2]). For a raster grid with m1 rows, m2 columns,
#' pixel size delta, call fourier_freq_grid(J1 = m2, J2 = m1, domain_bbox).
#'
#' Mode ordering: expand.grid(j1 = 1:J1, j2 = 1:J2), j1 varies fastest.
#' Column k of A from fourier_integrate_basis() corresponds to diagonal
#' entry k of Q from fourier_sar_Q() in fourier_prior.R.
#'
#' Normalisation: integral phi_{j1,j2}(s)^2 ds = 1 on [0,Lx]x[0,Ly].
#'
#' @param J1          Positive integer. Frequency count in x (columns).
#'   For full raster basis use J1 = m2.
#' @param J2          Positive integer. Frequency count in y (rows).
#'   Default J1. For full raster basis use J2 = m1.
#' @param domain_bbox Length-4 numeric (xmin, ymin, xmax, ymax).
#'
#' @return A list with: omega_mat (K x 2), indices (K x 2), norm_const (K),
#'   J1, J2, domain_bbox.
#'
#' @seealso fourier_integrate_basis, fourier_sar_Q (fourier_prior.R)
#'
#' @examples
#' # m1=128 rows, m2=170 cols, delta=330m
#' bbox <- c(304971, 4666232, 304971 + 170*330, 4666232 + 128*330)
#' fg   <- fourier_freq_grid(J1 = 170, J2 = 128, domain_bbox = bbox)
#' nrow(fg$omega_mat)  # 170 * 128 = 21760
#'
#' @export
fourier_freq_grid <- function(J1, J2 = J1, domain_bbox) {

  if (!is.numeric(J1) || length(J1) != 1L || J1 < 1L || J1 != as.integer(J1))
    stop("`J1` must be a positive integer.")
  if (!is.numeric(J2) || length(J2) != 1L || J2 < 1L || J2 != as.integer(J2))
    stop("`J2` must be a positive integer.")
  if (!is.numeric(domain_bbox) || length(domain_bbox) != 4L || any(is.na(domain_bbox)))
    stop("`domain_bbox` must be a length-4 numeric vector (xmin, ymin, xmax, ymax).")

  J1 <- as.integer(J1)
  J2 <- as.integer(J2)
  Lx <- domain_bbox[3L] - domain_bbox[1L]   # x = column direction
  Ly <- domain_bbox[4L] - domain_bbox[2L]   # y = row direction
  if (Lx <= 0 || Ly <= 0)
    stop("`domain_bbox` must have positive width and height.")

  # j1 indexes x (columns), j2 indexes y (rows), j1 varies fastest
  grid <- expand.grid(j1 = seq_len(J1), j2 = seq_len(J2))

  omega_mat <- cbind(
    omega_x = (grid$j1 - 1L) * pi / Lx,
    omega_y = (grid$j2 - 1L) * pi / Ly
  )

  # Physical L2 normalisation: integral phi_k^2 ds = 1 on [0,Lx]x[0,Ly]
  c1 <- ifelse(grid$j1 == 1L, 1 / sqrt(Lx), sqrt(2 / Lx))
  c2 <- ifelse(grid$j2 == 1L, 1 / sqrt(Ly), sqrt(2 / Ly))

  list(
    omega_mat   = omega_mat,
    indices     = as.matrix(grid),
    norm_const  = c1 * c2,
    J1          = J1,
    J2          = J2,
    domain_bbox = domain_bbox
  )
}


#' Exact Cosine Basis Integration Over Polygons
#'
#' Computes the n x K design matrix A_f where entry (i,k) is the
#' area-normalised integral of basis function k over polygon i:
#'
#'   A_f[i,k] = (1/|D_i|) * integral_{D_i} phi_k(s) ds
#'
#' Integration is exact: each polygon is triangulated and the closed-form
#' triangle integral of e^{i*omega.s} is accumulated analytically.
#'
#' Polygon coordinates are shifted to domain-local frame [0,Lx]x[0,Ly]
#' using freq_grid$domain_bbox before integration.
#'
#' @param polygons_sf  sf object in a projected CRS. One row per polygon.
#' @param freq_grid    List from fourier_freq_grid().
#' @param min_area     Non-negative numeric. Triangles with area <= min_area
#'   are skipped. Default 0.
#'
#' @return Numeric matrix n x K. Columns in j1-fastest order matching
#'   freq_grid$indices. Rows for polygons with no valid triangles are NA.
#'
#' @seealso fourier_freq_grid
#'
#' @export
fourier_integrate_basis <- function(polygons_sf, freq_grid, min_area = 0) {

  if (!inherits(polygons_sf, "sf"))
    stop("`polygons_sf` must be an sf object.")
  if (nrow(polygons_sf) == 0L)
    stop("`polygons_sf` must have at least one row.")
  assert_projected(polygons_sf, arg_name = "polygons_sf")

  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (!is.numeric(min_area) || length(min_area) != 1L ||
      is.na(min_area) || min_area < 0)
    stop("`min_area` must be a non-negative numeric scalar.")

  omega_mat  <- freq_grid$omega_mat
  norm_const <- freq_grid$norm_const
  r          <- nrow(omega_mat)
  n_poly     <- nrow(polygons_sf)
  geom       <- sf::st_geometry(polygons_sf)

  origin_x <- freq_grid$domain_bbox[1L]
  origin_y <- freq_grid$domain_bbox[2L]

  # cos(omega_x*sx)*cos(omega_y*sy)
  #   = 0.5*Re(e^{i(omega_x*sx - omega_y*sy)})
  #   + 0.5*Re(e^{i(omega_x*sx + omega_y*sy)})
  omega_minus <- cbind( omega_mat[, 1L], -omega_mat[, 2L])
  omega_plus  <- omega_mat

  A <- matrix(NA_real_, nrow = n_poly, ncol = r)

  for (i in seq_len(n_poly)) {

    tris <- tryCatch(
      triangulate_sf(geom[i], min_area = min_area),
      error = function(e) sf::st_sfc()
    )
    if (length(tris) == 0L) next

    areas      <- as.numeric(sf::st_area(tris))
    total_area <- sum(areas)
    if (total_area == 0) next

    int_minus <- complex(r)
    int_plus  <- complex(r)

    for (j in seq_along(tris)) {
      coords       <- get_triangle_coords(tris[[j]])
      coords[, 1L] <- coords[, 1L] - origin_x
      coords[, 2L] <- coords[, 2L] - origin_y
      int_minus <- int_minus +
        .fourier_triangle_integral(coords, omega_minus, areas[j])
      int_plus  <- int_plus  +
        .fourier_triangle_integral(coords, omega_plus,  areas[j])
    }

    A[i, ] <- norm_const * 0.5 * (Re(int_minus) + Re(int_plus)) / total_area
  }

  colnames(A) <- paste0("phi_", freq_grid$indices[, 1L], "_",
                        freq_grid$indices[, 2L])
  A
}
