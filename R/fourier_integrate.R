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
# Normalisation
# -------------
# Two normalisation conventions are supported via the `norm` argument:
#
#   norm = "physical"  (default, original behaviour)
#     phi_{j1,j2}(x,y) = c_{j1} * c_{j2} * cos(...) * cos(...)
#     c_{jk} = 1/sqrt(Lk) if jk=1,  sqrt(2/Lk) if jk>1
#     => integral phi_k^2 ds = 1 on [0,Lx]x[0,Ly]
#
#   norm = "unit_square"  (matches SpatialBasis::make_matern_fourier_basis)
#     c_{jk} = 1 if jk=1,  sqrt(2) if jk>1
#     => orthonormal on [0,1]^2 unit square embedding
#     Use this when combining with make_Q_cosine(..., norm="unit_square")
#     or with the SpatialBasis Matern-Fourier prior.
#
# Patch notes
# -----------
# - get_triangle_coords() replaced with sf::st_coordinates() since the
#   former is an internal unexported function. Verified equivalent output.
# - triangulate_sf() min_area argument removed (not supported in all versions).
# - norm argument added to support unit_square convention.
# - parallel argument added; parallelises over polygons via future.apply,
#   matching the pattern used in compute_overlap_fractions.


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


# Extract triangle vertex coordinates from a triangulated sf geometry.
# Replaces the unexported get_triangle_coords() from spatintegrate.
# sf::st_coordinates() returns a matrix with columns X, Y (plus ring indices);
# rows 1:3 are the three vertices (row 4 closes the ring).
.get_triangle_coords <- function(tri_geom) {
  sf::st_coordinates(tri_geom)[1:3, 1:2, drop = FALSE]
}


# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------

#' Build a Neumann Cosine Frequency Grid
#'
#' Constructs the mode index table, angular frequencies, and normalisation
#' constants for the 2D Neumann cosine basis on [0,Lx] x [0,Ly].
#'
#' @param J1          Positive integer. Frequency count in x (columns).
#' @param J2          Positive integer. Frequency count in y (rows). Default J1.
#' @param domain_bbox Length-4 numeric (xmin, ymin, xmax, ymax).
#' @param norm        Character. "physical" (default) or "unit_square".
#'   See file header for details.
#'
#' @return A list with: omega_mat (K x 2), indices (K x 2), norm_const (K),
#'   J1, J2, domain_bbox, norm.
#'
#' @export
fourier_freq_grid <- function(J1, J2 = J1, domain_bbox,
                              norm = c("physical", "unit_square")) {
  norm <- match.arg(norm)

  if (!is.numeric(J1) || length(J1) != 1L || J1 < 1L || J1 != as.integer(J1))
    stop("`J1` must be a positive integer.")
  if (!is.numeric(J2) || length(J2) != 1L || J2 < 1L || J2 != as.integer(J2))
    stop("`J2` must be a positive integer.")
  if (!is.numeric(domain_bbox) || length(domain_bbox) != 4L || any(is.na(domain_bbox)))
    stop("`domain_bbox` must be a length-4 numeric vector (xmin, ymin, xmax, ymax).")

  J1 <- as.integer(J1)
  J2 <- as.integer(J2)
  Lx <- domain_bbox[3L] - domain_bbox[1L]
  Ly <- domain_bbox[4L] - domain_bbox[2L]
  if (Lx <= 0 || Ly <= 0)
    stop("`domain_bbox` must have positive width and height.")

  grid <- expand.grid(j1 = seq_len(J1), j2 = seq_len(J2))

  omega_mat <- cbind(
    omega_x = (grid$j1 - 1L) * pi / Lx,
    omega_y = (grid$j2 - 1L) * pi / Ly
  )

  if (norm == "physical") {
    # Physical L2 norm: integral phi_k^2 ds = 1 on [0,Lx]x[0,Ly]
    c1 <- ifelse(grid$j1 == 1L, 1 / sqrt(Lx), sqrt(2 / Lx))
    c2 <- ifelse(grid$j2 == 1L, 1 / sqrt(Ly), sqrt(2 / Ly))
  } else {
    # Unit-square norm: matches SpatialBasis::make_matern_fourier_basis
    # phi_k values are O(1); use with make_Q_cosine(..., norm="unit_square")
    c1 <- ifelse(grid$j1 == 1L, 1, sqrt(2))
    c2 <- ifelse(grid$j2 == 1L, 1, sqrt(2))
  }

  list(
    omega_mat   = omega_mat,
    indices     = as.matrix(grid),
    norm_const  = c1 * c2,
    J1          = J1,
    J2          = J2,
    domain_bbox = domain_bbox,
    norm        = norm
  )
}


#' Exact Cosine Basis Integration Over Polygons
#'
#' Computes the n x K design matrix A_f where entry (i,k) is the
#' area-normalised integral of basis function k over polygon i:
#'
#'   A_f[i,k] = (1/|D_i|) * integral_{D_i} phi_k(s) ds
#'
#' @param polygons_sf  sf object in a projected CRS. One row per polygon.
#' @param freq_grid    List from fourier_freq_grid().
#' @param parallel     Logical. If \code{TRUE}, parallelises over polygons using
#'   the current \code{future::plan()}. Set a plan before calling, e.g.
#'   \code{future::plan(future::multisession, workers = 4)}. Default
#'   \code{FALSE}.
#'
#' @return Numeric matrix n x K.
#'
#' @export
fourier_integrate_basis <- function(polygons_sf, freq_grid, parallel = FALSE) {

  if (!inherits(polygons_sf, "sf"))
    stop("`polygons_sf` must be an sf object.")
  if (nrow(polygons_sf) == 0L)
    stop("`polygons_sf` must have at least one row.")
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")

  omega_mat  <- freq_grid$omega_mat
  norm_const <- freq_grid$norm_const
  r          <- nrow(omega_mat)
  n_poly     <- nrow(polygons_sf)
  geom       <- sf::st_geometry(polygons_sf)

  origin_x <- freq_grid$domain_bbox[1L]
  origin_y <- freq_grid$domain_bbox[2L]

  omega_minus <- cbind( omega_mat[, 1L], -omega_mat[, 2L])
  omega_plus  <- omega_mat

  # --- Per-polygon worker -----------------------------------------------------

  worker <- function(i) {
    row <- rep(NA_real_, r)

    tris <- tryCatch(
      triangulate_sf(geom[i]),        # no min_area: not supported in all versions
      error = function(e) sf::st_sfc()
    )
    if (length(tris) == 0L) return(row)

    areas      <- as.numeric(sf::st_area(tris))
    total_area <- sum(areas)
    if (total_area == 0) return(row)

    int_minus <- complex(r)
    int_plus  <- complex(r)

    for (j in seq_along(tris)) {
      coords       <- .get_triangle_coords(tris[[j]])
      coords[, 1L] <- coords[, 1L] - origin_x
      coords[, 2L] <- coords[, 2L] - origin_y
      int_minus <- int_minus +
        .fourier_triangle_integral(coords, omega_minus, areas[j])
      int_plus  <- int_plus  +
        .fourier_triangle_integral(coords, omega_plus,  areas[j])
    }

    norm_const * 0.5 * (Re(int_minus) + Re(int_plus)) / total_area
  }

  # --- Parallelise or run sequentially ----------------------------------------

  if (parallel) {
    rows <- future.apply::future_lapply(
      seq_len(n_poly),
      worker,
      future.globals = list(
        worker        = worker,
        geom          = geom,
        r             = r,
        omega_minus   = omega_minus,
        omega_plus    = omega_plus,
        norm_const    = norm_const,
        origin_x      = origin_x,
        origin_y      = origin_y,
        .sinc_e                    = .sinc_e,
        .fourier_triangle_integral = .fourier_triangle_integral,
        .get_triangle_coords       = .get_triangle_coords
      ),
      future.packages = "sf",
      future.seed     = TRUE
    )
  } else {
    rows <- lapply(seq_len(n_poly), worker)
  }

  # --- Assemble and return ----------------------------------------------------

  A <- matrix(unlist(rows, use.names = FALSE), nrow = n_poly, ncol = r, byrow = TRUE)

  colnames(A) <- paste0("phi_", freq_grid$indices[, 1L], "_",
                        freq_grid$indices[, 2L])
  A
}
