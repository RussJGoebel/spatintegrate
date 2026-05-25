# fourier_integrate.R
# Exact integration of the Neumann cosine basis over polygons.
#
# Matches the basis from the dissertation (Chapter 4):
#
#   phi_j(s) = prod_{k=1}^2 c_{j_k} * cos((j_k - 1) * pi * s_k / L_k)
#
# where s in [0, Lx] x [0, Ly] (domain-local coordinates),
# j = (j1, j2) with j1 in {1, ..., J1}, j2 in {1, ..., J2},
# and the normalisation constants are:
#
#   c_{j_k} = 1/sqrt(L_k)    when j_k = 1  (constant in that dimension)
#   c_{j_k} = sqrt(2/L_k)    when j_k > 1  (non-constant)
#
# This makes the basis orthonormal under the physical L2 inner product on
# [0,Lx] x [0,Ly]:
#
#   integral phi_{j1,j2}(s)^2 ds = 1   for all (j1, j2)
#
# and equivalently satisfies h^2 * Phi'Phi = I at the discrete level
# (where h = pixel size and Phi is the m x K matrix of basis values at
# cell centres).
#
# This normalisation matches make_Q_cosine() in sar_matern.R, which also
# uses the physical L2 norm. The two must be consistent for the matched
# SAR-Matern prior to work correctly.
#
# The corresponding Matern eigenvalues are (alpha = nu + d/2 = 2 for nu=1, d=2):
#   lambda_{j1,j2} = sigma2_M_scale * Lx * Ly /
#                    (kappa^2 + pi^2*((j1-1)^2/Lx^2 + (j2-1)^2/Ly^2))^2
#
# For a square grid (J1 = J2 = J, Lx = Ly = L) the normalisation simplifies to
#   c_{j_k} = 1/sqrt(L) or sqrt(2/L).
#
# Mode ordering: expand.grid(j1 = 1:J1, j2 = 1:J2), so j1 varies fastest.
# This matches make_Q_cosine() in sar_matern.R, guaranteeing that column k
# of A from fourier_integrate_basis() corresponds to diagonal entry k of Q
# from make_Q_cosine().
#
# Integration strategy
# --------------------
# Each polygon is triangulated. For each triangle T and each basis function
# phi_j, the exact integral is obtained from:
#
#   int_T cos(omega . s) ds = Re( int_T e^{i omega . s} ds )
#
# where omega_k = (j_k - 1) * pi / L_k maps back to the original CRS units
# (L_k = domain width in CRS units). The complex triangle integral is:
#
#   int_T e^{i omega . s} ds = 2|T| e^{i theta1} F(a, b)
#
# with theta1 = omega.s1, a = omega.(s2-s1), b = omega.(s3-s1), and F(a,b)
# derived from the barycentric change of variables (see .fourier_triangle_integral).
#
# IMPORTANT: triangle coordinates must be expressed relative to the domain
# origin (domain_bbox[1], domain_bbox[2]) before being passed to
# .fourier_triangle_integral. The basis functions are defined on
# [0, Lx] x [0, Ly], not on raw CRS coordinates. fourier_integrate_basis()
# applies this shift automatically using freq_grid$domain_bbox.
#
# Two functions:
#   .fourier_triangle_integral()  — internal, one triangle, all frequencies
#   fourier_integrate_basis()     — user-facing, returns n x (J1*J2) matrix A


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Numerically stable sinc_e(x) = (e^{ix} - 1) / (ix), sinc_e(0) = 1.
# Uses Taylor series for |x| < 1e-6 to avoid cancellation.
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
# triangle_coords : 3 x 2 numeric matrix, vertices s1, s2, s3 in domain-local
#                   coordinates (i.e. already shifted so domain origin = (0,0))
# omega_mat       : r x 2 numeric matrix, angular frequency vectors (rad / CRS unit)
# area            : scalar, area of T in CRS units^2
#
# Returns complex vector of length r.
# Re() gives cos integrals; Im() gives sin integrals.
.fourier_triangle_integral <- function(triangle_coords, omega_mat, area) {

  s1  <- triangle_coords[1L, ]
  d21 <- triangle_coords[2L, ] - s1
  d31 <- triangle_coords[3L, ] - s1

  theta1 <- as.vector(omega_mat %*% s1)
  a      <- as.vector(omega_mat %*% d21)
  b      <- as.vector(omega_mat %*% d31)

  # F(a, b) from barycentric change-of-variables:
  #   general (b != 0): [ e^{ib} sinc_e(a-b) - sinc_e(a) ] / (ib)
  #   b ~ 0, a != 0:    sinc_e(a) - (e^{ia} - sinc_e(a)) / (ia)
  #   b ~ 0, a ~ 0:     1/2   (so 2|T| * 1/2 = |T|, the area)
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

#' Build a Neumann Cosine Frequency Grid
#'
#' Constructs the index matrix and corresponding angular frequency vectors for
#' the 2D Neumann cosine basis on a rectangular domain, matching the
#' eigenfunctions of the Laplacian under Neumann boundary conditions.
#'
#' The basis functions are orthonormal under the physical L2 inner product on
#' \eqn{[0,L_x] \times [0,L_y]}:
#' \deqn{
#'   \phi_{\mathbf{j}}(s)
#'   = c^{(1)}_{j_1} c^{(2)}_{j_2}
#'     \cos\!\big((j_1 - 1)\pi\, s_x / L_x\big)
#'     \cos\!\big((j_2 - 1)\pi\, s_y / L_y\big)
#' }
#' where the normalisation constants are
#' \eqn{c^{(k)}_1 = 1/\sqrt{L_k}} (constant mode) and
#' \eqn{c^{(k)}_j = \sqrt{2/L_k}} for \eqn{j \geq 2} (non-constant modes).
#'
#' This ensures \eqn{\int \phi_j(s)^2 \, ds = 1} for all \eqn{j}, and at the
#' discrete level \eqn{h^2 \Phi^\top \Phi = I} where \eqn{h} is the pixel size.
#' This normalisation matches \code{make_Q_cosine()} in \code{sar_matern.R}.
#'
#' For a square grid (\eqn{J_1 = J_2 = J}) this is equivalent to the
#' previous single-\code{J} interface. For a rectangular \eqn{m_1 \times m_2}
#' raster grid use \code{J1 = m1}, \code{J2 = m2} to get exactly
#' \eqn{m_1 m_2} modes with no wasted columns.
#'
#' @section Mode ordering:
#' Modes are enumerated by \code{expand.grid(j1 = 1:J1, j2 = 1:J2)},
#' so \eqn{j_1} varies fastest. This matches the ordering used by
#' \code{make_Q_cosine()} in \code{sar_matern.R}, guaranteeing that column
#' \eqn{k} of the design matrix \eqn{A} from \code{\link{fourier_integrate_basis}}
#' corresponds to diagonal entry \eqn{k} of \eqn{Q} from \code{make_Q_cosine()}.
#'
#' @param J1 Positive integer. Number of frequency components in the \eqn{x}
#'   direction. For a raster grid use \code{J1 = m1}.
#' @param J2 Positive integer. Number of frequency components in the \eqn{y}
#'   direction. For a raster grid use \code{J2 = m2}.
#'   Default \code{J1} (square grid).
#' @param domain_bbox Length-4 numeric vector
#'   \eqn{(x_{\min}, y_{\min}, x_{\max}, y_{\max})} in the projected CRS.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{omega_mat}}{A \eqn{J_1 J_2 \times 2} matrix of angular
#'       frequencies in CRS units (radians per CRS unit).}
#'     \item{\code{indices}}{A \eqn{J_1 J_2 \times 2} integer matrix of
#'       \eqn{(j_1, j_2)} index pairs, \eqn{j_1} varying fastest.}
#'     \item{\code{norm_const}}{Length-\eqn{J_1 J_2} vector of physical-unit
#'       normalisation constants \eqn{c^{(1)}_{j_1} \cdot c^{(2)}_{j_2}}.}
#'     \item{\code{J1}, \code{J2}}{The frequency counts used.}
#'     \item{\code{domain_bbox}}{The supplied bounding box.}
#'   }
#'
#' @seealso [fourier_integrate_basis()], \code{make_Q_cosine()} in sar_matern.R
#'
#' @examples
#' # Square grid
#' bbox <- c(0, 0, 200e3, 200e3)
#' fg   <- fourier_freq_grid(J1 = 5, domain_bbox = bbox)
#' fg$indices    # 25 x 2 integer matrix
#' fg$omega_mat  # 25 x 2 angular frequency matrix
#'
#' # Rectangular raster grid: 89 rows x 150 cols, h = 330m
#' bbox_rect <- c(0, 0, 89 * 330, 150 * 330)
#' fg_rect   <- fourier_freq_grid(J1 = 89, J2 = 150, domain_bbox = bbox_rect)
#' nrow(fg_rect$omega_mat)  # 89 * 150 = 13350
#'
#' @export
fourier_freq_grid <- function(J1, J2 = J1, domain_bbox) {

  # ── input validation ───────────────────────────────────────────────────────
  if (!is.numeric(J1) || length(J1) != 1L || J1 < 1L || J1 != as.integer(J1))
    stop("`J1` must be a positive integer.")
  if (!is.numeric(J2) || length(J2) != 1L || J2 < 1L || J2 != as.integer(J2))
    stop("`J2` must be a positive integer.")
  if (!is.numeric(domain_bbox) || length(domain_bbox) != 4L ||
      any(is.na(domain_bbox)))
    stop("`domain_bbox` must be a length-4 numeric vector (xmin, ymin, xmax, ymax).")

  J1 <- as.integer(J1)
  J2 <- as.integer(J2)
  Lx <- domain_bbox[3L] - domain_bbox[1L]
  Ly <- domain_bbox[4L] - domain_bbox[2L]
  if (Lx <= 0 || Ly <= 0)
    stop("`domain_bbox` must have positive width and height.")

  # ── index grid: j1 in {1,..,J1}, j2 in {1,..,J2}, j1 varies fastest ───────
  grid <- expand.grid(j1 = seq_len(J1), j2 = seq_len(J2))

  # ── angular frequencies: omega_k = (j_k - 1) * pi / L_k ──────────────────
  omega_mat <- cbind(
    omega_x = (grid$j1 - 1L) * pi / Lx,
    omega_y = (grid$j2 - 1L) * pi / Ly
  )

  # ── normalisation constants (physical L2 norm on [0,Lx] x [0,Ly]) ─────────
  # c_1 = 1/sqrt(L)   for the constant mode (j=1)
  # c_j = sqrt(2/L)   for non-constant modes (j>1)
  # Ensures: integral phi_j(s)^2 ds = 1  and  h^2 * Phi'Phi = I
  c1 <- ifelse(grid$j1 == 1L, 1 / sqrt(Lx), sqrt(2 / Lx))
  c2 <- ifelse(grid$j2 == 1L, 1 / sqrt(Ly), sqrt(2 / Ly))
  norm_const <- c1 * c2

  list(
    omega_mat   = omega_mat,
    indices     = as.matrix(grid),
    norm_const  = norm_const,
    J1          = J1,
    J2          = J2,
    domain_bbox = domain_bbox
  )
}


# ------------------------------------------------------------------------------

#' Exact Cosine Basis Integration Over a Set of Polygons
#'
#' Computes an \eqn{n \times J_1 J_2} matrix \eqn{A} where each entry is the
#' normalised area-weighted average of one Neumann cosine basis function over
#' one polygon:
#'
#' \deqn{
#'   A_{i,\mathbf{j}}
#'   = \frac{1}{|D_i|}
#'     \int_{D_i}
#'     \phi_{\mathbf{j}}(s)\, ds
#'   = \frac{c^{(1)}_{j_1} c^{(2)}_{j_2}}{|D_i|}
#'     \int_{D_i}
#'     \cos\!\big((j_1-1)\pi s_x / L_x\big)\,
#'     \cos\!\big((j_2-1)\pi s_y / L_y\big)\, ds
#' }
#'
#' where \eqn{c^{(k)}_{j_k}} are the physical-unit normalisation constants from
#' \code{\link{fourier_freq_grid}}, and \eqn{s_x, s_y} are coordinates relative
#' to the domain origin \eqn{(x_{\min}, y_{\min})} from
#' \code{freq_grid$domain_bbox}.
#'
#' Integration is exact up to floating-point precision: each polygon is
#' triangulated and the closed-form integral of \eqn{e^{i\omega\cdot s}} over
#' each triangle is accumulated analytically. No quadrature error is introduced.
#'
#' @section Coordinate shifting:
#' Triangle coordinates are shifted to domain-local frame
#' \eqn{[0, L_x] \times [0, L_y]} before integration, using the origin stored
#' in \code{freq_grid$domain_bbox}. This means \code{polygons_sf} may be in
#' any absolute projected CRS — only the bounding box passed to
#' \code{fourier_freq_grid()} needs to match the extent of the data.
#'
#' @section Note on separability:
#' The 2D cosine basis is separable:
#' \eqn{\phi_{j_1,j_2}(s) = \phi_{j_1}(s_x)\,\phi_{j_2}(s_y)}.
#' The integral over an arbitrary polygon is \emph{not} separable (the
#' polygon shape couples the two dimensions), so integration is performed
#' jointly over both dimensions using triangle decomposition.
#'
#' @param polygons_sf An \code{sf} object in a projected CRS. One row per
#'   polygon.
#' @param freq_grid A list returned by [fourier_freq_grid()].
#' @param min_area Non-negative numeric. Triangles with area \code{<= min_area}
#'   are skipped. Default \code{0}.
#'
#' @return A numeric matrix of dimensions \eqn{n \times J_1 J_2}. Columns
#'   correspond to \eqn{(j_1, j_2)} index pairs in the order returned by
#'   \code{freq_grid$indices} (\eqn{j_1} varying fastest). Rows for polygons
#'   producing no valid triangles contain \code{NA}.
#'
#' @seealso [fourier_freq_grid()]
#'
#' @examples
#' \dontrun{
#' bbox <- as.vector(sf::st_bbox(polygons_sf))
#' fg   <- fourier_freq_grid(J1 = 89, J2 = 150, domain_bbox = bbox)
#' A    <- fourier_integrate_basis(polygons_sf, fg)
#' dim(A)  # nrow(polygons_sf) x 13350
#' }
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

  # Domain origin — coordinates are shifted to [0, Lx] x [0, Ly] before
  # integration so that the basis functions cos((j-1)*pi*s/L) are evaluated
  # at the correct phase. Without this shift, raw CRS coordinates (e.g. UTM
  # values in the hundreds of thousands) produce cosines at completely wrong
  # phases, causing spatial artifacts in the fitted field.
  origin_x <- freq_grid$domain_bbox[1L]
  origin_y <- freq_grid$domain_bbox[2L]

  # For the 2D Neumann cosine basis, each basis function is a product:
  #   cos(omega_x * s_x) * cos(omega_y * s_y)
  #
  # This is NOT the same as Re(e^{i(omega_x*s_x + omega_y*s_y)}).
  #
  # The trig identity gives:
  #   cos(A) cos(B) = 0.5 * [cos(A - B) + cos(A + B)]
  #                 = 0.5 * Re(e^{i(omega_x*s_x - omega_y*s_y)}) +
  #                   0.5 * Re(e^{i(omega_x*s_x + omega_y*s_y)})
  #
  # So we call .fourier_triangle_integral twice per triangle:
  #   once with omega_minus = (omega_x, -omega_y)
  #   once with omega_plus  = (omega_x, +omega_y)
  # and average the real parts.
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
      coords <- get_triangle_coords(tris[[j]])

      # Shift to domain-local frame [0, Lx] x [0, Ly].
      # Areas are translation-invariant so areas[j] is unaffected.
      coords[, 1L] <- coords[, 1L] - origin_x
      coords[, 2L] <- coords[, 2L] - origin_y

      int_minus <- int_minus + .fourier_triangle_integral(coords, omega_minus, areas[j])
      int_plus  <- int_plus  + .fourier_triangle_integral(coords, omega_plus,  areas[j])
    }

    # Area-normalise, apply trig identity average, apply normalisation constant
    A[i, ] <- norm_const * 0.5 * (Re(int_minus) + Re(int_plus)) / total_area
  }

  # Column names: phi_{j1}_{j2}
  colnames(A) <- paste0("phi_", freq_grid$indices[, 1], "_", freq_grid$indices[, 2])
  A
}
