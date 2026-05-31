# fourier_prior.R
#
# Prior precision constructors for the Fourier / cosine basis, compatible with
# fastblm::fit_fastblm.
#
# fastblm convention
# ------------------
# Posterior precision: A'A + (1/phi) * Q
# Prior:              beta ~ N(0, phi * sigma2e * Q^{-1})
# Therefore:          phi = sigma2b / sigma2e  (larger phi = weaker prior)
#
# Two normalisation conventions (must match fourier_freq_grid norm argument)
# --------------------------------------------------------------------------
#
#   norm = "physical"  (default)
#     Basis is L2-normalised on [0,Lx]x[0,Ly]: integral phi_k^2 ds = 1.
#     Q diagonal entry k:
#       q_k = n_cells * (kappa2 + omega_sq_k)^alpha
#     where n_cells = (Lx/delta)*(Ly/delta), omega_sq_k in physical rad/m units.
#     phi = sigma2_M / sigma2e  where sigma2_M is the Matern field variance.
#
#   norm = "unit_square"  (matches SpatialBasis::make_matern_fourier_basis)
#     Basis is normalised on [0,1]^2: c_k = 1 (k=0) or sqrt(2) (k>0).
#     Q diagonal entry k:
#       q_k = TAU * (kappa2_norm + lambda_norm_k)^alpha
#     where kappa2_norm and lambda_norm are in [0,1]^2 units,
#     TAU = 1 / (sigma2_SAR * (2d/(rho*h))^2).
#     phi = 1  (TAU is already absorbed into Q).
#
#     Verified: q_k identical between this and SpatialBasis::compute_precision
#     for all K modes (ratio sd < 1e-6). See verify_integration_paths3.R.
#
# SAR -> Matern parameter matching
# ---------------------------------
# get_matern_parameters_from_SAR(rho, sar_precision, h, d=2):
#   kappa2_norm = 2d*(1-rho) / (rho*h^2)          [unit-square units]
#   TAU         = 1 / (sigma2_SAR * (2d/(rho*h))^2)
#
# where h = 1/n_cells_per_side (normalised pixel size on [0,1]^2).
#
# SAR-Fourier exact spectral matching (physical units, norm="physical")
# ---------------------------------------------------------------------
# Q_f[k] = n_cells * (1 - rho * lambda_W[k])^2
# lambda_W[k] = 0.5*(cos(omega_x[k]*delta) + cos(omega_y[k]*delta))
# phi_fourier = sar_phi_to_fourier_phi(phi_sar, rho, delta, freq_grid)
#
# Workflow (unit_square, recommended for matching SpatialBasis)
# -------------------------------------------------------------
#   h      <- 1 / n_cells_side
#   params <- get_matern_parameters_from_SAR(rho, sar_precision, h)
#   fg     <- fourier_freq_grid(J1, J2, domain_bbox, norm = "unit_square")
#   Q      <- fourier_matern_Q_sb(fg, kappa2_norm = params$kappa2,
#                                 TAU = params$matern_precision, alpha = 2L)
#   A      <- fourier_integrate_basis(polygons_sf, fg)
#   fit    <- fit_fastblm(y, A, Q$Q, phi = 1,
#                         solver = "woodbury", Q_inv = function(v) v / Q$q_diag)
#
# Workflow (physical, original fastblm workflow)
# -----------------------------------------------
#   fg      <- fourier_freq_grid(J1, J2, domain_bbox, norm = "physical")
#   Q       <- fourier_sar_Q(fg, rho, delta)          # or fourier_matern_Q()
#   phi     <- sar_phi_to_fourier_phi(phi_sar, rho, delta, fg)
#   q_diag  <- diag(Q)
#   fit     <- fit_fastblm(y, A, Q, phi,
#                          solver = "woodbury", Q_inv = function(v) v / q_diag)


# ------------------------------------------------------------------------------
# SAR -> Matern parameter matching (unit-square frame)
# ------------------------------------------------------------------------------

#' Convert SAR parameters to matched Matern parameters (unit-square frame)
#'
#' All quantities live in the normalised [0,1]^2 domain.
#' Use h = 1/n_cells_side (normalised pixel size).
#'
#' @param rho           SAR autocorrelation in (0,1).
#' @param sar_precision SAR prior precision tau (= 1/sigma2_SAR).
#' @param h             Normalised pixel size (1/n_cells_side).
#' @param d             Dimension (default 2).
#'
#' @return Named list: matern_precision (TAU), kappa2 (in [0,1]^2 units).
#'
#' @export
get_matern_parameters_from_SAR <- function(rho, sar_precision, h, d = 2) {
  if (!is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("`rho` must be in (0, 1).")
  if (!is.numeric(sar_precision) || sar_precision <= 0)
    stop("`sar_precision` must be positive.")
  if (!is.numeric(h) || h <= 0)
    stop("`h` must be positive.")

  sar_variance    <- 1 / sar_precision
  matern_variance <- sar_variance * (2 * d / (rho * h))^2
  kappa2          <- 2 * d * (1 - rho) / rho / h^2

  list(
    matern_precision = 1 / matern_variance,   # TAU
    kappa2           = kappa2                 # in [0,1]^2 units
  )
}


# ------------------------------------------------------------------------------
# Prior precision: unit-square convention (matches SpatialBasis)
# ------------------------------------------------------------------------------

#' Build Matern-Fourier prior precision matrix (unit-square normalisation)
#'
#' Constructs the diagonal prior Q matched to SpatialBasis::make_matern_fourier_prior.
#' Use with fourier_freq_grid(..., norm = "unit_square") and phi = 1 in fit_fastblm.
#'
#' Diagonal entries:
#'   q_k = TAU * (kappa2_norm + pi^2*(kx_k^2 + ky_k^2))^alpha
#'
#' where kx_k, ky_k are the 0-indexed frequency integers stored in freq_grid.
#'
#' Verified identical to SpatialBasis::compute_precision(make_matern_fourier_prior(...))
#' for all K modes (ratio range [1,1], sd < 1e-6). See verify_integration_paths3.R.
#'
#' @param freq_grid    List from fourier_freq_grid(..., norm = "unit_square").
#' @param kappa2_norm  kappa^2 in [0,1]^2 units. From get_matern_parameters_from_SAR().
#' @param TAU          Matern precision. From get_matern_parameters_from_SAR().
#' @param alpha        Smoothness integer (default 2).
#'
#' @return Named list: Q (diagonal sparse Matrix), q_diag (numeric vector).
#'
#' @seealso get_matern_parameters_from_SAR, fourier_freq_grid
#' @export
fourier_matern_Q_sb <- function(freq_grid, kappa2_norm, TAU, alpha = 2L) {
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (is.null(freq_grid$norm) || freq_grid$norm != "unit_square")
    warning("`freq_grid` was not built with norm='unit_square'; results may be incorrect.")
  if (!is.numeric(kappa2_norm) || kappa2_norm <= 0)
    stop("`kappa2_norm` must be positive.")
  if (!is.numeric(TAU) || TAU <= 0)
    stop("`TAU` must be positive.")

  # kx, ky: 0-indexed frequency integers recovered from indices (j - 1)
  kx <- freq_grid$indices[, 1L] - 1L
  ky <- freq_grid$indices[, 2L] - 1L

  lambda_norm <- pi^2 * (kx^2 + ky^2)
  q_diag      <- TAU * (kappa2_norm + lambda_norm)^alpha

  list(
    Q      = Matrix::Diagonal(x = q_diag),
    q_diag = q_diag
  )
}


# ------------------------------------------------------------------------------
# Parameter conversions (physical units)
# ------------------------------------------------------------------------------

#' Convert SAR phi to Fourier phi (physical units)
#'
#' phi_fourier = phi_sar * 16 * Lx * Ly / (rho^2 * delta^2)
#'
#' @param phi_sar  Positive numeric. phi from SAR CV fit.
#' @param rho      Numeric in (0,1).
#' @param delta    Positive numeric. Pixel size in CRS units.
#' @param freq_grid List from fourier_freq_grid().
#'
#' @return Positive numeric. phi for the Fourier model.
#' @export
sar_phi_to_fourier_phi <- function(phi_sar, rho, delta, freq_grid) {
  if (!is.numeric(phi_sar) || phi_sar <= 0)
    stop("`phi_sar` must be a positive scalar.")
  if (!is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("`rho` must be in (0, 1).")
  if (!is.numeric(delta) || delta <= 0)
    stop("`delta` must be a positive scalar.")
  if (!is.list(freq_grid) || is.null(freq_grid$domain_bbox))
    stop("`freq_grid` must be a list from fourier_freq_grid().")

  Lx <- freq_grid$domain_bbox[3L] - freq_grid$domain_bbox[1L]
  Ly <- freq_grid$domain_bbox[4L] - freq_grid$domain_bbox[2L]
  phi_sar * 16 * Lx * Ly / (rho^2 * delta^2)
}


#' Convert Fourier phi back to SAR phi
#' @export
fourier_phi_to_sar_phi <- function(phi_fourier, rho, delta, freq_grid) {
  if (!is.numeric(phi_fourier) || phi_fourier <= 0)
    stop("`phi_fourier` must be a positive scalar.")
  if (!is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("`rho` must be in (0, 1).")
  if (!is.numeric(delta) || delta <= 0)
    stop("`delta` must be a positive scalar.")
  if (!is.list(freq_grid) || is.null(freq_grid$domain_bbox))
    stop("`freq_grid` must be a list from fourier_freq_grid().")

  Lx <- freq_grid$domain_bbox[3L] - freq_grid$domain_bbox[1L]
  Ly <- freq_grid$domain_bbox[4L] - freq_grid$domain_bbox[2L]
  phi_fourier * rho^2 * delta^2 / (16 * Lx * Ly)
}


# ------------------------------------------------------------------------------
# Prior precision: physical convention (original)
# ------------------------------------------------------------------------------

#' Build matched SAR-Fourier prior precision matrix (physical normalisation)
#'
#' Uses exact discrete Laplacian eigenvalues. Use with norm="physical" freq_grid
#' and phi = sar_phi_to_fourier_phi(phi_sar, rho, delta, freq_grid).
#'
#' q_k = (1 - rho * lambda_W[k])^2
#' lambda_W[k] = 0.5 * (cos(omega_x[k]*delta) + cos(omega_y[k]*delta))
#'
#' @param freq_grid List from fourier_freq_grid().
#' @param rho       Numeric in (0,1).
#' @param delta     Positive numeric. Pixel size in CRS units.
#'
#' @return Named list: Q (diagonal sparse Matrix), q_diag (numeric vector).
#' @export
fourier_sar_Q <- function(freq_grid, rho, delta) {
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (!is.numeric(rho) || length(rho) != 1L || rho <= 0 || rho >= 1)
    stop("`rho` must be a scalar in (0, 1).")
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0)
    stop("`delta` must be a positive scalar.")

  omega_x  <- freq_grid$omega_mat[, 1L]
  omega_y  <- freq_grid$omega_mat[, 2L]
  lambda_W <- 0.5 * (cos(omega_x * delta) + cos(omega_y * delta))
  q_diag   <- (1 - rho * lambda_W)^2

  list(
    Q      = Matrix::Diagonal(x = q_diag),
    q_diag = q_diag
  )
}


#' Build Matern-Fourier prior precision matrix (physical normalisation)
#'
#' q_k = n_cells * (kappa2 + omega_sq_k)^alpha
#' Use with phi = sigma2_M / sigma2e.
#'
#' @param freq_grid List from fourier_freq_grid().
#' @param kappa2    kappa^2 in physical CRS units (1/m^2).
#' @param delta     Positive numeric. Pixel size in CRS units.
#' @param alpha     Smoothness integer (default 2).
#'
#' @return Named list: Q (diagonal sparse Matrix), q_diag (numeric vector).
#' @export
fourier_matern_Q <- function(freq_grid, kappa2, delta, alpha = 2L) {
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (!is.numeric(kappa2) || kappa2 <= 0)
    stop("`kappa2` must be positive.")
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0)
    stop("`delta` must be a positive scalar.")

  Lx      <- freq_grid$domain_bbox[3L] - freq_grid$domain_bbox[1L]
  Ly      <- freq_grid$domain_bbox[4L] - freq_grid$domain_bbox[2L]
  n_cells <- (Lx / delta) * (Ly / delta)

  omega_sq <- freq_grid$omega_mat[, 1L]^2 + freq_grid$omega_mat[, 2L]^2
  q_diag   <- n_cells * (kappa2 + omega_sq)^alpha

  list(
    Q      = Matrix::Diagonal(x = q_diag),
    q_diag = q_diag
  )
}
