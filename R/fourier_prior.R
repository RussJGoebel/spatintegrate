# fourier_prior.R
#
# Prior precision constructors for the Fourier basis, compatible with
# fastblm::fit_fastblm and fastblm::tune_cv.
#
# fastblm convention
# ------------------
# Posterior precision: A'A + (1/phi) * Q
# Prior:              beta ~ N(0, phi * sigma2e * Q^{-1})
# Therefore:          phi = sigma2b / sigma2e  (larger phi = weaker prior)
#
# SAR-Fourier matching
# --------------------
# The SAR model (I - rho*W)y = e has precision Q_SAR = (I-rhoW)'(I-rhoW).
# Its eigenfunctions are the Neumann cosines; eigenvalue for mode k is:
#
#   (1 - rho * lambda_W[k])^2
#
# where lambda_W[k] = 0.5*(cos(omega_x[k]*delta) + cos(omega_y[k]*delta))
# is the exact discrete Laplacian eigenvalue at the physical frequency.
#
# The matched Fourier prior precision diagonal is:
#
#   Q_f[k] = n_cells * (1 - rho * lambda_W[k])^2
#
# where n_cells = (Lx/delta) * (Ly/delta) = m1*m2.
# The factor n_cells accounts for the physical L2 normalisation of the
# basis (integral phi_k^2 ds = 1) vs the discrete SAR convention.
#
# The same phi from the SAR fit transfers directly to the Fourier model
# after converting via sar_phi_to_fourier_phi():
#
#   phi_fourier = phi_sar * 16 * (Lx * Ly) / (rho^2 * delta^2)
#
# This converts the SAR signal-to-noise ratio to the equivalent Matern
# signal-to-noise ratio, accounting for the sigma_SAR -> sigma_M conversion.
# sigma2e cancels in the conversion so it is not needed.
#
# Workflow
# --------
#   # 1. Fit SAR model, get phi_hat
#   cv      <- tune_cv(y, A, Q_fun_sar, ...)
#   phi_sar <- cv$phi
#
#   # 2. Build frequency grid (J1=m2 cols, J2=m1 rows)
#   fg <- fourier_freq_grid(J1=m2, J2=m1, domain_bbox=bbox)
#
#   # 3. Build matched Q and convert phi
#   Q_f         <- fourier_sar_Q(fg, rho=rho_hat, delta=delta)
#   phi_fourier <- sar_phi_to_fourier_phi(phi_sar, rho=rho_hat, delta=delta)
#
#   # 4. Fit Fourier model
#   q_diag <- diag(Q_f)
#   fit_f  <- fit_fastblm(y, A_f, Q_f, phi=phi_fourier,
#                         solver="woodbury", Q_inv=function(v) v/q_diag)


# ------------------------------------------------------------------------------
# Parameter conversions
# ------------------------------------------------------------------------------

#' Convert SAR phi to Fourier phi
#'
#' Converts the signal-to-noise ratio phi from a SAR model to the equivalent
#' phi for the matched Fourier-Matern model. The conversion uses the
#' SAR-to-Matern amplitude relationship:
#'
#'   sigma2_M = sigma2_SAR * (4 / (rho * delta^2))^2
#'
#' Since phi = sigma2b / sigma2e and sigma2e cancels, this gives:
#'
#'   phi_fourier = phi_sar * 16 * (Lx * Ly) / (rho^2 * delta^2)
#'
#' @param phi_sar Positive numeric. phi from SAR CV fit.
#' @param rho     Numeric in (0,1). SAR spatial autoregression parameter.
#' @param delta   Positive numeric. Pixel size in CRS units.
#'
#' @return Positive numeric. phi for the Fourier model.
#'
#' @seealso fourier_sar_Q
#'
#' @examples
#' phi_fourier <- sar_phi_to_fourier_phi(phi_sar=185, rho=0.95, delta=330, freq_grid=fg)
#'
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
#'
#' Inverse of sar_phi_to_fourier_phi().
#'
#' @param phi_fourier Positive numeric. phi from Fourier model.
#' @param rho         Numeric in (0,1). SAR spatial autoregression parameter.
#' @param delta       Positive numeric. Pixel size in CRS units.
#'
#' @return Positive numeric. Equivalent SAR phi.
#'
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


#' Convert SAR Parameters to Matern Parameters
#'
#' @param rho      Numeric in (0,1). SAR spatial autoregression parameter.
#' @param sigma2_sar Positive numeric. SAR marginal field variance (= phi*sigma2e).
#' @param delta    Positive numeric. Pixel size in CRS units.
#'
#' @return Named list: kappa2, range, sigma2_M, eff_range.
#'
#' @export
sar_to_matern <- function(rho, sigma2_sar, delta) {
  if (!is.numeric(rho)       || rho       <= 0 || rho >= 1)
    stop("`rho` must be in (0, 1).")
  if (!is.numeric(sigma2_sar)|| sigma2_sar <= 0)
    stop("`sigma2_sar` must be a positive scalar.")
  if (!is.numeric(delta)     || delta <= 0)
    stop("`delta` must be a positive scalar.")

  kappa2   <- 4 * (1 - rho) / (rho * delta^2)
  # sigma_M^2 = sigma_SAR^2 * (kappa^2 + 4/delta^2)^2 / delta^2
  # The /delta^2 converts from discrete (vector) variance to
  # continuous (function) variance.
  sigma2_M <- sigma2_sar * (kappa2 + 4 / delta^2)^2 / delta^2
  list(
    kappa2    = kappa2,
    range     = sqrt(8 / kappa2),
    sigma2_M  = sigma2_M,
    eff_range = sqrt(8 / kappa2)
  )
}


#' Convert Matern Parameters to SAR Parameters
#'
#' Inverse of sar_to_matern().
#'
#' @param range    Positive numeric. Matern range (CRS units).
#' @param sigma2_M Positive numeric. Matern marginal variance.
#' @param delta    Positive numeric. Pixel size in CRS units.
#'
#' @return Named list: rho, sigma2_sar.
#'
#' @export
matern_to_sar <- function(range, sigma2_M, delta) {
  if (!is.numeric(range)   || range   <= 0) stop("`range` must be positive.")
  if (!is.numeric(sigma2_M)|| sigma2_M <= 0) stop("`sigma2_M` must be positive.")
  if (!is.numeric(delta)   || delta   <= 0) stop("`delta` must be positive.")

  kappa2     <- 8 / range^2
  rho        <- (4 / delta^2) / (kappa2 + 4 / delta^2)
  sigma2_sar <- sigma2_M * delta^2 / (kappa2 + 4 / delta^2)^2
  list(rho = rho, sigma2_sar = sigma2_sar)
}


# ------------------------------------------------------------------------------
# Prior precision constructors
# ------------------------------------------------------------------------------

#' Build Matched SAR-Fourier Prior Precision Matrix
#'
#' Constructs the diagonal prior precision matrix Q_f for the Fourier basis
#' matched to a SAR(rook) prior. Uses the exact discrete Laplacian eigenvalues
#' rather than the continuous Matern approximation, giving exact spectral
#' matching at all frequencies.
#'
#' The diagonal entries are:
#'
#'   Q_f[k] = n_cells * (1 - rho * lambda_W[k])^2
#'
#' where:
#'   n_cells    = (Lx/delta) * (Ly/delta)  =  m1 * m2
#'   lambda_W[k] = 0.5 * (cos(omega_x[k]*delta) + cos(omega_y[k]*delta))
#'
#' lambda_W[k] is the exact eigenvalue of the row-normalised rook weight
#' matrix W at the continuous frequency (omega_x[k], omega_y[k]).
#'
#' Use with phi = sar_phi_to_fourier_phi(phi_sar, rho, delta).
#'
#' @param freq_grid List from fourier_freq_grid(J1=m2, J2=m1, domain_bbox).
#' @param rho       Numeric in (0,1). SAR spatial autoregression parameter.
#' @param delta     Positive numeric. Pixel size in CRS units.
#'
#' @return A diagonal sparse Matrix of class ddiMatrix, dimension K x K.
#'
#' @seealso sar_phi_to_fourier_phi, fourier_freq_grid
#'
#' @examples
#' bbox <- c(0, 0, 170*330, 128*330)
#' fg   <- fourier_freq_grid(J1=170, J2=128, domain_bbox=bbox)
#' Q_f  <- fourier_sar_Q(fg, rho=0.95, delta=330)
#' # phi_f <- sar_phi_to_fourier_phi(phi_sar, rho=0.95, delta=330, freq_grid=fg)
#' # fit   <- fit_fastblm(y, A_f, Q_f, phi=phi_f, solver="woodbury",
#' #                      Q_inv=function(v) v/diag(Q_f))
#'
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

  # Exact discrete Laplacian eigenvalue of row-normalised rook W at frequency k
  lambda_W <- 0.5 * (cos(omega_x * delta) + cos(omega_y * delta))

  # Shape only -- no n_cells factor. Amplitude is carried by phi_fourier
  # from sar_phi_to_fourier_phi().
  q_diag <- (1 - rho * lambda_W)^2

  Matrix::Diagonal(x = q_diag)
}


#' Build Matern-Fourier Prior Precision Matrix
#'
#' Constructs the diagonal prior precision matrix Q_f for the Fourier basis
#' under a Matern(nu=1) prior. Uses the continuous spectral eigenvalues.
#'
#' The diagonal entries are:
#'
#'   Q_f[k] = n_cells * (kappa^2 + omega_sq[k])^2
#'
#' where n_cells = (Lx/delta)*(Ly/delta) and kappa^2 = 8/range^2.
#'
#' Use with phi = phi_M where phi_M = sigma2_M / sigma2e.
#' To convert from a SAR fit: phi_M = sar_phi_to_fourier_phi(phi_sar, rho, delta).
#'
#' @param freq_grid List from fourier_freq_grid().
#' @param range     Positive numeric. Matern range (CRS units).
#' @param delta     Positive numeric. Pixel size in CRS units.
#'
#' @return A diagonal sparse Matrix of class ddiMatrix, dimension K x K.
#'
#' @seealso fourier_sar_Q, sar_phi_to_fourier_phi
#'
#' @export
fourier_matern_Q <- function(freq_grid, range, delta) {
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (!is.numeric(range) || length(range) != 1L || range <= 0)
    stop("`range` must be a positive scalar.")
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0)
    stop("`delta` must be a positive scalar.")

  Lx      <- freq_grid$domain_bbox[3L] - freq_grid$domain_bbox[1L]
  Ly      <- freq_grid$domain_bbox[4L] - freq_grid$domain_bbox[2L]
  n_cells <- (Lx / delta) * (Ly / delta)

  kappa2   <- 8 / range^2
  omega_sq <- freq_grid$omega_mat[, 1L]^2 + freq_grid$omega_mat[, 2L]^2

  q_diag <- n_cells * (kappa2 + omega_sq)^2

  Matrix::Diagonal(x = q_diag)
}
