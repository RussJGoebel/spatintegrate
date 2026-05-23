# fourier_prior.R
# Prior precision constructors for the Fourier basis, compatible with
# fastblm::tune_ml.
#
# The Fourier basis diagonalises any stationary covariance, so Q is always
# diagonal. Two prior families are supported:
#
#   Matern (nu = 1):  parameters (range, sigma2)
#   SAR (rook, nu=1): parameters (rho, tau2, delta)
#
# Both families are internally the same model — the SAR constructor converts
# to Matern parameters first via exact formulas derived from matching the
# discrete Laplacian (rook neighbors) to the SPDE (kappa^2 - Delta)u = W.
#
# Conversion formulas (nu = 1, d = 2):
#
# From the discrete Laplacian / SPDE matching (Chapter 5, Table 1):
#
#   kappa^2      = (2d/h^2) * (1-rho)/rho = 4*(1-rho)/(rho*delta^2)
#   sigma_SAR^2  = sigma_M^2 * h^2 / (kappa^2 + 4/h^2)^2
#                = sigma_M^2 * delta^2 / (kappa^2 + 4/delta^2)^2
#
# and the inverse:
#
#   range   = sqrt(8 / kappa^2)   [kappa^2 = 8/range^2 convention]
#   rho     = (4/delta^2) / (kappa^2 + 4/delta^2)
#   sigma_M^2 = sigma_SAR^2 * (kappa^2 + 4/delta^2)^2 / delta^2
#
# Note: sigma_SAR^2 here is the marginal variance of y (the grid vector),
# which differs from sigma_M^2 by a factor of delta^2 due to the
# vector vs function norm difference (see dissertation Chapter 5).


# ------------------------------------------------------------------------------
# Parameter conversions
# ------------------------------------------------------------------------------

#' Convert Matern Parameters to SAR Parameters
#'
#' Converts Matern(nu=1) parameters to equivalent SAR (rook neighbors) parameters
#' on a regular grid with spacing \code{delta}, using the exact correspondence
#' between the discrete Laplacian and the SPDE (kappa^2 - Delta)u = W.
#'
#' @param range Positive numeric. Matern range (CRS units).
#' @param sigma2 Positive numeric. Matern marginal variance.
#' @param delta Positive numeric. Grid spacing (CRS units).
#'
#' @return Named list with elements \code{rho} and \code{tau2}.
#' @seealso [sar_to_matern()]
#' @export
matern_to_sar <- function(range, sigma2, delta) {
  if (!is.numeric(range)  || range  <= 0) stop("`range` must be a positive scalar.")
  if (!is.numeric(sigma2) || sigma2 <= 0) stop("`sigma2` must be a positive scalar.")
  if (!is.numeric(delta)  || delta  <= 0) stop("`delta` must be a positive scalar.")
  kappa2   <- 8 / range^2
  rho      <- (4 / delta^2) / (kappa2 + 4 / delta^2)
  # sigma_SAR^2 = sigma_M^2 * delta^2 / (kappa^2 + 4/delta^2)^2
  sigma_sar2 <- sigma2 * delta^2 / (kappa2 + 4 / delta^2)^2
  list(rho = rho, sigma_sar2 = sigma_sar2)
}


#' Convert SAR Parameters to Matern Parameters
#'
#' Converts SAR (rook neighbors, nu=1) parameters to equivalent Matern(nu=1)
#' parameters, using the exact discrete-Laplacian correspondence.
#'
#' @param rho Numeric in (0, 1). SAR spatial autoregression parameter.
#' @param tau2 Positive numeric. SAR conditional variance.
#' @param delta Positive numeric. Grid spacing (CRS units).
#'
#' @return Named list with elements \code{range} and \code{sigma2}.
#' @seealso [matern_to_sar()]
#' @export
sar_to_matern <- function(rho, sigma_sar2, delta) {
  if (!is.numeric(rho)       || rho       <= 0 || rho       >= 1) stop("`rho` must be in (0, 1).")
  if (!is.numeric(sigma_sar2)|| sigma_sar2 <= 0)                  stop("`sigma_sar2` must be a positive scalar.")
  if (!is.numeric(delta)     || delta      <= 0)                  stop("`delta` must be a positive scalar.")
  kappa2   <- 4 * (1 - rho) / (rho * delta^2)
  range    <- sqrt(8 / kappa2)
  # sigma_M^2 = sigma_SAR^2 * (kappa^2 + 4/delta^2)^2 / delta^2
  sigma_m2 <- sigma_sar2 * (kappa2 + 4 / delta^2)^2 / delta^2
  list(range = range, sigma_m2 = sigma_m2)
}


# ------------------------------------------------------------------------------
# Matern eigenvalues for the Neumann cosine basis (dissertation Chapter 4)
#
# For basis index j = (j1, j2) on physical domain [0,Lx] x [0,Ly]:
#
#   lambda_j = sigma2 * (kappa^2 + omega_x^2 + omega_y^2)^{-alpha}
#
# where omega_k = (j_k - 1) * pi / L_k  is the Laplacian eigenvalue on [0,L_k],
# alpha = nu + d/2 = nu + 1  (d=2), and kappa^2 = 8 / range^2.
#
# This matches the tex formula lambda_j = sigma2*(kappa^2 + pi^2*sum(j_k-1)^2)^{-alpha}
# on the unit square [0,1]^2, since there omega_k = (j_k-1)*pi.
# On a physical domain [0,L]^2 the correct Laplacian eigenvalue is
# omega_k^2 = pi^2*(j_k-1)^2/L^2, already encoded in freq_grid$omega_mat.
#
# omega_sq: length-r vector ||omega_k||^2 from freq_grid$omega_mat
# ------------------------------------------------------------------------------

.matern_eigenvalues <- function(omega_sq, range, sigma2, nu) {
  kappa2 <- 8 / range^2
  alpha  <- nu + 1L          # nu + d/2, d = 2
  sigma2 * (kappa2 + omega_sq)^(-alpha)
}



# ------------------------------------------------------------------------------

#' Build a Q_fun for a Matern(nu=1) Prior over a Fourier Basis
#'
#' Returns a callable object compatible with fastblm::tune_ml() that constructs
#' the diagonal prior precision matrix Q for a Fourier basis under a
#' Matern(nu=1) covariance.
#'
#' Parameters tuned by tune_ml: log_range, log_sigma2.
#'
#' @param freq_grid A list from [fourier_freq_grid()].
#' @param nu Positive numeric. Matern smoothness. Fixed at construction time.
#'   Default 1 (matches SAR correspondence and dissertation Chapter 4).
#' @param dc_precision Positive numeric. Precision at the constant term
#'   (j1=j2=1). Default 1e-6 (near-flat prior on mean).
#'
#' @return Callable object of class fourier_matern_Q_fun.
#' @seealso [fourier_sar_Q_fun()], [summarize_prior()], [matern_to_sar()]
#' @export
fourier_matern_Q_fun <- function(freq_grid, nu = 1, dc_precision = 1e-6) {
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (!is.numeric(nu) || nu <= 0)
    stop("`nu` must be a positive numeric scalar.")
  if (!is.numeric(dc_precision) || dc_precision <= 0)
    stop("`dc_precision` must be a positive numeric scalar.")

  omega_sq <- freq_grid$omega_mat[, 1L]^2 + freq_grid$omega_mat[, 2L]^2
  is_dc    <- omega_sq == 0   # j1=j2=1: constant term

  fn <- function(theta) {
    range  <- exp(theta[["log_range"]])
    sigma2 <- exp(theta[["log_sigma2"]])
    eigs        <- .matern_eigenvalues(omega_sq, range, sigma2, nu)
    eigs[is_dc] <- 1 / dc_precision
    q_diag    <- 1 / eigs
    Q         <- Matrix::Diagonal(x = q_diag)
    log_det_Q <- sum(log(q_diag))
    list(Q = Q, log_det_Q = log_det_Q)
  }

  structure(fn,
            class        = "fourier_matern_Q_fun",
            freq_grid    = freq_grid,
            nu           = nu,
            dc_precision = dc_precision,
            param_names  = c("log_range", "log_sigma2")
  )
}


#' Build a Q_fun for a SAR(rook, nu=1) Prior over a Fourier Basis
#'
#' Returns a callable object compatible with fastblm::tune_ml() that constructs
#' the diagonal prior precision matrix Q for a Fourier basis under a SAR
#' (rook neighbors) prior. Internally converts SAR to Matern(nu=1) via exact
#' formulas, so fitted parameters can be reported in both parameterizations
#' via summarize_prior().
#'
#' Parameters tuned by tune_ml: logit_rho, log_tau2.
#'
#' @param freq_grid A list from [fourier_freq_grid()].
#' @param delta Positive numeric. Grid spacing (CRS units) of the SAR model
#'   being matched.
#' @param dc_precision Positive numeric. Precision at the constant term. Default 1e-6.
#'
#' @section Parameters tuned by tune_ml:
#' \describe{
#'   \item{\code{logit_rho}}{Logit of rho in (0,1). Controls spatial range.}
#'   \item{\code{log_tau2}}{Log of the SAR conditional variance tau2.
#'     Same tau2 as in Q_SAR = (I - rho*W)\'(I - rho*W) / tau2. Using the
#'     same (rho, tau2, delta) in both models gives exactly matched priors,
#'     since the Fourier basis diagonalises the SAR precision on a regular grid.}
#' }
#'
#' @return Callable object of class fourier_sar_Q_fun.
#' @seealso [fourier_matern_Q_fun()], [summarize_prior()], [sar_to_matern()]
#' @export
fourier_sar_Q_fun <- function(freq_grid, delta, dc_precision = NULL) {
  if (!is.list(freq_grid) || is.null(freq_grid$omega_mat))
    stop("`freq_grid` must be a list from fourier_freq_grid().")
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0)
    stop("`delta` must be a positive numeric scalar.")

  # Fourier frequencies for the Neumann cosine basis
  # omega_k = (j_k - 1) * pi / L_k
  # Eigenvalues of rook W at these frequencies:
  #   lambda_W(omega) = 0.5 * (cos(omega_x * delta) + cos(omega_y * delta))
  # This is exact for a periodic grid; for Neumann BC it is an approximation
  # that becomes exact in the interior as grid size grows.
  omega_x  <- freq_grid$omega_mat[, 1L]
  omega_y  <- freq_grid$omega_mat[, 2L]
  lambda_W <- 0.5 * (cos(omega_x * delta) + cos(omega_y * delta))

  # Number of SAR grid cells: n = (L/delta)^2
  # The Fourier basis is orthonormal on the continuous domain,
  # but SAR coefficients are field values at grid cells.
  # To match field variances: q_fourier_k = n * (1-rho*lambda_W_k)^2 / tau2
  Lx <- freq_grid$domain_bbox[3L] - freq_grid$domain_bbox[1L]
  Ly <- freq_grid$domain_bbox[4L] - freq_grid$domain_bbox[2L]
  n_cells <- (Lx / delta) * (Ly / delta)

  fn <- function(theta) {
    rho  <- 1 / (1 + exp(-theta[["logit_rho"]]))
    tau2 <- exp(theta[["log_tau2"]])
    # SAR precision eigenvalues scaled by n to match field variance:
    # Var(v(s)) under Fourier prior = Var(y_j) under SAR prior
    q_diag    <- n_cells * (1 - rho * lambda_W)^2 / tau2
    Q         <- Matrix::Diagonal(x = q_diag)
    log_det_Q <- sum(log(q_diag))
    list(Q = Q, log_det_Q = log_det_Q)
  }

  structure(fn,
            class       = "fourier_sar_Q_fun",
            freq_grid   = freq_grid,
            delta       = delta,
            param_names = c("logit_rho", "log_tau2")
  )
}


# ------------------------------------------------------------------------------
# Print and summarize
# ------------------------------------------------------------------------------

#' @export
print.fourier_matern_Q_fun <- function(x, ...) {
  freqs <- attr(x, "freqs")
  cat("Fourier Matern(nu=1) Q_fun\n")
  cat("  frequencies :", nrow(freqs), "(", 2L * nrow(freqs), "columns in A)\n")
  cat("  tune via    : log_range, log_sigma2\n")
  cat("  use summarize_prior(Q_fun, fit$theta) after fitting\n")
  invisible(x)
}

#' @export
print.fourier_sar_Q_fun <- function(x, ...) {
  fg <- attr(x, "freq_grid")
  cat("Fourier SAR(rook, nu=1) Q_fun\n")
  cat("  delta       :", attr(x, "delta"), "(CRS units)\n")
  cat("  basis fns   :", nrow(fg$omega_mat), "\n")
  cat("  tune via    : logit_rho, log_sigma2\n")
  cat("  use summarize_prior(Q_fun, fit$theta) after fitting\n")
  invisible(x)
}


#' Summarize a Fitted Fourier Prior
#'
#' Given a Q_fun and fitted theta from fastblm::tune_ml(), prints parameters
#' in both Matern(nu=1) and SAR(rook) parameterizations.
#'
#' @param Q_fun A fourier_matern_Q_fun or fourier_sar_Q_fun object.
#' @param theta Named numeric vector of fitted parameters from fit$theta.
#' @param delta Grid spacing for SAR equivalents. Required for
#'   fourier_matern_Q_fun; ignored for fourier_sar_Q_fun.
#'
#' @return Invisibly returns a list with elements matern and sar.
#' @export
summarize_prior <- function(Q_fun, theta, delta = NULL) {
  UseMethod("summarize_prior")
}

#' @export
summarize_prior.fourier_matern_Q_fun <- function(Q_fun, theta, delta = NULL) {
  range  <- exp(theta[["log_range"]])
  sigma2 <- exp(theta[["log_sigma2"]])

  cat("Fitted Fourier prior — Matern(nu=1)\n")
  cat("  --- Matern ---\n")
  cat(sprintf("  range  : %.1f  (CRS units)\n", range))
  cat(sprintf("  sigma2 : %.4g\n", sigma2))

  sar <- if (!is.null(delta)) matern_to_sar(range, sigma2, delta) else NULL

  if (!is.null(sar)) {
    cat(sprintf("  --- SAR equivalent (delta = %.1f) ---\n", delta))
    cat(sprintf("  rho    : %.4f\n", sar$rho))
    cat(sprintf("  tau2   : %.4g\n", sar$tau2))
  } else {
    cat("  (supply delta to see SAR equivalent)\n")
  }

  invisible(list(matern = list(range = range, sigma2 = sigma2), sar = sar))
}

#' @export
summarize_prior.fourier_sar_Q_fun <- function(Q_fun, theta, delta = NULL) {
  rho  <- 1 / (1 + exp(-theta[["logit_rho"]]))
  tau2 <- exp(theta[["log_tau2"]])
  d    <- attr(Q_fun, "delta")

  # Convert tau2 to sigma_sar2 for the Matern conversion
  # (tau2 and sigma_sar2 are the same parameter: both are the SAR marginal variance)
  matern <- sar_to_matern(rho, tau2, d)

  cat("Fitted Fourier prior — SAR(rook, nu=1)\n")
  cat(sprintf("  --- SAR (delta = %.1f) ---\n", d))
  cat(sprintf("  rho      : %.4f\n", rho))
  cat(sprintf("  tau2     : %.4g\n", tau2))
  cat("  --- Matern equivalent ---\n")
  cat(sprintf("  range    : %.1f  (CRS units)\n", matern$range))
  cat(sprintf("  sigma_M^2: %.4g\n", matern$sigma_m2))

  invisible(list(
    sar    = list(rho = rho, tau2 = tau2, delta = d),
    matern = list(range = matern$range, sigma_m2 = matern$sigma_m2)
  ))
}
