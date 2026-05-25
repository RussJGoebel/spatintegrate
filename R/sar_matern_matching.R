# sar_matern.R
#
# SAR <-> Matérn parameter matching (Table 2.1) and cosine prior builder.
#
# Three independent concerns:
#   1. Parameter conversion : sar_to_matern() / matern_to_sar()
#   2. Prior precision      : make_Q_cosine()
#   3. Design matrix        : fourier_integrate_basis() in spatintegrate (separate)
#
# Matching formulas (rook adjacency, d=2, square pixels of side h):
#
#   sar_to_matern:
#     kappa2         = 4*(1-rho) / (rho * h^2)
#     sigma2_M_scale = sigma2_SAR_scale * (kappa2 + 4/h^2)^2
#                    = sigma2_SAR_scale * (4/(rho*h^2))^2
#
#   matern_to_sar:
#     rho             = (4/h^2) / (kappa2 + 4/h^2)
#     sigma2_SAR_scale = sigma2_M_scale / (kappa2 + 4/h^2)^2
#
# sigma2_SAR_scale and sigma2_M_scale are the SPDE/covariance shape parameters,
# independent of phi and sigma2e. In the package:
#   sigma2_SAR_scale = phi * sigma2e / tau   (= sigma2b / tau)
#
# The cosine basis phi_{j1,j2} is normalised so that h^2 * Phi'Phi = I
# (i.e. the physical L2 inner product). This means the prior on cosine
# coefficients beta needs a compensating factor of Lx*Ly = m1*m2*h^2 relative
# to the Table 2.1 sigma2_M_scale. make_Q_cosine handles this internally.
#
# h can be supplied directly or derived from domain_bbox + m1 (+ m2).
# Units of h are arbitrary but must match domain_bbox.
# Pixels must be square: Lx/m1 == Ly/m2.

# ── internal helper ───────────────────────────────────────────────────────────

#' Resolve pixel size h from explicit value or domain_bbox + grid dims
#' @keywords internal
.resolve_h <- function(h, domain_bbox, m1, m2) {
  if (!is.null(h)) {
    if (!is.numeric(h) || length(h) != 1L || h <= 0)
      stop("`h` must be a positive scalar.")
    if (!is.null(domain_bbox) && !is.null(m1)) {
      h_check <- (domain_bbox[3L] - domain_bbox[1L]) / m1
      if (abs(h - h_check) / h > 1e-4)
        warning(sprintf(
          "`h` (%.6g) is inconsistent with `domain_bbox` and `m1` (implied h=%.6g).",
          h, h_check))
    }
    return(h)
  }
  if (!is.null(domain_bbox) && !is.null(m1)) {
    if (!is.numeric(domain_bbox) || length(domain_bbox) != 4L)
      stop("`domain_bbox` must be a length-4 numeric vector (xmin,ymin,xmax,ymax).")
    if (!is.numeric(m1) || length(m1) != 1L || m1 < 2L)
      stop("`m1` must be an integer >= 2.")
    Lx <- domain_bbox[3L] - domain_bbox[1L]
    Ly <- domain_bbox[4L] - domain_bbox[2L]
    if (Lx <= 0 || Ly <= 0)
      stop("`domain_bbox` must have positive width and height.")
    hx <- Lx / m1
    if (!is.null(m2)) {
      hy <- Ly / m2
      if (abs(hx - hy) / hx > 1e-3)
        warning(sprintf(
          "Pixels are not square: Lx/m1=%.6g, Ly/m2=%.6g. SAR matching assumes square pixels.",
          hx, hy))
    }
    return(hx)
  }
  stop("Supply either `h` or both `domain_bbox` and `m1`.")
}

# ── 1. SAR -> Matern ──────────────────────────────────────────────────────────

#' Convert SAR parameters to matched Matérn parameters
#'
#' Uses the spectral matching from Table 2.1 to convert SAR prior shape
#' parameters to Matérn SPDE parameters.
#'
#' \code{sigma2_SAR_scale} is the prior covariance shape, independent of
#' \eqn{\phi} and \eqn{\sigma^2_e}. In the package convention:
#' \code{sigma2_SAR_scale = phi * sigma2e / tau}.
#'
#' @param rho              SAR autocorrelation \eqn{\rho \in (0,1)}.
#' @param sigma2_SAR_scale SAR prior covariance scale (Table 2.1 \eqn{\sigma^2_{\rm SAR}}).
#' @param h                Pixel side length in CRS units. Supply or derive via
#'   \code{domain_bbox} + \code{m1}.
#' @param domain_bbox      Length-4 numeric \eqn{(x_{\min},y_{\min},x_{\max},y_{\max})}.
#' @param m1               Number of grid rows.
#' @param m2               Number of grid cols (for square-pixel check).
#'
#' @return Named list: \code{kappa2}, \code{sigma2_M_scale}, \code{eff_range}, \code{h}.
#'
#' @seealso \code{\link{matern_to_sar}}, \code{\link{make_Q_cosine}}
#' @export
sar_to_matern <- function(rho, sigma2_SAR_scale,
                          h = NULL, domain_bbox = NULL, m1 = NULL, m2 = NULL) {
  if (rho <= 0 || rho >= 1)
    stop("`rho` must be in (0, 1).")
  if (!is.numeric(sigma2_SAR_scale) || length(sigma2_SAR_scale) != 1L ||
      sigma2_SAR_scale <= 0)
    stop("`sigma2_SAR_scale` must be a positive scalar.")

  h <- .resolve_h(h, domain_bbox, m1, m2)

  kappa2         <- 4 * (1 - rho) / (rho * h^2)
  c_operator     <- kappa2 + 4 / h^2            # = 4/(rho*h^2)
  sigma2_M_scale <- sigma2_SAR_scale * c_operator^2
  eff_range      <- sqrt(8 / kappa2)

  list(
    kappa2         = kappa2,
    sigma2_M_scale = sigma2_M_scale,
    eff_range      = eff_range,
    h              = h
  )
}

# ── 2. Matern -> SAR ──────────────────────────────────────────────────────────

#' Convert Matérn parameters to matched SAR parameters
#'
#' Inverse of \code{\link{sar_to_matern}}.
#'
#' @param kappa2         Squared range in \eqn{1/(\text{CRS units})^2}.
#' @param sigma2_M_scale Matérn SPDE amplitude scale.
#' @param h              Pixel side length. Supply or derive via
#'   \code{domain_bbox} + \code{m1}.
#' @param domain_bbox    Length-4 numeric bounding box.
#' @param m1             Number of grid rows.
#' @param m2             Number of grid cols.
#'
#' @return Named list: \code{rho}, \code{sigma2_SAR_scale}, \code{eff_range}, \code{h}.
#'
#' @seealso \code{\link{sar_to_matern}}, \code{\link{make_Q_cosine}}
#' @export
matern_to_sar <- function(kappa2, sigma2_M_scale,
                          h = NULL, domain_bbox = NULL, m1 = NULL, m2 = NULL) {
  if (!is.numeric(kappa2) || length(kappa2) != 1L || kappa2 <= 0)
    stop("`kappa2` must be a positive scalar.")
  if (!is.numeric(sigma2_M_scale) || length(sigma2_M_scale) != 1L ||
      sigma2_M_scale <= 0)
    stop("`sigma2_M_scale` must be a positive scalar.")

  h <- .resolve_h(h, domain_bbox, m1, m2)

  c_operator       <- kappa2 + 4 / h^2
  rho              <- (4 / h^2) / c_operator
  sigma2_SAR_scale <- sigma2_M_scale / c_operator^2
  eff_range        <- sqrt(8 / kappa2)

  if (rho <= 0 || rho >= 1)
    stop(sprintf("Implied rho=%.6f is outside (0,1). Check kappa2 and h.", rho))

  list(
    rho              = rho,
    sigma2_SAR_scale = sigma2_SAR_scale,
    eff_range        = eff_range,
    h                = h
  )
}

# ── 3. Cosine prior precision ─────────────────────────────────────────────────

#' Build the diagonal precision matrix for the Matérn cosine prior
#'
#' Constructs the \eqn{K \times K} diagonal precision matrix \eqn{Q} for the
#' physical-unit cosine basis matched to a Matérn field.
#'
#' The cosine basis \eqn{\phi_{j_1,j_2}} is normalised so that
#' \eqn{h^2 \Phi^\top \Phi = I_K} (physical L2 inner product). The diagonal
#' precision entry for mode \eqn{(j_1, j_2)} is:
#' \deqn{
#'   Q_{j_1,j_2} =
#'   \frac{(\kappa^2 + \lambda^\Delta_{j_1,j_2})^2}
#'        {\sigma^2_{M,\rm scale} \cdot L_x L_y}
#' }
#' where \eqn{L_x = m_1 h}, \eqn{L_y = m_2 h}, and
#' \eqn{\lambda^\Delta_{j_1,j_2} = \pi^2[(j_1-1)^2/L_x^2 + (j_2-1)^2/L_y^2]}.
#' The \eqn{L_x L_y} factor compensates for the physical normalisation of
#' \eqn{\Phi}, ensuring that the implied field variance matches the SAR prior.
#'
#' Mode ordering matches \code{fourier_integrate_basis()} from spatintegrate:
#' \eqn{j_1} and \eqn{j_2} run from 1 to \eqn{m_1} and \eqn{m_2} respectively,
#' with \eqn{j = 1} giving the constant term.
#'
#' @param kappa2         Squared range in \eqn{1/(\text{CRS units})^2}.
#' @param sigma2_M_scale Matérn SPDE amplitude scale. From \code{\link{sar_to_matern}}.
#' @param m1             Number of grid rows (max j1 index).
#' @param m2             Number of grid cols (max j2 index). Default \code{m1}.
#' @param h              Pixel side length in CRS units. Supply or derive via
#'   \code{domain_bbox} + \code{m1}.
#' @param domain_bbox    Length-4 numeric bounding box.
#' @param K              Modes to retain (by decreasing prior variance).
#'   Default \code{m1 * m2}.
#'
#' @return Named list:
#'   \describe{
#'     \item{\code{Q}}{Diagonal sparse \eqn{K \times K} \code{Matrix}.}
#'     \item{\code{modes}}{Data frame: \code{j1}, \code{j2}, \code{lam_D},
#'       \code{prior_var}, \code{q} — ordered by decreasing prior variance.
#'       Column order matches \code{fourier_integrate_basis()}.}
#'     \item{\code{h}}{Pixel size used.}
#'   }
#'
#' @seealso \code{\link{sar_to_matern}}, \code{\link{matern_to_sar}}
#' @export
make_Q_cosine <- function(kappa2, sigma2_M_scale, m1, m2 = m1,
                          h = NULL, domain_bbox = NULL, K = NULL) {
  if (!is.numeric(kappa2) || length(kappa2) != 1L || kappa2 <= 0)
    stop("`kappa2` must be a positive scalar.")
  if (!is.numeric(sigma2_M_scale) || length(sigma2_M_scale) != 1L ||
      sigma2_M_scale <= 0)
    stop("`sigma2_M_scale` must be a positive scalar.")
  if (!is.numeric(m1) || m1 < 1L)
    stop("`m1` must be a positive integer.")
  if (!is.numeric(m2) || m2 < 1L)
    stop("`m2` must be a positive integer.")

  h <- .resolve_h(h, domain_bbox, m1, m2)

  Lx <- m1 * h
  Ly <- m2 * h

  # Physical normalisation factor: Phi is normalised so h^2*Phi'Phi = I.
  # The prior variance of beta_k in field units needs compensating by Lx*Ly.
  sigma2_M_phys <- sigma2_M_scale * Lx * Ly

  if (is.null(K)) K <- m1 * m2
  K <- min(as.integer(K), m1 * m2)

  # Enumerate all m1*m2 modes (j1, j2 from 1, matching fourier_freq_grid)
  # Enumerate with j1 varying fastest to match fourier_freq_grid(J1,J2):
  # expand.grid(j1 = 1:J1, j2 = 1:J2) -- j1 is the inner loop.
  modes <- do.call(rbind, lapply(seq_len(m2), function(j2)
    do.call(rbind, lapply(seq_len(m1), function(j1) {
      lam_D     <- pi^2 * ((j1 - 1L)^2 / Lx^2 + (j2 - 1L)^2 / Ly^2)
      q         <- (kappa2 + lam_D)^2 / sigma2_M_phys
      prior_var <- 1 / q
      data.frame(j1 = as.integer(j1), j2 = as.integer(j2),
                 lam_D = lam_D, prior_var = prior_var, q = q)
    }))))

  modes <- modes[order(modes$q), ]
  rownames(modes) <- NULL
  modes <- modes[seq_len(K), ]

  list(
    Q     = Matrix::Diagonal(x = modes$q),
    modes = modes,
    h     = h
  )
}
