# sar_matern.R
#
# SAR <-> Matérn parameter matching and cosine prior builder.
#
# Two h conventions are supported:
#
#   Physical h (CRS units, e.g. metres):
#     kappa2_phys = 4*(1-rho) / (rho * h_phys^2)
#     Used by sar_to_matern() / matern_to_sar() / make_Q_cosine()
#
#   Normalised h (unit-square frame, h_norm = 1/n_cells_side):
#     kappa2_norm = 4*(1-rho) / (rho * h_norm^2)  [same formula, different h]
#     kappa2_phys = kappa2_norm / L^2              [L = max(Lx, Ly)]
#     Used by get_matern_parameters_from_SAR() / fourier_matern_Q_sb()
#
# The two are related by:
#   kappa2_norm = kappa2_phys * L^2
#   TAU         = 1 / (sigma2_M_scale * Lx * Ly)
#
# Verified (verify_integration_paths3.R):
#   fourier_matern_Q_sb() matches SpatialBasis::compute_precision() exactly
#   (ratio sd < 1e-6 across all K modes).
#
# Recommended workflow (unit-square, for SpatialBasis compatibility)
# ------------------------------------------------------------------
#   h_norm <- 1 / n_cells_side
#   p      <- get_matern_parameters_from_SAR(rho, sar_precision, h = h_norm)
#   fg     <- fourier_freq_grid(J1, J2, domain_bbox, norm = "unit_square")
#   Q      <- fourier_matern_Q_sb(fg, kappa2_norm = p$kappa2, TAU = p$matern_precision)
#   A      <- fourier_integrate_basis(polygons_sf, fg)
#   fit    <- fit_fastblm(y, A, Q$Q, phi = 1,
#                         solver = "woodbury", Q_inv = function(v) v / Q$q_diag)
#
# Recommended workflow (physical, original fastblm workflow)
# ----------------------------------------------------------
#   p      <- sar_to_matern(rho, sigma2_SAR_scale, h = h_phys)
#   fg     <- fourier_freq_grid(J1, J2, domain_bbox, norm = "physical")
#   Q      <- make_Q_cosine(p$kappa2, p$sigma2_M_scale, m1, m2, h = h_phys)
#   phi    <- sigma2_M / sigma2e   # from CV or SAR fit
#   fit    <- fit_fastblm(y, A, Q$Q, phi,
#                         solver = "woodbury", Q_inv = function(v) v / Q$q_diag)


# ------------------------------------------------------------------------------
# Internal helper
# ------------------------------------------------------------------------------

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
          "`h` (%.6g) inconsistent with `domain_bbox`/`m1` (implied h=%.6g).",
          h, h_check))
    }
    return(h)
  }
  if (!is.null(domain_bbox) && !is.null(m1)) {
    if (!is.numeric(domain_bbox) || length(domain_bbox) != 4L)
      stop("`domain_bbox` must be length-4 (xmin,ymin,xmax,ymax).")
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
          "Pixels not square: Lx/m1=%.6g, Ly/m2=%.6g. SAR matching assumes square pixels.",
          hx, hy))
    }
    return(hx)
  }
  stop("Supply either `h` or both `domain_bbox` and `m1`.")
}


# ------------------------------------------------------------------------------
# 1. Unit-square convention  (SpatialBasis-compatible)
# ------------------------------------------------------------------------------

#' Convert SAR parameters to Matern parameters (unit-square frame)
#'
#' All quantities in the normalised [0,1]^2 domain.
#' Use h = 1/n_cells_side (normalised pixel size).
#'
#' Matching formulas (rook adjacency, d=2):
#'   kappa2_norm     = 4*(1-rho) / (rho * h_norm^2)
#'   matern_variance = sigma2_SAR * (2d/(rho*h_norm))^2
#'   TAU             = 1 / matern_variance
#'
#' kappa2_norm relates to physical kappa2 by: kappa2_phys = kappa2_norm / L^2
#' where L = max(Lx, Ly) is the square embedding side length.
#'
#' @param rho           SAR autocorrelation in (0,1).
#' @param sar_precision SAR prior precision tau (= 1/sigma2_SAR).
#' @param h             Normalised pixel size = 1/n_cells_side.
#' @param d             Dimension (default 2).
#'
#' @return Named list: matern_precision (TAU), kappa2 (unit-square units),
#'   sar_variance, matern_variance.
#'
#' @seealso fourier_matern_Q_sb, sar_to_matern
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
    matern_precision = 1 / matern_variance,
    kappa2           = kappa2,
    sar_variance     = sar_variance,
    matern_variance  = matern_variance
  )
}


#' Convert Matern parameters back to SAR parameters (unit-square frame)
#'
#' Inverse of get_matern_parameters_from_SAR().
#'
#' @param matern_precision TAU from get_matern_parameters_from_SAR().
#' @param kappa2_norm      kappa^2 in [0,1]^2 units.
#' @param h                Normalised pixel size = 1/n_cells_side.
#' @param d                Dimension (default 2).
#'
#' @return Named list: sar_precision, rho.
#' @export
get_SAR_parameters_from_matern <- function(matern_precision, kappa2_norm, h, d = 2) {
  if (!is.numeric(matern_precision) || matern_precision <= 0)
    stop("`matern_precision` must be positive.")
  if (!is.numeric(kappa2_norm) || kappa2_norm <= 0)
    stop("`kappa2_norm` must be positive.")
  if (!is.numeric(h) || h <= 0)
    stop("`h` must be positive.")

  rho          <- 2 * d / (kappa2_norm * h^2 + 2 * d)
  sar_variance <- (1 / matern_precision) / (2 * d / (rho * h))^2

  list(
    sar_precision = 1 / sar_variance,
    rho           = rho
  )
}


# ------------------------------------------------------------------------------
# 2. Physical convention  (original fastblm workflow)
# ------------------------------------------------------------------------------

#' Convert SAR parameters to matched Matérn parameters (physical units)
#'
#' Uses spectral matching in physical CRS units.
#' sigma2_SAR_scale = phi * sigma2e / tau (field variance shape parameter).
#'
#' @param rho              SAR autocorrelation in (0,1).
#' @param sigma2_SAR_scale SAR prior covariance scale (= 1/tau when phi*sigma2e=1).
#' @param h                Pixel side length in CRS units.
#' @param domain_bbox      Length-4 numeric (xmin,ymin,xmax,ymax). Alternative to h.
#' @param m1               Number of grid rows. Used with domain_bbox.
#' @param m2               Number of grid cols. Used for square-pixel check.
#'
#' @return Named list: kappa2, sigma2_M_scale, eff_range, h.
#'
#' @seealso matern_to_sar, make_Q_cosine, get_matern_parameters_from_SAR
#' @export
sar_to_matern <- function(rho, sigma2_SAR_scale,
                          h = NULL, domain_bbox = NULL, m1 = NULL, m2 = NULL) {
  if (rho <= 0 || rho >= 1)
    stop("`rho` must be in (0, 1).")
  if (!is.numeric(sigma2_SAR_scale) || length(sigma2_SAR_scale) != 1L ||
      sigma2_SAR_scale <= 0)
    stop("`sigma2_SAR_scale` must be a positive scalar.")

  h          <- .resolve_h(h, domain_bbox, m1, m2)
  kappa2     <- 4 * (1 - rho) / (rho * h^2)
  c_op       <- kappa2 + 4 / h^2      # = 4/(rho*h^2)
  sigma2_M   <- sigma2_SAR_scale * c_op^2
  eff_range  <- sqrt(8 / kappa2)

  list(kappa2 = kappa2, sigma2_M_scale = sigma2_M, eff_range = eff_range, h = h)
}


#' Convert Matérn parameters to matched SAR parameters (physical units)
#'
#' Inverse of sar_to_matern().
#'
#' @param kappa2         kappa^2 in physical units (1/m^2).
#' @param sigma2_M_scale Matérn SPDE amplitude scale.
#' @param h              Pixel side length in CRS units.
#' @param domain_bbox    Length-4 numeric bounding box.
#' @param m1             Number of grid rows.
#' @param m2             Number of grid cols.
#'
#' @return Named list: rho, sigma2_SAR_scale, eff_range, h.
#' @export
matern_to_sar <- function(kappa2, sigma2_M_scale,
                          h = NULL, domain_bbox = NULL, m1 = NULL, m2 = NULL) {
  if (!is.numeric(kappa2) || length(kappa2) != 1L || kappa2 <= 0)
    stop("`kappa2` must be a positive scalar.")
  if (!is.numeric(sigma2_M_scale) || length(sigma2_M_scale) != 1L ||
      sigma2_M_scale <= 0)
    stop("`sigma2_M_scale` must be a positive scalar.")

  h        <- .resolve_h(h, domain_bbox, m1, m2)
  c_op     <- kappa2 + 4 / h^2
  rho      <- (4 / h^2) / c_op
  sigma2_S <- sigma2_M_scale / c_op^2
  eff_range <- sqrt(8 / kappa2)

  if (rho <= 0 || rho >= 1)
    stop(sprintf("Implied rho=%.6f outside (0,1). Check kappa2 and h.", rho))

  list(rho = rho, sigma2_SAR_scale = sigma2_S, eff_range = eff_range, h = h)
}


# ------------------------------------------------------------------------------
# 3. Bridge: convert between conventions
# ------------------------------------------------------------------------------

#' Convert physical Matern parameters to unit-square TAU and kappa2_norm
#'
#' Given physical kappa2 and sigma2_M_scale from sar_to_matern(), and the
#' domain bounding box, returns the unit-square parameters for fourier_matern_Q_sb().
#'
#' Relationships:
#'   kappa2_norm = kappa2_phys * L^2        (L = max(Lx,Ly))
#'   TAU         = 1 / (sigma2_M_scale * Lx * Ly)
#'
#' @param kappa2_phys    Physical kappa^2 from sar_to_matern().
#' @param sigma2_M_scale Matern scale from sar_to_matern().
#' @param domain_bbox    Length-4 numeric (xmin,ymin,xmax,ymax).
#'
#' @return Named list: kappa2_norm, TAU, L, Lx, Ly.
#' @export
physical_to_unit_square <- function(kappa2_phys, sigma2_M_scale, domain_bbox) {
  Lx <- domain_bbox[3L] - domain_bbox[1L]
  Ly <- domain_bbox[4L] - domain_bbox[2L]
  L  <- max(Lx, Ly)
  list(
    kappa2_norm = kappa2_phys * L^2,
    TAU         = 1 / (sigma2_M_scale * Lx * Ly),
    L = L, Lx = Lx, Ly = Ly
  )
}


#' Convert unit-square parameters to physical Matern parameters
#'
#' Inverse of physical_to_unit_square().
#'
#' @param kappa2_norm Unit-square kappa^2.
#' @param TAU         Unit-square precision.
#' @param domain_bbox Length-4 numeric (xmin,ymin,xmax,ymax).
#'
#' @return Named list: kappa2_phys, sigma2_M_scale.
#' @export
unit_square_to_physical <- function(kappa2_norm, TAU, domain_bbox) {
  Lx <- domain_bbox[3L] - domain_bbox[1L]
  Ly <- domain_bbox[4L] - domain_bbox[2L]
  L  <- max(Lx, Ly)
  list(
    kappa2_phys    = kappa2_norm / L^2,
    sigma2_M_scale = 1 / (TAU * Lx * Ly)
  )
}


# ------------------------------------------------------------------------------
# 4. Prior precision (physical convention)
# ------------------------------------------------------------------------------

#' Build diagonal prior precision matrix for Matern cosine basis (physical units)
#'
#' Diagonal entry for mode (j1, j2):
#'   q_k = (kappa2 + lambda_k)^2 / (sigma2_M_scale * Lx * Ly)
#' where lambda_k = pi^2 * ((j1-1)^2/Lx^2 + (j2-1)^2/Ly^2).
#'
#' Use with fourier_freq_grid(..., norm="physical") and phi = sigma2_M/sigma2e.
#' Mode ordering matches fourier_integrate_basis(): j1 varies fastest.
#'
#' When K < m1*m2 modes are requested, the K modes with largest prior variance
#' (smallest q) are retained and sorted by decreasing prior variance.
#'
#' @param kappa2         kappa^2 in physical units. From sar_to_matern().
#' @param sigma2_M_scale Matern scale. From sar_to_matern().
#' @param m1             Number of grid rows (J1 in fourier_freq_grid).
#' @param m2             Number of grid cols (J2). Default m1.
#' @param h              Pixel side length. Supply or derive via domain_bbox+m1.
#' @param domain_bbox    Length-4 numeric bounding box.
#' @param K              Modes to retain. Default m1*m2 (all).
#'
#' @return Named list: Q (diagonal sparse Matrix), q_diag, modes (data frame).
#'
#' @seealso sar_to_matern, fourier_matern_Q_sb
#' @export
make_Q_cosine <- function(kappa2, sigma2_M_scale, m1, m2 = m1,
                          h = NULL, domain_bbox = NULL, K = NULL) {
  if (!is.numeric(kappa2) || length(kappa2) != 1L || kappa2 <= 0)
    stop("`kappa2` must be a positive scalar.")
  if (!is.numeric(sigma2_M_scale) || length(sigma2_M_scale) != 1L ||
      sigma2_M_scale <= 0)
    stop("`sigma2_M_scale` must be a positive scalar.")
  if (!is.numeric(m1) || m1 < 1L) stop("`m1` must be a positive integer.")
  if (!is.numeric(m2) || m2 < 1L) stop("`m2` must be a positive integer.")

  h  <- .resolve_h(h, domain_bbox, m1, m2)
  Lx <- m1 * h
  Ly <- m2 * h
  sigma2_M_phys <- sigma2_M_scale * Lx * Ly

  full_basis <- is.null(K) || as.integer(K) >= m1 * m2
  if (is.null(K)) K <- m1 * m2
  K    <- min(as.integer(K), m1 * m2)
  sort <- !full_basis   # sort only when truncating

  # j1 varies fastest to match fourier_freq_grid(J1=m1, J2=m2) ordering
  grid  <- expand.grid(j1 = seq_len(m1), j2 = seq_len(m2))
  lam_D <- pi^2 * ((grid$j1 - 1L)^2 / Lx^2 + (grid$j2 - 1L)^2 / Ly^2)
  q_vec <- (kappa2 + lam_D)^2 / sigma2_M_phys

  modes <- data.frame(
    j1        = as.integer(grid$j1),
    j2        = as.integer(grid$j2),
    lam_D     = lam_D,
    prior_var = 1 / q_vec,
    q         = q_vec
  )

  if (sort) modes <- modes[order(modes$q), ]
  modes <- modes[seq_len(K), ]
  rownames(modes) <- NULL

  list(
    Q      = Matrix::Diagonal(x = modes$q),
    q_diag = modes$q,
    modes  = modes,
    h      = h
  )
}
