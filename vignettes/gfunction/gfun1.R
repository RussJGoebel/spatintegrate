# g_sensitivity.R
#
# Sensitivity analysis for the choice of g (sensor footprint response).
#
# Strategy:
#   - Define a bivariate Gaussian g centered on each footprint centroid
#   - For PLOTTING: evaluate g over the original footprint extent + padding
#   - For INTEGRATION: buffer each footprint by 3*sigma so g is negligible
#     at the integration boundary, then pass to integrate_basis()
#   - Compare resulting A matrices across several sigma scale factors
#
# OCO-2 footprint dimensions (approx): 1290m across-track x 2250m along-track
# Base sigmas set so the footprint edge is at ~1 sigma:
#   sigma_across = 1290 / 2 = 645m
#   sigma_along  = 2250 / 2 = 1125m
#
# Orientation convention:
#   make_demo_footprint constructs the footprint with the long axis along y,
#   then rotates CCW by angle_deg. make_gaussian_g takes angle_rad as the
#   CCW rotation of the footprint frame — we negate it so the Gaussian's
#   long axis aligns with the footprint's long axis in world coordinates.

library(sf)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Gaussian g factory
#
#    Returns a basis_fn compatible with integrate_basis():
#      input  : n x 2 matrix of projected (x, y) coordinates
#      output : n x 1 matrix of unnormalized Gaussian weights
#
#    sigma_across : std dev in the across-track direction (short axis, metres)
#    sigma_along  : std dev in the along-track direction  (long axis,  metres)
#    angle_rad    : CCW rotation of the footprint frame from x-axis (radians)
#                   Pass the NEGATED footprint angle if the footprint long
#                   axis starts along y (see demo below).
# ------------------------------------------------------------------------------
make_gaussian_g <- function(cx, cy, sigma_across, sigma_along, angle_rad = 0) {
  force(cx); force(cy)
  force(sigma_across); force(sigma_along); force(angle_rad)
  function(coords) {
    dx <- coords[, 1] - cx
    dy <- coords[, 2] - cy

    # Rotate world coords into footprint-aligned frame
    cos_a <- cos(angle_rad)
    sin_a <- sin(angle_rad)
    u <-  cos_a * dx + sin_a * dy   # across-track
    v <- -sin_a * dx + cos_a * dy   # along-track

    w <- exp(-0.5 * ((u / sigma_across)^2 + (v / sigma_along)^2))
    matrix(w, ncol = 1L)
  }
}

# ------------------------------------------------------------------------------
# 2. Buffer a footprint polygon for numerical integration
#    Buffer radius = n_sigma * max(sigma_across, sigma_along)
#    so g < exp(-n_sigma^2 / 2) at the boundary (~0.01% for n_sigma = 3)
# ------------------------------------------------------------------------------
buffer_footprint <- function(poly_sf, sigma_across, sigma_along, n_sigma = 3) {
  sf::st_buffer(poly_sf, dist = n_sigma * max(sigma_across, sigma_along))
}

# ------------------------------------------------------------------------------
# 3. Plot g over a footprint
#    poly_sf : ORIGINAL footprint (not buffered) — sets plot extent and outline
#    g_fn    : closure from make_gaussian_g()
#    pad     : padding beyond footprint bbox (metres) to show Gaussian tails
# ------------------------------------------------------------------------------
plot_g_on_footprint <- function(poly_sf, g_fn, n_grid = 150, pad = 600,
                                title = "Gaussian g over footprint") {
  bbox <- sf::st_bbox(poly_sf)
  xs   <- seq(bbox["xmin"] - pad, bbox["xmax"] + pad, length.out = n_grid)
  ys   <- seq(bbox["ymin"] - pad, bbox["ymax"] + pad, length.out = n_grid)
  grid <- expand.grid(x = xs, y = ys)

  grid$g <- as.numeric(g_fn(as.matrix(grid)))

  fp_coords <- as.data.frame(sf::st_coordinates(poly_sf))

  ggplot(grid, aes(x = x, y = y, fill = g)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma", name = "g(s)") +
    geom_polygon(
      data    = fp_coords,
      aes(x = X, y = Y, fill = NULL),
      colour  = "white", fill = NA, linewidth = 0.8
    ) +
    coord_equal() +
    labs(title    = title,
         subtitle = "White outline = footprint boundary",
         x        = "Easting (m)",
         y        = "Northing (m)") +
    theme_minimal()
}

# ------------------------------------------------------------------------------
# 4. Build a g-weighted A matrix
#
#    For each sounding i:
#      - make a Gaussian g centered on the sounding centroid
#      - buffer the sounding polygon for numerical integration
#      - call integrate_basis() to get g-weighted averages over fine-grid cells
#      - normalise so the row sums to 1 (consistent with uniform-weight A)
#
#    soundings_sf         : sf of footprint polygons (projected CRS)
#    fine_grid_sf         : sf of latent grid cells  (same CRS)
#    sigma_across, sigma_along : base sigmas in metres
#    angle_col            : column in soundings_sf with per-sounding rotation
#                           angle in radians (negated convention); NULL = 0
#    n_sigma              : buffer multiplier (default 3)
#    ...                  : passed to integrate_basis()
# ------------------------------------------------------------------------------
build_gaussian_A <- function(soundings_sf, fine_grid_sf,
                             sigma_across, sigma_along,
                             angle_col = NULL,
                             n_sigma   = 3,
                             ...) {
  # source("integrate.R"); source("qmc_utils.R")  # if not already loaded

  centroids <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(soundings_sf)))
  n         <- nrow(soundings_sf)
  rows      <- vector("list", n)

  for (i in seq_len(n)) {
    cx    <- centroids[i, 1]
    cy    <- centroids[i, 2]
    angle <- if (!is.null(angle_col)) soundings_sf[[angle_col]][i] else 0

    g_fn     <- make_gaussian_g(cx, cy, sigma_across, sigma_along, angle)
    buf_poly <- buffer_footprint(soundings_sf[i, ], sigma_across, sigma_along, n_sigma)

    # integrate_basis returns an m x 1 matrix of g-weighted averages
    # over each fine-grid cell from the perspective of sounding i
    raw       <- integrate_basis(basis_fn = g_fn, polygons_sf = fine_grid_sf, ...)
    r         <- as.numeric(raw)
    rows[[i]] <- r / sum(r, na.rm = TRUE)
  }

  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------
# 5. Demo
# ------------------------------------------------------------------------------

# Construct a rotated rectangular footprint.
# Long axis (height) starts along y, then rotate CCW by angle_deg.
make_demo_footprint <- function(cx = 0, cy = 0,
                                width     = 1290,
                                height    = 2250,
                                angle_deg = 15) {
  a  <- angle_deg * pi / 180
  dx <- width  / 2
  dy <- height / 2

  corners_local <- rbind(
    c(-dx, -dy), c(dx, -dy), c(dx, dy), c(-dx, dy), c(-dx, -dy)
  )

  # CCW rotation: long axis (y) rotates toward x
  cos_a <- cos(a); sin_a <- sin(a)
  rot   <- matrix(c(cos_a, sin_a, -sin_a, cos_a), 2, 2, byrow = TRUE)

  corners_world <- t(rot %*% t(corners_local))
  corners_world[, 1] <- corners_world[, 1] + cx
  corners_world[, 2] <- corners_world[, 2] + cy

  sf::st_sf(
    id       = 1L,
    geometry = sf::st_sfc(
      sf::st_polygon(list(corners_world)),
      crs = 32619   # UTM 19N — Boston
    )
  )
}

# --- Run demo -----------------------------------------------------------------

angle_deg <- 15
angle_rad <- angle_deg * pi / 180

# Base sigmas: footprint edge ~ 1 sigma
sigma_across <- 1290 / 4   # 645 m  (short axis)
sigma_along  <- 2250 / 4   # 1125 m (long axis)

fp       <- make_demo_footprint(angle_deg = angle_deg)
centroid <- sf::st_coordinates(sf::st_centroid(fp))

# Negate angle_rad: the footprint long axis starts along y before rotation,
# so the inverse rotation (world -> footprint frame) uses the negative angle
g_base <- make_gaussian_g(
  cx           = centroid[1],
  cy           = centroid[2],
  sigma_across = sigma_across,
  sigma_along  = sigma_along,
  angle_rad    = -angle_rad
)

# Plot 1: base sigma, g shown over original footprint extent + padding
p1 <- plot_g_on_footprint(
  fp, g_base,
  title = sprintf(
    "Gaussian g  (sigma_across = %.0fm, sigma_along = %.0fm)",
    sigma_across, sigma_along
  )
)
print(p1)

# Plot 2: three scale factors for supplement comparison
scale_factors <- c(0.5, 1.0, 2.0)

plots <- lapply(scale_factors, function(s) {
  sa <- sigma_across * s
  sl <- sigma_along  * s
  g  <- make_gaussian_g(centroid[1], centroid[2], sa, sl, -angle_rad)
  plot_g_on_footprint(
    fp, g,
    title = sprintf(
      "scale = %.1fx  (sigma_across = %.0f, sigma_along = %.0f)", s, sa, sl
    )
  )
})

for (p in plots) print(p)

# To combine into one figure with patchwork:
# library(patchwork)
# wrap_plots(plots, nrow = 1) +
#   plot_annotation(title = "Sensitivity to sigma scale factor")
