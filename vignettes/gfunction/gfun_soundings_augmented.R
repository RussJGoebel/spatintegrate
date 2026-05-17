# build_g_A_matrix.R
#
# Builds the g-weighted aggregation matrix A using real OCO-2 sounding
# geometries from goebel2026::soundings_augmented.
#
# Each row A[i, j] = integral of g_i(s) over fine-grid cell j, normalised
# so rows sum to 1. This replaces the uniform area-intersection A from
# compute_overlap_fractions() with a smooth Gaussian-weighted version.
#
# The sensitivity function g_i is defined via an affine map from each
# sounding's quadrilateral to a canonical [-1,1]^2 space, where an
# isotropic Gaussian with std dev tau is evaluated. tau is a single
# global parameter:
#   tau = 1/3  =>  sounding edge sits at ~1 sigma  (recommended default)
#   tau -> inf =>  g approaches uniform (recovers area-intersection A)
#   tau -> 0   =>  g concentrates at sounding centroid
#
# Parallelised over soundings via future.apply (workers = 16).
#
# Requires: goebel2026, sf, Matrix, future, future.apply

library(goebel2026)
library(spatintegrate)
library(sf)
library(Matrix)
library(future)
library(future.apply)

# source("integrate.R")   # if not exported from goebel2026
# source("qmc_utils.R")

# ------------------------------------------------------------------------------
# make_affine_g()
#
# Returns a basis_fn for integrate_basis() that evaluates an affine-mapped
# Gaussian over a single sounding quadrilateral.
#
# sounding_sf : sf object, exactly one row, quadrilateral geometry
# tau         : std dev in canonical [-1,1]^2 space (global scale parameter)
# ------------------------------------------------------------------------------
make_affine_g <- function(sounding_sf, tau) {

  coords <- sf::st_coordinates(sounding_sf)
  xy     <- unique(coords[, c("X", "Y")])

  if (nrow(xy) != 4L)
    stop("Expected 4 unique vertices, got ", nrow(xy))

  cx <- mean(xy[, "X"])
  cy <- mean(xy[, "Y"])

  xy_c   <- sweep(xy, 2, c(cx, cy), "-")
  angles <- atan2(xy_c[, "Y"], xy_c[, "X"])
  ord    <- order(angles)
  xy_ord <- xy_c[ord, , drop = FALSE]

  canon <- matrix(c(
    1, -1,
    1,  1,
    -1,  1,
    -1, -1
  ), ncol = 2, byrow = TRUE)

  M_t   <- solve(crossprod(canon), crossprod(canon, xy_ord))
  M_inv <- solve(t(M_t))

  force(tau); force(cx); force(cy); force(M_inv)

  function(pts) {
    dx <- pts[, 1] - cx
    dy <- pts[, 2] - cy
    uv <- cbind(dx, dy) %*% t(M_inv)
    w  <- exp(-0.5 * rowSums(uv^2) / tau^2)
    matrix(w, ncol = 1L)
  }
}

# ------------------------------------------------------------------------------
# buffer_sounding()
#
# Buffers a sounding polygon so g is negligible at the integration boundary.
# Buffer radius = buf_factor * sqrt(area).
# ------------------------------------------------------------------------------
buffer_sounding <- function(sounding_sf, buf_factor = 3) {
  char_len <- sqrt(as.numeric(sf::st_area(sounding_sf)))
  sf::st_buffer(sounding_sf, dist = buf_factor * char_len)
}

# ------------------------------------------------------------------------------
# build_g_A()
#
# Builds the full n x m g-weighted aggregation matrix, parallelised over
# soundings using future.apply with whatever plan is set at call time.
#
# soundings_sf  : sf of sounding polygons (projected CRS), n rows
# fine_grid_sf  : sf of latent grid cells (same CRS), m rows
# tau           : global Gaussian scale parameter (default 1/3)
# buf_factor    : buffer multiplier for integration domain (default 3)
# sparse        : return a sparse Matrix (default TRUE)
# n_per_triangle: QMC points per triangle (default 16)
# verbose       : print worker count and tau on entry
# ------------------------------------------------------------------------------
build_g_A <- function(soundings_sf,
                      fine_grid_sf,
                      tau            = 1/3,
                      buf_factor     = 3,
                      sparse         = TRUE,
                      n_per_triangle = 16L,
                      verbose        = TRUE) {

  n <- nrow(soundings_sf)
  m <- nrow(fine_grid_sf)

  if (verbose) message(sprintf(
    "Building g-weighted A: n=%d soundings, m=%d cells, tau=%.3f, workers=%d",
    n, m, tau, future::nbrOfWorkers()
  ))

  # Pull geometry out once — passed as a global to avoid serialising the
  # full sf object (with attributes) on every worker call
  fine_geom <- sf::st_geometry(fine_grid_sf)

  # Assign package-internal functions as local variables so future's
  # auto-detection can find them, and so workers receive them explicitly
  integrate_basis             <- spatintegrate:::integrate_basis
  .integrate_one_polygon      <- spatintegrate:::.integrate_one_polygon
  .extract_polygon_pieces     <- spatintegrate:::.extract_polygon_pieces
  .get_triangle_coords        <- spatintegrate:::.get_triangle_coords
  map_unit_square_to_triangle <- spatintegrate:::map_unit_square_to_triangle
  generate_qmc_unit_square    <- spatintegrate:::generate_qmc_unit_square

  rows <- future.apply::future_lapply(
    seq_len(n),
    function(i) {
      s_i      <- soundings_sf[i, ]
      g_fn     <- make_affine_g(s_i, tau = tau)
      buf_poly <- buffer_sounding(s_i, buf_factor = buf_factor)

      hits <- unlist(sf::st_intersects(buf_poly, fine_geom, sparse = TRUE))

      if (length(hits) == 0L)
        return(list(j = integer(0), x = numeric(0)))

      raw <- integrate_basis(
        basis_fn       = g_fn,
        polygons_sf    = fine_grid_sf[hits, ],
        n_per_triangle = n_per_triangle,
        parallel_plan  = "sequential"  # parallelism is over soundings
      )

      r <- as.numeric(raw)
      s <- sum(r, na.rm = TRUE)
      if (s > 0) r <- r / s

      list(j = hits, x = r)
    },
    future.seed     = TRUE,
    future.packages = c("sf", "goebel2026", "spatintegrate"),
    future.globals  = list(
      soundings_sf                 = soundings_sf,
      fine_grid_sf                 = fine_grid_sf,
      fine_geom                    = fine_geom,
      tau                          = tau,
      buf_factor                   = buf_factor,
      n_per_triangle               = n_per_triangle,
      make_affine_g                = make_affine_g,
      buffer_sounding              = buffer_sounding,
      integrate_basis              = integrate_basis,
      .integrate_one_polygon       = .integrate_one_polygon,
      .extract_polygon_pieces      = .extract_polygon_pieces,
      .get_triangle_coords         = .get_triangle_coords,
      map_unit_square_to_triangle  = map_unit_square_to_triangle,
      generate_qmc_unit_square     = generate_qmc_unit_square
    )
  )

  # Assemble sparse matrix from per-row (j, x) lists
  i_idx <- rep(seq_len(n), vapply(rows, function(r) length(r$j), integer(1)))
  j_idx <- unlist(lapply(rows, `[[`, "j"))
  x_val <- unlist(lapply(rows, `[[`, "x"))

  A <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                            dims = c(n, m))
  if (!sparse) A <- as.matrix(A)
  A
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

# Set up parallel workers
future::plan(future::multisession, workers = 16)

# Load and project data
soundings_proj <- sf::st_transform(soundings_augmented, crs = 32619)
fine_grid_proj <- sf::st_transform(target_grid,         crs = 32619)

# Uniform baseline
A_uniform <- compute_overlap_fractions(soundings_proj, fine_grid_proj)

# g-weighted A matrices across tau values
tau_values <- c(1/5, 1/3, 1/2)

A_list <- lapply(tau_values, function(tau) {
  message(sprintf("\n--- tau = %.3f ---", tau))
  build_g_A(
    soundings_sf   = soundings_proj,
    fine_grid_sf   = fine_grid_proj,
    tau            = tau,
    buf_factor     = 3,
    n_per_triangle = 16L,
    verbose        = TRUE
  )
})
names(A_list) <- paste0("tau_", round(tau_values, 3))

# Reset plan
future::plan(future::sequential)

# Sanity check: row sums should be ~1
for (nm in names(A_list)) {
  rs <- Matrix::rowSums(A_list[[nm]])
  message(sprintf("%s: row sum range [%.4f, %.4f]", nm, min(rs), max(rs)))
}

# Save
saveRDS(
  list(A_uniform = A_uniform, A_g = A_list, tau_values = tau_values),
  file = "A_matrices_g_sensitivity.rds"
)

message("Done. Saved to A_matrices_g_sensitivity.rds")
