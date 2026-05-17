# build_g_A_matrix.R
#
# Builds a g-weighted aggregation matrix using the updated spatintegrate
# integrate_basis() with weight_fn support (PR2) and batch triangulation (PR1).
#
# For each sounding i:
#   1. Find hit cells via st_intersects
#   2. Compute intersection polygons D_i ∩ D_j
#   3. Call integrate_basis(const_basis, inter_sf, weight_fn = g_fn)
#      which returns the g-weighted average of 1 over each intersection —
#      i.e. effectively just the g-weighted area fraction
#   4. Multiply by intersection area, row-normalise
#
# tau is the single global scale parameter for the affine Gaussian g.

library(goebel2026)
library(spatintegrate)
library(sf)
library(Matrix)
library(future)
library(future.apply)

# ------------------------------------------------------------------------------
# make_affine_g_fn()
#
# Returns a weight_fn compatible with the new integrate_basis() weight_fn arg:
#   input  : n x 2 coordinate matrix
#   output : length-n numeric vector of Gaussian weights
# ------------------------------------------------------------------------------
make_affine_g_fn <- function(sounding_sf, tau) {
  coords <- sf::st_coordinates(sounding_sf)
  xy     <- unique(coords[, c("X", "Y")])
  if (nrow(xy) != 4L) stop("Expected 4 unique vertices, got ", nrow(xy))

  cx <- mean(xy[, "X"]); cy <- mean(xy[, "Y"])
  xy_c   <- sweep(xy, 2, c(cx, cy), "-")
  xy_ord <- xy_c[order(atan2(xy_c[, "Y"], xy_c[, "X"])), , drop = FALSE]
  canon  <- matrix(c(1,-1, 1,1, -1,1, -1,-1), ncol = 2, byrow = TRUE)
  M_inv  <- solve(t(solve(crossprod(canon), crossprod(canon, xy_ord))))

  force(tau); force(cx); force(cy); force(M_inv)

  # Returns a plain numeric vector (weight_fn contract)
  function(pts) {
    dx <- pts[, 1] - cx; dy <- pts[, 2] - cy
    uv <- cbind(dx, dy) %*% t(M_inv)
    as.numeric(exp(-0.5 * rowSums(uv^2) / tau^2))
  }
}

# ------------------------------------------------------------------------------
# build_g_A()
# ------------------------------------------------------------------------------
build_g_A <- function(soundings_sf, fine_grid_sf,
                      tau            = 1/3,
                      n_per_triangle = 16L) {

  n <- nrow(soundings_sf)
  m <- nrow(fine_grid_sf)
  message(sprintf("Building g-weighted A: n=%d, m=%d, tau=%.3f, workers=%d",
                  n, m, tau, future::nbrOfWorkers()))

  sounding_geom <- sf::st_geometry(soundings_sf)
  fine_geom     <- sf::st_geometry(fine_grid_sf)
  crs           <- sf::st_crs(soundings_sf)
  touches       <- sf::st_intersects(soundings_sf, fine_grid_sf, sparse = TRUE)

  # Constant basis: integrate_basis with weight_fn=g gives g-weighted average
  # of 1 over each polygon = effectively the normalised g mass per cell
  const_basis <- function(coords) matrix(1, nrow = nrow(coords), ncol = 1L)

  worker <- function(i) {
    js <- as.integer(touches[[i]])
    if (length(js) == 0L)
      return(list(j = integer(0), x = numeric(0)))

    inter_geoms <- suppressWarnings(
      sf::st_intersection(sounding_geom[i], fine_geom[js])
    )
    areas <- as.numeric(sf::st_area(inter_geoms))
    keep  <- !sf::st_is_empty(inter_geoms) & areas > 0
    if (!any(keep))
      return(list(j = integer(0), x = numeric(0)))

    inter_sf   <- sf::st_sf(geometry = inter_geoms[keep], crs = crs)
    areas_keep <- areas[keep]
    js_keep    <- js[keep]

    g_fn <- make_affine_g_fn(soundings_sf[i, ], tau = tau)

    # integrate_basis with weight_fn: returns g-weighted average of 1
    # over each intersection polygon — uses batch triangulation (PR1)
    # and weight_fn weighting (PR2)
    g_ij <- as.numeric(integrate_basis(
      basis_fn       = const_basis,
      polygons_sf    = inter_sf,
      weight_fn      = g_fn,
      n_per_triangle = n_per_triangle
    ))

    weights <- g_ij * areas_keep
    s       <- sum(weights, na.rm = TRUE)
    if (s > 0) weights <- weights / s

    list(j = js_keep, x = weights)
  }

  rows <- future.apply::future_lapply(
    seq_len(n),
    worker,
    future.seed     = TRUE,
    future.packages = c("sf", "spatintegrate"),
    future.globals  = list(
      soundings_sf   = soundings_sf,
      fine_geom      = fine_geom,
      crs            = crs,
      touches        = touches,
      tau            = tau,
      n_per_triangle = n_per_triangle,
      const_basis    = const_basis,
      make_affine_g_fn = make_affine_g_fn,
      integrate_basis  = integrate_basis
    )
  )

  i_idx <- rep(seq_len(n), vapply(rows, function(r) length(r$j), integer(1)))
  j_idx <- unlist(lapply(rows, `[[`, "j"))
  x_val <- unlist(lapply(rows, `[[`, "x"))

  Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(n, m))
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
fine_grid_proj <- spatintegrate::ensure_projected(goebel2026::target_grid)

future::plan(future::multisession, workers = 16)

# Test on 10 soundings
message("--- test: 10 soundings ---")
system.time({
  A_test <- build_g_A(soundings_proj[1:10, ], fine_grid_proj,
                      tau = 1/3, n_per_triangle = 16L)
})

rs <- Matrix::rowSums(A_test)
message(sprintf("Row sum range: [%.4f, %.4f]", min(rs), max(rs)))

# Full run
tau_values <- c(1/5, 1/3, 1/2)

A_list <- lapply(tau_values, function(tau) {
  message(sprintf("\n--- tau = %.3f ---", tau))
  build_g_A(soundings_proj, fine_grid_proj, tau = tau, n_per_triangle = 16L)
})
names(A_list) <- paste0("tau_", round(tau_values, 3))

future::plan(future::sequential)

A_uniform <- spatintegrate::compute_overlap_fractions(soundings_proj, fine_grid_proj)

for (nm in names(A_list)) {
  rs   <- Matrix::rowSums(A_list[[nm]])
  diff <- norm(A_list[[nm]] - A_uniform, type = "F")
  message(sprintf("%s: row sums [%.4f, %.4f]  Frobenius diff = %.4f",
                  nm, min(rs), max(rs), diff))
}

saveRDS(
  list(A_uniform = A_uniform, A_g = A_list, tau_values = tau_values),
  file = "A_matrices_g_sensitivity.rds"
)
message("Done.")
