# compute_overlap_fractions.R
# Layer: Change-of-support aggregation matrix via exact area intersection
#
# One user-facing function:
#
#   compute_overlap_fractions(soundings, fine_grid, ...) — sparse n x m matrix
#
# Computes A[i, j] = |D_i ∩ D_j^(ℓ)| / |D_i|, the fraction of sounding i's
# area covered by fine-grid cell j. Assumes uniform sensor sensitivity (g_i = 1).
#
# For non-uniform sensitivity weighting, use integrate_basis() with a custom
# basis_fn instead.
#
# Requires: sf, Matrix, future.apply
# CRS: all inputs must be in a projected CRS. Use ensure_projected() first.


# ------------------------------------------------------------------------------

#' Compute Overlap Fractions Between Soundings and a Fine Grid
#'
#' Returns a sparse \eqn{n \times m} matrix \eqn{A} where \eqn{n} is the
#' number of coarse observations (soundings) and \eqn{m} is the number of
#' fine-grid cells, with entries
#'
#' \deqn{A_{ij} = \frac{|D_i \cap D_j^{(\ell)}|}{|D_i|}}
#'
#' This is the exact change-of-support discretization under uniform sensor
#' sensitivity (\eqn{g_i = 1}). Rows sum to 1 for soundings fully contained
#' within the fine grid, and to less than 1 for soundings that partially
#' extend beyond it.
#'
#' For non-uniform sensitivity (\eqn{g_i \neq 1}), use [integrate_basis()]
#' with a sensitivity-weighted basis function instead.
#'
#' @param soundings An \code{sf} object of coarse observation polygons
#'   (\eqn{n} rows). Must be in a projected CRS.
#' @param fine_grid An \code{sf} object of fine-resolution latent grid cells
#'   (\eqn{m} rows). Must be in the same projected CRS as \code{soundings}.
#' @param min_area Non-negative numeric. Intersection pieces with area
#'   \eqn{\leq} \code{min_area} (in CRS units squared) are treated as
#'   degenerate (e.g. line or point touches) and set to zero. Default
#'   \code{0}, dropping only exactly-zero-area intersections.
#' @param parallel Logical. If \code{TRUE}, parallelises over soundings using
#'   the current \code{future::plan()}. Set a plan before calling, e.g.
#'   \code{future::plan(future::multisession, workers = 4)}. Default
#'   \code{FALSE}.
#' @param sparse Logical. If \code{TRUE} (default), returns a
#'   \code{\link[Matrix]{sparseMatrix}}. If \code{FALSE}, returns a dense
#'   \code{base::matrix}.
#'
#' @return A numeric matrix of dimensions \eqn{n \times m}. Most entries are
#'   zero (each sounding overlaps only a small fraction of fine-grid cells).
#'   Rows corresponding to soundings with no overlap with the fine grid are
#'   all zero.
#'
#' @seealso [integrate_basis()], [ensure_projected()], [ensure_same_crs()]
#'
#' @examples
#' \dontrun{
#' fine_grid <- sf::st_sf(
#'   id = 1:4,
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(1,0, 2,0, 2,1, 1,1, 1,0), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(0,1, 1,1, 1,2, 0,2, 0,1), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(1,1, 2,1, 2,2, 1,2, 1,1), ncol=2, byrow=TRUE)))
#'   ), crs = 32619
#' )
#' soundings <- sf::st_sf(
#'   id = 1:2,
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(matrix(c(0,0, 2,0, 2,1, 0,1, 0,0), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(0,1, 2,1, 2,2, 0,2, 0,1), ncol=2, byrow=TRUE)))
#'   ), crs = 32619
#' )
#'
#' A <- compute_overlap_fractions(soundings, fine_grid)
#' # Each row sums to 1: soundings exactly tile the fine grid
#' rowSums(A)  # c(1, 1)
#'
#' # Parallel
#' future::plan(future::multisession, workers = 4)
#' A <- compute_overlap_fractions(soundings, fine_grid, parallel = TRUE)
#' future::plan(future::sequential)
#' }
#'
#' @export
compute_overlap_fractions <- function(soundings,
                                      fine_grid,
                                      min_area = 0,
                                      parallel = FALSE,
                                      sparse   = TRUE) {

  # --- Validate inputs --------------------------------------------------------

  if (!inherits(soundings, "sf")) stop("`soundings` must be an sf object.")
  if (!inherits(fine_grid, "sf")) stop("`fine_grid` must be an sf object.")
  if (nrow(soundings) == 0L)      stop("`soundings` must have at least one row.")
  if (nrow(fine_grid) == 0L)      stop("`fine_grid` must have at least one row.")
  if (!is.numeric(min_area) || length(min_area) != 1L ||
      is.na(min_area) || min_area < 0) {
    stop("`min_area` must be a non-negative numeric scalar.")
  }

  assert_projected(soundings, arg_name = "soundings")
  assert_projected(fine_grid, arg_name = "fine_grid")
  ensure_same_crs(
    sf::st_crs(soundings),
    sf::st_crs(fine_grid),
    context = "compute_overlap_fractions"
  )

  # --- Precompute shared quantities -------------------------------------------

  n              <- nrow(soundings)
  m              <- nrow(fine_grid)
  sounding_geom  <- sf::st_geometry(soundings)
  fine_geom      <- sf::st_geometry(fine_grid)
  sounding_areas <- as.numeric(sf::st_area(soundings))

  # Sparse candidate index: for each sounding, which fine cells does it touch?
  # This avoids n*m intersection calls — only overlapping pairs are computed.
  touches <- sf::st_intersects(soundings, fine_grid, sparse = TRUE)

  # --- Per-sounding worker ----------------------------------------------------

  worker <- function(i) {
    row   <- numeric(m)
    j_idx <- touches[[i]]

    if (length(j_idx) == 0L) return(row)

    inter <- suppressWarnings(
      sf::st_intersection(sounding_geom[i], fine_geom[j_idx])
    )

    if (length(inter) == 0L) return(row)

    areas <- as.numeric(sf::st_area(inter))

    # Drop degenerate intersections (point/line touches have zero area)
    keep <- areas > min_area
    if (!any(keep)) return(row)

    row[j_idx[keep]] <- areas[keep] / sounding_areas[i]
    row
  }

  # --- Parallelise or run sequentially ----------------------------------------

  if (parallel) {
    rows <- future.apply::future_lapply(
      seq_len(n),
      worker,
      future.globals = list(
        worker         = worker,
        sounding_geom  = sounding_geom,
        fine_geom      = fine_geom,
        touches        = touches,
        sounding_areas = sounding_areas,
        min_area       = min_area,
        m              = m
      ),
      future.packages = "sf",
      future.seed     = TRUE
    )
  } else {
    rows <- lapply(seq_len(n), worker)
  }

  # --- Assemble and return ----------------------------------------------------

  A <- matrix(unlist(rows, use.names = FALSE), nrow = n, ncol = m, byrow = TRUE)
  if (sparse) Matrix::Matrix(A, sparse = TRUE) else A
}
