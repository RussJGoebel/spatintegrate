# integrate.R
# Layers 4 and 5: polygon integration for spatintegrate
#
# Two functions:
#
#   integrate_basis()        — user-facing, parallelises over polygons  [Layer 5]
#   .integrate_one_polygon() — internal per-polygon integration loop    [Layer 4]
#
# Internal helpers:
#
#   .extract_polygon_pieces() — pull flat POLYGON sfc from any geometry
#   .get_triangle_coords()    — extract 3x2 matrix from a triangle sfg
#   .infer_k()                — infer basis dimension from a pilot call


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Extract individual POLYGON pieces from any sfc, dropping empties and slivers.
.extract_polygon_pieces <- function(sfc, min_area = 0) {
  if (length(sfc) == 0L) return(sf::st_sfc())

  # Some geometry types can't be collection-extracted — guard with tryCatch
  pieces <- tryCatch(
    sf::st_collection_extract(sfc, "POLYGON") |> sf::st_cast("POLYGON"),
    error = function(e) sf::st_sfc()
  )

  if (length(pieces) == 0L) return(sf::st_sfc())

  areas  <- as.numeric(sf::st_area(pieces))
  pieces[areas > min_area]
}


# Extract a 3x2 coordinate matrix from a single triangle sfg.
# Strips the closing duplicate vertex sf appends to polygon rings.
.get_triangle_coords <- function(tri_sfg) {
  coords <- sf::st_coordinates(tri_sfg)

  # Exterior ring: L2 == 1 when present, otherwise all rows are exterior
  if ("L2" %in% colnames(coords)) {
    ext <- coords[coords[, "L2"] == 1L, c("X", "Y"), drop = FALSE]
  } else {
    ext <- coords[, c("X", "Y"), drop = FALSE]
  }

  # Remove closing duplicate
  if (nrow(ext) >= 2L && isTRUE(all(ext[1L, ] == ext[nrow(ext), ]))) {
    ext <- ext[-nrow(ext), , drop = FALSE]
  }

  unname(ext)
}


# Infer k by calling basis_fn on a single dummy point.
.infer_k <- function(basis_fn) {
  dummy  <- matrix(c(0, 0), nrow = 1L)
  result <- tryCatch(
    basis_fn(dummy),
    error = function(e) {
      stop(
        "Failed to infer basis dimension from basis_fn. ",
        "basis_fn must accept an n x 2 numeric matrix and return an n x k matrix. ",
        "Original error: ", conditionMessage(e)
      )
    }
  )
  if (!is.matrix(result) || nrow(result) != 1L) {
    stop(
      "basis_fn must return a matrix with one row per input point. ",
      "Got an object of class: ", paste(class(result), collapse = ", "), "."
    )
  }
  ncol(result)
}


# ------------------------------------------------------------------------------

#' Integrate a Basis Function Over a Single Polygon
#'
#' Computes one row of the integration matrix by evaluating \code{basis_fn}
#' over a single polygon, optionally intersecting with strata first. Returns
#' a length-k numeric vector — the area-weighted average of \code{basis_fn}
#' over the polygon.
#'
#' This is the per-polygon worker called by [integrate_basis()]. It can also
#' be called directly for debugging or custom parallelism.
#'
#' @param polygon_sfc A length-1 \code{sfc} of type POLYGON or MULTIPOLYGON,
#'   in a projected CRS.
#' @param basis_fn A function \code{function(coords)} where \code{coords} is
#'   an \eqn{n \times 2} numeric matrix of (x, y) coordinates in the same CRS
#'   as \code{polygon_sfc}. Must return an \eqn{n \times k} numeric matrix.
#' @param k Integer. Number of basis functions (columns returned by
#'   \code{basis_fn}). Typically inferred by [integrate_basis()].
#' @param strata_sf Optional \code{sf} object of strata polygons. If supplied,
#'   the polygon is intersected with each stratum before triangulating, which
#'   guarantees QMC points never straddle a stratum boundary. Must be in the
#'   same CRS as \code{polygon_sfc}.
#' @param qmc An \eqn{n \times 2} matrix of QMC points in \eqn{[0,1]^2},
#'   from [generate_qmc_unit_square()]. Generated once and shared across all
#'   polygons.
#' @param min_area Non-negative numeric. Triangles with area \code{<= min_area}
#'   are skipped. Default \code{0}.
#'
#' @return A length-k numeric vector. Returns a vector of \code{NA_real_} if
#'   the polygon produces no valid triangles (e.g. zero-area input).
#'
#' @seealso [integrate_basis()], [generate_qmc_unit_square()]
#'
#' @export
.integrate_one_polygon <- function(polygon_sfc, basis_fn, k, strata_sf = NULL,
                                   qmc, min_area = 0) {

  na_result <- rep(NA_real_, k)

  # --- Get pieces to triangulate ----------------------------------------------
  # Without strata: just the polygon itself.
  # With strata: one piece per stratum that intersects the polygon.

  if (is.null(strata_sf)) {
    pieces <- list(polygon_sfc)
  } else {
    strata_geom <- sf::st_geometry(strata_sf)
    pieces <- lapply(seq_along(strata_geom), function(i) {
      inter <- suppressWarnings(
        sf::st_intersection(polygon_sfc, strata_geom[i])
      )
      if (length(inter) == 0L || all(sf::st_is_empty(inter))) return(NULL)
      extracted <- .extract_polygon_pieces(inter, min_area = 0)
      if (length(extracted) == 0L) return(NULL)
      extracted
    })
    pieces <- Filter(Negate(is.null), pieces)
  }

  if (length(pieces) == 0L) return(na_result)

  # --- Triangulate all pieces -------------------------------------------------

  all_tris <- unlist(lapply(pieces, function(piece) {
    tris <- tryCatch(
      sf::st_triangulate_constrained(piece) |>
        sf::st_collection_extract("POLYGON") |>
        sf::st_cast("POLYGON"),
      error = function(e) sf::st_sfc()
    )
    if (length(tris) == 0L) return(NULL)
    # Return as plain list of sfg objects for easy iteration
    lapply(seq_along(tris), function(j) tris[[j]])
  }), recursive = FALSE)

  if (length(all_tris) == 0L) return(na_result)

  # --- Compute areas and drop degenerates ------------------------------------

  areas <- sapply(all_tris, function(tri) {
    as.numeric(sf::st_area(sf::st_sfc(tri, crs = sf::st_crs(polygon_sfc))))
  })

  keep       <- areas > min_area
  all_tris   <- all_tris[keep]
  areas      <- areas[keep]

  if (length(all_tris) == 0L) return(na_result)

  total_area <- sum(areas)
  if (total_area == 0) return(na_result)

  # --- Integration loop -------------------------------------------------------
  # One basis_fn call per triangle. phi_i is evaluated and immediately
  # discarded — memory stays bounded at n_per_triangle x k regardless of
  # how many triangles the polygon has.

  A_row <- numeric(k)

  for (j in seq_along(all_tris)) {
    coords <- .get_triangle_coords(all_tris[[j]])

    if (nrow(coords) != 3L) next   # malformed triangle — skip

    pts   <- map_unit_square_to_triangle(qmc, coords)
    phi   <- basis_fn(pts)                              # n_per_triangle x k

    if (!is.matrix(phi) || ncol(phi) != k) next        # basis_fn misbehaved

    weight <- areas[j] / total_area
    A_row  <- A_row + weight * colMeans(phi)
  }

  A_row
}


# ------------------------------------------------------------------------------

#' Integrate a Basis Function Over a Set of Polygons
#'
#' Computes an \eqn{n \times k} integration matrix \eqn{A} where
#' \eqn{A_{ij}} is the area-weighted average of basis function \eqn{j} over
#' polygon \eqn{i}:
#'
#' \deqn{A_{ij} = \frac{1}{|D_i|} \int_{D_i} \phi_j(s)\, ds}
#'
#' Integration is performed by triangulating each polygon (optionally
#' intersecting with strata first), then applying QMC quadrature with
#' \code{n_per_triangle} Sobol points per triangle.
#'
#' @section Strata:
#' When \code{strata_sf} is supplied, each polygon is intersected with each
#' stratum before triangulating. This guarantees that no QMC sample point
#' straddles a stratum boundary, which is essential when \code{basis_fn} has
#' discontinuities at stratum edges (e.g. a covariance kernel with zero
#' correlation across regions of different type).
#'
#' @section Sampling density:
#' Every triangle receives exactly \code{n_per_triangle} QMC points regardless
#' of its area. Sampling density is therefore controlled by triangulation
#' fineness — regions with more/finer triangles are sampled more densely.
#' To increase density in a specific region, supply finer strata there.
#'
#' @section Parallelism:
#' Set \code{parallel_plan = "ambient"} to parallelise over polygons using
#' the current \code{future::plan()}. Set a plan before calling:
#' \code{future::plan(future::multisession, workers = 4)}.
#'
#' @param basis_fn A function \code{function(coords)} where \code{coords} is
#'   an \eqn{n \times 2} numeric matrix of projected (x, y) coordinates.
#'   Must return an \eqn{n \times k} numeric matrix. The value of \eqn{k} is
#'   inferred automatically from a pilot call.
#' @param polygons_sf An \code{sf} object whose rows are the polygons to
#'   integrate over. Must be in a projected CRS.
#' @param strata_sf Optional \code{sf} object of strata polygons. Must be in
#'   the same projected CRS as \code{polygons_sf}.
#' @param n_per_triangle Positive integer. Number of QMC points per triangle.
#'   Default \code{16}. Higher values improve accuracy at the cost of more
#'   \code{basis_fn} evaluations.
#' @param min_area Non-negative numeric. Triangles with area \code{<= min_area}
#'   are skipped. Default \code{0}.
#' @param parallel_plan One of \code{"sequential"} or \code{"ambient"}.
#'   \code{"sequential"} runs in the current process. \code{"ambient"} uses
#'   the current \code{future::plan()}. Default \code{"sequential"}.
#'
#' @return A numeric matrix of dimensions \eqn{n \times k}, where \eqn{n} is
#'   \code{nrow(polygons_sf)} and \eqn{k} is inferred from \code{basis_fn}.
#'   Rows corresponding to polygons that produced no valid triangles contain
#'   \code{NA}.
#'
#' @seealso [ensure_projected()], [generate_qmc_unit_square()]
#'
#' @examples
#' # Build a small set of projected polygons to integrate over
#' polygons_sf <- sf::st_sf(
#'   id       = 1:3,
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(matrix(c(0,0, 2,0, 2,1, 0,1, 0,0), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(2,0, 4,0, 4,1, 2,1, 2,0), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(0,1, 2,1, 2,2, 0,2, 0,1), ncol=2, byrow=TRUE)))
#'   ),
#'   crs = 32619
#' )
#'
#' # A constant basis: returns 1 everywhere (k = 1)
#' # The integral of 1 over each polygon equals 1 (area-weighted average)
#' const_basis <- function(coords) matrix(1, nrow = nrow(coords), ncol = 1)
#' A <- integrate_basis(const_basis, polygons_sf)
#' # A is a 3x1 matrix of 1s
#'
#' # A two-dimensional basis: x-coordinate and y-coordinate averages (k = 2)
#' coord_basis <- function(coords) coords   # returns the coords themselves
#' A2 <- integrate_basis(coord_basis, polygons_sf)
#' # A2[i, 1] is the mean x-coordinate of polygon i
#' # A2[i, 2] is the mean y-coordinate of polygon i
#'
#' \dontrun{
#' # With strata — QMC points will never straddle stratum boundaries
#' strata_sf <- sf::st_sf(
#'   geometry = sf::st_sfc(
#'     sf::st_polygon(list(matrix(c(-1,-1, 2,-1, 2,3, -1,3, -1,-1), ncol=2, byrow=TRUE))),
#'     sf::st_polygon(list(matrix(c(2,-1,  5,-1, 5,3,  2,3,  2,-1), ncol=2, byrow=TRUE)))
#'   ),
#'   crs = 32619
#' )
#' A_strat <- integrate_basis(const_basis, polygons_sf, strata_sf = strata_sf)
#'
#' # Parallel execution
#' future::plan(future::multisession, workers = 4)
#' A_par <- integrate_basis(const_basis, polygons_sf, parallel_plan = "ambient")
#' future::plan(future::sequential)  # reset
#' }
#'
#' @export
integrate_basis <- function(basis_fn,
                            polygons_sf,
                            strata_sf        = NULL,
                            n_per_triangle   = 16L,
                            min_area         = 0,
                            parallel_plan    = c("sequential", "ambient")) {

  parallel_plan <- match.arg(parallel_plan)

  # --- Validate inputs --------------------------------------------------------

  if (!is.function(basis_fn)) {
    stop("`basis_fn` must be a function.")
  }
  if (!inherits(polygons_sf, "sf")) {
    stop("`polygons_sf` must be an sf object.")
  }
  if (nrow(polygons_sf) == 0L) {
    stop("`polygons_sf` must have at least one row.")
  }
  if (!is.numeric(n_per_triangle) || length(n_per_triangle) != 1L ||
      n_per_triangle < 1L || n_per_triangle != as.integer(n_per_triangle)) {
    stop("`n_per_triangle` must be a positive integer.")
  }
  if (!is.numeric(min_area) || length(min_area) != 1L ||
      is.na(min_area) || min_area < 0) {
    stop("`min_area` must be a non-negative numeric scalar.")
  }

  assert_projected(polygons_sf, arg_name = "polygons_sf")

  if (!is.null(strata_sf)) {
    if (!inherits(strata_sf, "sf")) {
      stop("`strata_sf` must be an sf object or NULL.")
    }
    assert_projected(strata_sf, arg_name = "strata_sf")
    ensure_same_crs(
      sf::st_crs(polygons_sf),
      sf::st_crs(strata_sf),
      context = "integrate_basis"
    )
  }

  # --- Infer k from basis_fn --------------------------------------------------

  k <- .infer_k(basis_fn)

  # --- Generate QMC points once -----------------------------------------------
  # Same sequence reused across all polygons and triangles.

  n_per_triangle <- as.integer(n_per_triangle)
  qmc            <- generate_qmc_unit_square(n_per_triangle)

  # --- Integrate over each polygon --------------------------------------------

  n       <- nrow(polygons_sf)
  geom    <- sf::st_geometry(polygons_sf)

  worker <- function(i) {
    .integrate_one_polygon(
      polygon_sfc = geom[i],
      basis_fn    = basis_fn,
      k           = k,
      strata_sf   = strata_sf,
      qmc         = qmc,
      min_area    = min_area
    )
  }

  if (parallel_plan == "sequential") {
    rows <- lapply(seq_len(n), worker)
  } else {
    rows <- future.apply::future_lapply(
      seq_len(n),
      worker,
      future.globals  = list(
        .integrate_one_polygon = .integrate_one_polygon,
        .extract_polygon_pieces = .extract_polygon_pieces,
        .get_triangle_coords   = .get_triangle_coords,
        map_unit_square_to_triangle = map_unit_square_to_triangle,
        geom     = geom,
        basis_fn = basis_fn,
        strata_sf = strata_sf,
        qmc      = qmc,
        k        = k,
        min_area = min_area
      ),
      future.packages = "sf",
      future.seed     = TRUE
    )
  }

  # --- Assemble result --------------------------------------------------------

  matrix(
    unlist(rows, use.names = FALSE),
    nrow  = n,
    ncol  = k,
    byrow = TRUE
  )
}
