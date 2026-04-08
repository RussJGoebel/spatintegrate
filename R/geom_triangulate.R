# geom_triangulate.R
# Layer 2: Polygon triangulation for spatintegrate
#
# Two functions:
#
#   triangulate_sf(polygon_sfc)       — triangulate a polygon into interior triangles
#   get_triangle_coords(triangle_sfg) — extract 3x2 coordinate matrix from a triangle
#
# Contract: all inputs must be in a projected CRS.
# Use ensure_projected() to convert geographic CRS before calling these functions.
#
# Requires GEOS >= 3.10 for sf::st_triangulate_constrained().
# Check with: sf::sf_extSoftVersion()[["GEOS"]]
#
# TODO (quality): for minimum-angle guarantees (eliminating thin sliver triangles
# from strata intersection pieces), consider sfdct::ct_triangulate() which wraps
# Shewchuk's Triangle library via RTriangle and supports Steiner point insertion.
# The swap would be isolated to triangulate_sf() — nothing else in the stack changes.
# See: https://cran.r-project.org/package=sfdct


# ------------------------------------------------------------------------------

#' Triangulate a Projected Polygon into Interior Triangles
#'
#' Decomposes a polygon (or multipolygon) into a set of triangles using
#' constrained Delaunay triangulation via \code{sf::st_triangulate_constrained()}.
#' All returned triangles are guaranteed to lie inside the input polygon —
#' constrained triangulation respects polygon boundary edges by construction,
#' so exterior triangles cannot occur.
#'
#' @section CRS requirement:
#' Input must be in a **projected CRS** (units of meters or feet). Geographic
#' CRS (lon/lat) produces incorrect results because the barycentric mapping in
#' [map_unit_square_to_triangle()] is Euclidean. Use [ensure_projected()] to
#' convert before calling this function.
#'
#' @section Requirements:
#' Requires GEOS >= 3.10. Check your version with
#' \code{sf::sf_extSoftVersion()[["GEOS"]]}. GEOS 3.10 was released in 2021
#' and is available in all current sf binary distributions.
#'
#' @section Quality note:
#' Constrained Delaunay triangulation maximises the minimum angle across all
#' triangles — the provably optimal shape quality from a fixed vertex set. For
#' additional quality control (guaranteed minimum angle via Steiner point
#' insertion to eliminate thin slivers), see the TODO comment in the source
#' for the \pkg{sfdct} upgrade path.
#'
#' @param x An \code{sf}, \code{sfc}, or \code{sfg} object of type POLYGON or
#'   MULTIPOLYGON. Must be in a projected CRS.
#' @param min_area Non-negative numeric scalar. Triangles with area
#'   \code{<= min_area} (in CRS units squared) are dropped as degenerate.
#'   Default \code{0}, dropping only exactly-zero-area triangles.
#'
#' @return An \code{sfc} of POLYGON geometries, each a triangle with exactly
#'   3 unique exterior vertices. All returned triangles are interior to the
#'   input polygon. The CRS is preserved from the input. Returns an empty
#'   \code{sfc} with a warning if no interior triangles are found.
#'
#' @seealso [get_triangle_coords()], [ensure_projected()],
#'   [map_unit_square_to_triangle()]
#'
#' @examples
#' \dontrun{
#' sq <- sf::st_sfc(sf::st_polygon(list(matrix(
#'   c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE
#' ))), crs = 32619)
#' tris <- triangulate_sf(sq)
#' length(tris)  # 2
#'
#' # L-shape — non-convex, interior triangles only
#' l_shape <- sf::st_sfc(sf::st_polygon(list(matrix(
#'   c(0,0, 2,0, 2,1, 1,1, 1,2, 0,2, 0,0), ncol = 2, byrow = TRUE
#' ))), crs = 32619)
#' tris <- triangulate_sf(l_shape)
#' }
#'
#' @export
triangulate_sf <- function(x, min_area = 0) {

  # --- Input validation -------------------------------------------------------

  if (!inherits(x, c("sf", "sfc", "sfg"))) {
    stop("`x` must be an sf, sfc, or sfg object.")
  }
  if (!is.numeric(min_area) || length(min_area) != 1L ||
      is.na(min_area) || min_area < 0) {
    stop("`min_area` must be a non-negative numeric scalar.")
  }

  # Coerce to sfc
  if (inherits(x, "sf"))  x <- sf::st_geometry(x)
  if (inherits(x, "sfg")) x <- sf::st_sfc(x)

  # CRS check — must be projected
  crs <- sf::st_crs(x)
  if (is.na(crs)) {
    stop(
      "Input has no CRS. Assign a projected CRS before calling triangulate_sf()."
    )
  }
  if (isTRUE(sf::st_is_longlat(x))) {
    stop(
      "Input is in a geographic CRS (lon/lat). ",
      "triangulate_sf() requires a projected CRS (units of meters). ",
      "Use ensure_projected() to convert first."
    )
  }

  # Geometry type check
  geom_types <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop(
      "`x` must contain only POLYGON or MULTIPOLYGON geometries. ",
      "Found: ", paste(geom_types, collapse = ", "), "."
    )
  }

  # --- Flatten multipolygons --------------------------------------------------

  if (any(geom_types == "MULTIPOLYGON")) {
    x <- sf::st_cast(x, "POLYGON")
  }

  # --- Triangulate ------------------------------------------------------------

  # Constrained Delaunay: boundary edges are respected by construction,
  # all output triangles are interior — no post-processing filter needed.
  # Requires GEOS >= 3.10.
  tris <- tryCatch(
    sf::st_triangulate_constrained(x) |>
      sf::st_collection_extract("POLYGON") |>
      sf::st_cast("POLYGON"),
    error = function(e) {
      stop(
        "sf::st_triangulate_constrained() failed: ", conditionMessage(e), "\n",
        "Requires GEOS >= 3.10. Check: sf::sf_extSoftVersion()[[\"GEOS\"]]"
      )
    }
  )

  if (length(tris) == 0L) {
    warning(
      "Triangulation produced no triangles. ",
      "Check that the input polygon has non-zero area."
    )
    return(sf::st_sfc(crs = crs))
  }

  # --- Drop degenerate triangles ----------------------------------------------

  if (min_area > 0) {
    areas <- as.numeric(sf::st_area(tris))
    tris  <- tris[areas > min_area]

    if (length(tris) == 0L) {
      warning(
        "All triangles have area <= min_area (", min_area, "). ",
        "Try reducing min_area."
      )
      return(sf::st_sfc(crs = crs))
    }
  }

  tris
}


# ------------------------------------------------------------------------------

#' Extract Vertex Coordinates from a Triangle Geometry
#'
#' Given a single triangle \code{sfg} or single-element \code{sfc}, returns
#' the three vertex coordinates as a \eqn{3 \times 2} numeric matrix suitable
#' for passing to [map_unit_square_to_triangle()].
#'
#' This function is called once per triangle inside the integration loop.
#' It strips the closing duplicate vertex that \code{sf} appends to polygon
#' coordinate sequences and validates that exactly three unique vertices remain.
#'
#' @param x A single triangle geometry: either an \code{sfg} of type POLYGON,
#'   or a length-1 \code{sfc}. Must have exactly 3 unique exterior vertices
#'   (i.e. be a triangle, not a general polygon).
#'
#' @return A \eqn{3 \times 2} numeric matrix of vertex coordinates
#'   \eqn{(x, y)}, one vertex per row, with no closing duplicate. Row names
#'   and column names are stripped (\code{unname()}).
#'
#' @seealso [triangulate_sf()], [map_unit_square_to_triangle()]
#'
#' @examples
#' \dontrun{
#' tris   <- triangulate_sf(some_polygon_sfc)
#' coords <- get_triangle_coords(tris[[1]])
#' dim(coords)  # 3 x 2
#' }
#'
#' @export
get_triangle_coords <- function(x) {

  # Coerce sfc to sfg
  if (inherits(x, "sfc")) {
    if (length(x) != 1L) {
      stop("`x` must be a single geometry: an sfg or a length-1 sfc.")
    }
    x <- x[[1L]]
  }

  if (!inherits(x, "sfg")) {
    stop("`x` must be an sfg or a length-1 sfc.")
  }

  coords <- sf::st_coordinates(x)

  # sf encodes polygon rings with an L1 (ring index) and L2 (polygon index)
  # column. The exterior ring is L2 == 1. For a simple triangle sfg produced
  # by triangulate_sf(), there is only one ring.
  if ("L2" %in% colnames(coords)) {
    ext <- coords[coords[, "L2"] == 1L, c("X", "Y"), drop = FALSE]
  } else if ("L1" %in% colnames(coords)) {
    ext <- coords[coords[, "L1"] == 1L, c("X", "Y"), drop = FALSE]
  } else {
    ext <- coords[, c("X", "Y"), drop = FALSE]
  }

  # Remove closing duplicate (sf repeats the first vertex to close the ring)
  if (nrow(ext) >= 2L && isTRUE(all(ext[1L, ] == ext[nrow(ext), ]))) {
    ext <- ext[-nrow(ext), , drop = FALSE]
  }

  if (nrow(ext) != 3L) {
    stop(sprintf(
      "Expected a triangle with 3 unique exterior vertices; got %d. ",
      nrow(ext)
    ))
  }

  unname(ext)
}
