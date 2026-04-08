#' Square grid over (buffered) geometry, preserving the input CRS
#'
#' @description
#' Create a regular grid of **square** polygons with cell size `res` over either a
#' buffered bounding rectangle or a buffered convex hull of `x`. The function does
#' **not** reproject: all distances are interpreted in the units of `st_crs(x)`.
#' If `x` is geographic (lon/lat), `res` is in **degrees**.
#'
#' @param x An `sf` or `sfc` geometry with a CRS.
#' @param res Positive numeric scalar. Square cell size in **CRS units**
#'   (degrees for lon/lat; meters/feet for projected).
#' @param buffer Nonnegative numeric scalar. Buffer distance applied to the masking
#'   geometry (rectangle or hull) in **CRS units**. Default `0`.
#' @param shape One of `"rectangle"` (buffered bbox; fastest) or `"hull"` (buffered convex hull).
#' @param align_origin Optional numeric length-2 `(x0, y0)` to align the grid lattice.
#'   If `NULL`, the grid anchors at the (buffered) extent's lower-left corner.
#' @param clip Logical. If `TRUE` and `shape = "hull"`, drop cells outside the buffered hull.
#'
#' @return An `sf` polygon grid with columns:
#' \itemize{
#'   \item `row`, `col` — 1-based lattice indices (stable under clipping).
#'   \item Attributes: `cellsize`, `origin`, `buffer`, `shape`, `crs_used`.
#' }
#'
#' @note On geographic CRS with `sf::sf_use_s2()` TRUE, `st_buffer()` is geodesic (meters)
#' while `res` is in degrees. Either disable s2 (`sf::sf_use_s2(FALSE)`), set `buffer = 0`,
#' or reproject beforehand if you need metric consistency.
#'
#' @examples
#' \dontrun{
#' g <- make_square_grid_in_crs(my_geom, res = 330, buffer = 0, shape = "rectangle")
#' g2 <- make_square_grid_in_crs(my_geom, res = 250, buffer = 500,
#'                               shape = "hull", align_origin = c(0, 0))
#' }
#' @export
make_square_grid_in_crs <- function(
    x,
    res,
    buffer = 0,
    shape = c("rectangle", "hull"),
    align_origin = NULL,
    clip = TRUE
) {
  stopifnot(inherits(x, c("sf", "sfc")))
  crs_in <- sf::st_crs(x)
  if (is.na(crs_in)) stop("`x` must have a CRS; assign one before calling.")
  if (!is.numeric(res) || length(res) != 1L || !is.finite(res) || res <= 0) {
    stop("`res` must be a positive finite scalar (in CRS units).")
  }
  if (!is.numeric(buffer) || length(buffer) != 1L || !is.finite(buffer) || buffer < 0) {
    stop("`buffer` must be a nonnegative finite scalar (in CRS units).")
  }

  shape <- match.arg(shape)
  crs_proj <- crs_in

  u <- sf::st_union(x)
  mask_geom <- if (shape == "hull") {
    sf::st_buffer(sf::st_convex_hull(u), dist = buffer)
  } else {
    rect <- sf::st_as_sfc(sf::st_bbox(u), crs = crs_proj)
    if (buffer > 0) sf::st_buffer(rect, dist = buffer) else rect
  }

  bb <- sf::st_bbox(mask_geom)
  if (is.null(align_origin)) {
    x0 <- as.numeric(bb["xmin"]); y0 <- as.numeric(bb["ymin"])
  } else {
    stopifnot(is.numeric(align_origin), length(align_origin) == 2L)
    x0 <- align_origin[1]; y0 <- align_origin[2]
  }

  eps    <- res * 1e-9
  width  <- as.numeric(bb["xmax"] - x0)
  height <- as.numeric(bb["ymax"] - y0)
  n_cols <- max(1L, ceiling((width  + eps) / res))
  n_rows <- max(1L, ceiling((height + eps) / res))

  grid <- sf::st_make_grid(
    mask_geom,
    cellsize = c(res, res),
    offset   = c(x0, y0),
    n        = c(n_cols, n_rows),
    square   = TRUE
  )
  if (length(grid) == 0L) stop("st_make_grid() produced no cells; check `res`, `buffer`, and the extent.")

  grid_sf <- sf::st_as_sf(grid)
  cents   <- sf::st_centroid(sf::st_geometry(grid_sf))
  xy      <- sf::st_coordinates(cents)

  grid_sf$col <- 1L + floor((xy[, "X"] - x0) / res + 1e-12)
  grid_sf$row <- 1L + floor((xy[, "Y"] - y0) / res + 1e-12)

  if (clip && shape == "hull") {
    keep_mat <- sf::st_intersects(grid_sf, mask_geom, sparse = FALSE)
    if (!is.matrix(keep_mat) || ncol(keep_mat) == 0L) {
      stop("No intersection between grid and mask; check `buffer`/`shape`.")
    }
    grid_sf <- grid_sf[keep_mat[, 1], , drop = FALSE]
  }

  attr(grid_sf, "cellsize") <- res
  attr(grid_sf, "origin")   <- c(x0 = x0, y0 = y0)
  attr(grid_sf, "buffer")   <- buffer
  attr(grid_sf, "shape")    <- shape
  attr(grid_sf, "crs_used") <- crs_proj$wkt

  grid_sf
}
