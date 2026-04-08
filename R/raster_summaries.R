#' Summarize a SpatRaster over polygons (zonal mean)
#'
#' @description
#' For each polygon in an sf object, compute the mean of all raster cells
#' from a terra SpatRaster that fall inside it. Returns the polygons with
#' a new column (or columns) containing those means.
#'
#' @param r A terra SpatRaster (one or more layers).
#' @param sf_polygons An sf POLYGON / MULTIPOLYGON object. Each row is a zone.
#' @param id_col Optional name of a column in sf_polygons that uniquely
#'   identifies polygons. If NULL, the function will generate an internal ID.
#' @param stats_col_prefix Prefix for the output columns. Default "mean_".
#'
#' @return The same sf_polygons with new columns like "mean_<layername>".
#'   Polygons outside the raster extent get NA for summary columns.
#'
#' @details
#' - Reprojects polygons to raster CRS if needed (using st_transform).
#' - Uses terra::extract(..., fun = mean, na.rm = TRUE).
#' - Handles multilayer rasters.
#' - Preserves all input polygons and their original row order; those outside
#'   raster extent get NA.
#'
#' @examples
#' \dontrun{
#' grid <- summarize_raster_mean_over_polygons(albedo_rast, target_grid)
#' }
#'
#' @export
summarize_raster_mean_over_polygons <- function(r,
                                                sf_polygons,
                                                id_col = NULL,
                                                stats_col_prefix = "mean_") {
  if (is.null(id_col)) {
    sf_polygons$.tmp_poly_id <- seq_len(nrow(sf_polygons))
    id_col <- ".tmp_poly_id"
    drop_tmp_id <- TRUE
  } else {
    if (!id_col %in% names(sf_polygons)) {
      stop("id_col '", id_col, "' not found in sf_polygons.")
    }
    drop_tmp_id <- FALSE
  }

  # Align CRS
  rast_wkt <- terra::crs(r, proj = TRUE)
  poly_crs <- sf::st_crs(sf_polygons)
  if (!is.null(poly_crs) && !is.na(poly_crs$wkt) && rast_wkt != "") {
    if (poly_crs$wkt != rast_wkt) {
      sf_polygons <- sf::st_transform(sf_polygons, crs = rast_wkt)
    }
  }

  poly_vect <- terra::vect(sf_polygons)

  extracted_df <- terra::extract(r, poly_vect, fun = mean, na.rm = TRUE)
  extracted_df[[id_col]] <- sf_polygons[[id_col]]

  layer_cols <- setdiff(names(extracted_df), c("ID", id_col))
  layer_names <- names(r)
  if (is.null(layer_names) || any(is.na(layer_names)) || any(layer_names == "")) {
    layer_names <- layer_cols
  }

  new_names <- paste0(stats_col_prefix, layer_names)
  names(extracted_df)[match(layer_cols, names(extracted_df))] <- new_names

  summary_df <- extracted_df[, c(id_col, new_names), drop = FALSE]

  # left_join preserves all original polygons; unmatched get NA
  sf_out <- dplyr::left_join(sf_polygons, summary_df, by = id_col)

  if (drop_tmp_id) {
    sf_out$.tmp_poly_id <- NULL
  }

  return(sf_out)
}


#' Summarize a categorical SpatRaster over an sf grid
#'
#' @description
#' For each cell in an sf grid, computes the count and proportion of fine-resolution
#' raster pixels belonging to each category, plus the dominant (modal) class.
#' Works by rasterizing the sf grid onto the SpatRaster's own grid, then
#' tabulating category membership — which is efficient when the raster is at
#' much finer resolution than the grid (e.g., 10m landcover over a 330m grid).
#'
#' @param SpatRaster_object A single-layer terra SpatRaster with integer/categorical values.
#' @param sf_grid An sf POLYGON object representing the target grid. Each row is one cell.
#'
#' @return `sf_grid` with additional columns:
#' \itemize{
#'   \item `n_<label>` — count of fine-resolution pixels of each class within the cell.
#'   \item `proportion_<label>` — proportion of pixels of each class within the cell.
#'   \item `dominant_class` — the class with the highest proportion in each cell.
#'   \item `pixel` — integer index identifying each grid cell (1-based row index).
#' }
#' Grid cells outside the raster extent are retained with NA landcover values.
#' Original row order is preserved.
#'
#' @details
#' Column names use the raw raster values as labels (e.g., `n_1`, `proportion_6`).
#' Use \code{\link{rename_and_fill_proportions}} afterward to map numeric codes
#' to human-readable names, and to fill implicit NAs with 0.
#'
#' @examples
#' \dontrun{
#' grid <- summarize_raster_class_representation_over_grid(landcover_rast, target_grid)
#' }
#'
#' @export
summarize_raster_class_representation_over_grid <- function(SpatRaster_object, sf_grid) {

  sf_grid$pixel <- seq_len(nrow(sf_grid))

  message("Converting sf object into SpatRaster...")
  grid_rast <- terra::rasterize(sf_grid, SpatRaster_object, field = "pixel")

  message("Creating tibble to summarize classes...")
  t <- dplyr::tibble(
    pixel = as.integer(terra::values(grid_rast)),
    label = as.integer(terra::values(SpatRaster_object))
  )

  # drop rows where either pixel or label is NA
  t <- dplyr::filter(t, !is.na(pixel), !is.na(label))

  missing_pixels <- !(sf_grid$pixel %in% t$pixel)
  number_of_missing_pixels <- sum(missing_pixels)
  if (number_of_missing_pixels > 0) {
    message(
      number_of_missing_pixels,
      " pixels in the sf object don't overlap with the SpatRaster object; ",
      "they will be retained with NA landcover values."
    )
  }

  message("Summarizing classes by count...")
  t <- dplyr::count(t, pixel, label)

  message("Summarizing classes by proportion...")
  t <- dplyr::group_by(t, pixel)
  t <- dplyr::mutate(t, proportion = n / sum(n))

  message("Computing dominant class...")
  dominant_class <- dplyr::summarize(t, dominant_class = label[which.max(proportion)])

  message("Adding everything to original sf object.")
  t <- tidyr::pivot_wider(t, names_from = label, values_from = c(n, proportion))

  # right_join onto sf_grid preserves all original grid rows;
  # st_as_sf() called after each join to prevent dplyr from dropping the sf class;
  # ungroup() ensures no residual grouping interferes with downstream sf operations;
  # sort by pixel to restore original row order
  sf_grid <- dplyr::right_join(t, sf_grid, by = "pixel")
  sf_grid <- sf::st_as_sf(sf_grid)
  sf_grid <- dplyr::left_join(sf_grid, dominant_class, by = "pixel")
  sf_grid <- sf::st_as_sf(sf_grid)
  sf_grid <- dplyr::ungroup(sf_grid)
  sf_grid <- sf_grid[order(sf_grid$pixel), ]

  return(sf_grid)
}
