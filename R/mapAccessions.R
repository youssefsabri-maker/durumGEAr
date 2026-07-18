#' Map Accession or Unit Locations
#'
#' A lightweight base-graphics scatter of collection locations, optionally
#' coloured by a grouping variable (e.g. Country or Cluster). Mirrors the
#' \code{mapAccessions} utility of the \pkg{icardaFIGSr} package but avoids any
#' heavy mapping dependency, so it works in a minimal environment.
#'
#' @param data A data frame with longitude and latitude columns.
#' @param lon,lat Column names for coordinates. Defaults \code{"Longitude"},
#'   \code{"Latitude"}.
#' @param color_by Optional column name to colour points by. Default \code{NULL}.
#' @param main Plot title.
#' @param pch,cex Point style and size.
#' @param legend Logical; draw a legend when \code{color_by} is set and has
#'   fewer than 25 levels. Default TRUE.
#'
#' @return Invisibly, the data frame used for plotting. Called for its side
#'   effect (a plot).
#'
#' @examples
#' data(durumUnits)
#' mapAccessions(durumUnits, color_by = "Cluster",
#'               main = "Durum wheat effective units")
#'
#' @seealso \code{\link{collapseUnits}}
#' @export
mapAccessions <- function(data, lon = "Longitude", lat = "Latitude",
                          color_by = NULL, main = "Accession locations",
                          pch = 19, cex = 0.6, legend = TRUE) {
  stopifnot(all(c(lon, lat) %in% names(data)))
  x <- data[[lon]]; y <- data[[lat]]
  if (!is.null(color_by) && color_by %in% names(data)) {
    g <- factor(data[[color_by]])
    pal <- grDevices::hcl.colors(nlevels(g), palette = "Dark 3")
    col <- pal[as.integer(g)]
  } else { g <- NULL; col <- "steelblue" }
  graphics::plot(x, y, col = col, pch = pch, cex = cex,
                 xlab = "Longitude", ylab = "Latitude", main = main)
  graphics::grid(col = "grey90")
  if (!is.null(g) && legend && nlevels(g) < 25)
    graphics::legend("topright", legend = levels(g),
                     col = grDevices::hcl.colors(nlevels(g), palette = "Dark 3"),
                     pch = pch, cex = 0.7, bg = "white", ncol = 2)
  invisible(data)
}
