#' Find seeds for the vectorization
#'
#' Drive along a contour line of the raster using \link{track_line} to find seeds that will be used
#' to vectorized the network
#'
#' @param map SpatRaster. Raster of probabilities/conductivities. Either a binary map or a probability
#' map the former being a special case of the later (see references)
#' @param contour an sfc_POLYGON that correspond to the contour of the raster. If the raster is
#' rectangular the default corresponds to the bounding box.
#' @param ... propagated to \link{track_line}
#'
#' @references
#' Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
#' (2022) Vectorial and topologically valid segmentation of forestry road networks from ALS data. In prep.\cr\cr
#' Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
#' (2022) Correction, update, and enhancement of vectorial forestry road maps using ALS data, a
#' pathfinder, and seven metrics. Journal of Applied Earth Observation and Geoinformation.
#'
#' @export
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' map <- system.file("extdata", "network2.tif", package = "vecnet")
#' map <- rast(map)
#'
#' seeds = init_seeds(map, min_conductivity = 0.8)
#'
#' plot(map, col = viridis::inferno(25))
#' plot(seeds, add = TRUE, col = "green", lwd = 3)
#' }
init_seeds = function(map, contour = sf::st_as_sfc(sf::st_bbox(map)), ...)
{
  network = sf::st_buffer(contour, dist = -75)
  polygon = sf::st_difference(contour, network)
  network = sf::st_cast(network, "LINESTRING")
  bb = sf::st_cast(network, "POINT")
  seed = sf::st_linestring(sf::st_coordinates(bb[1:2]))
  l = sf::st_length(seed)
  r = 10/l
  seed = lwgeom::st_linesubstring(seed, 0, r)
  seed = sf::st_sfc(seed)
  sf::st_crs(seed) <- sf::st_crs(map)

  ans = track_line(seed, map, network, sightline = 50, seed_mode = TRUE, ...)

  return(ans$seeds)
}

