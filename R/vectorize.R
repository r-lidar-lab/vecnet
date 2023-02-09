#' Vectorize a raster representing a segmented roads network
#'
#' Vectorize a raster representing a segmented roads network. The raster can be either a binary
#' map or a non binary map (see examples). The algorithm is expected to be robust to gaps, and
#' false positive (see references).
#'
#' @param seeds  \code{sfc} a set of seeds to start driving the network (see references and examples).
#' @param map raster a probability raster.  Either a binary map or a propability map the
#' former being a special case of the later (see references).
#' @param network a \code{sf/sfc} with an already existing network such as the algorithm
#' can stop when vectorizing an already vectorized part of the network.
#' @param max_sinuosity numeric. False positive may be detected and discarded based on the fact that
#' they have likely an irrelevantly high sinuosity.
#' @param min_length numeric. False positive may be detected and discarded based on the fact that
#' they have likely an very short length.
#' @param verbose boolean.
#' @param display boolean. display a map in pseudo real time of the advancement. For debbuging purpose mainly.
#' @param ... propagated to \link{track_line}
#'
#' @references
#' Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
#' (2022) Vectorial and topologically valid segmentation of forestry road networks from ALS data. In prep.\cr\cr
#' Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
#' (2022) Correction, update, and enhancement of vectorial forestry road maps using ALS data, a
#' pathfinder, and seven metrics. Journal of Applied Earth Observation and Geoinformation.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' library(sf)
#'
#' # =======================
#' # On a perfect binary map
#' # =======================
#'
#' # Generate a binary map from a vector for the example
#' network <- system.file("extdata", "network.shp", package = "vecnet")
#' network <- buffer(vect(network), 1)
#' map <- rast(network, resolution = 2)
#' map <- rasterize(network, map)
#' map[is.na(map)] = 0.1
#' seeds <- system.file("extdata", "seeds.shp", package = "vecnet")
#' seeds <- st_read(seeds) |> st_geometry()
#'
#' # Takes ~1 min to run
#' res <- vectorize_network(map, seeds, min_length = 150,  min_conductivity = 0.3, display = TRUE)
#'
#' # notice that some parts of the network are missing because they are not reachable
#' # with the seed used in this example.
#' col = sample(rainbow(length(res)))
#' plot(network, border = "gray")
#' plot(res, add = TRUE, col = col, lwd = 2)
#' starts = lwgeom::st_startpoint(res)
#' ends = lwgeom::st_endpoint(res)
#' plot(starts, add = TRUE, col = "black", pch = 19, cex = 0.5)
#' plot(ends, add = TRUE, col = "black", pch = 19, cex = 0.5)
#'
#' # =============================
#' # On an imperfect non binary map
#' # =============================
#'
#' map <- system.file("extdata", "network.tif", package = "vecnet")
#' seeds <- system.file("extdata", "seeds.shp", package = "vecnet")
#' map <- rast(map)
#' seeds <- st_read(seeds, quiet = TRUE) |> st_geometry()
#' plot(map, col = viridis::inferno(25), smooth = TRUE)
#' plot(seeds, add = T, col = "cyan", lwd = 4)
#'
#' res <- vectorize_network(map, seeds, min_conductivity = 0.6, display = TRUE)
#'
#' plot(map, col = viridis::inferno(25), smooth = TRUE, alpha = 0.5)
#' col = sample(rainbow(length(res)))
#' plot(res, add = TRUE, col = col, lwd = 2)
#' starts = lwgeom::st_startpoint(res)
#' ends = lwgeom::st_endpoint(res)
#' plot(starts, add = TRUE, col = "black", pch = 19, cex = 0.5)
#' plot(ends, add = TRUE, col = "black", pch = 19, cex = 0.5)
#' }
#' @export
vectorize_network = function(map, seeds, network = NULL, max_sinuosity = 1.8, min_length = 400, verbose = FALSE, display = FALSE, ...)
{
  level <- 1
  col <- c("#FF0000", "#FF6600", "#FFCC00", "#CCFF00", "#66FF00", "#00FF00",
           "#00FF66", "#00FFCC", "#00CCFF", "#0066FF", "#0000FF", "#6600FF",
           "#CC00FF", "#FF00CC", "#FF0066") #rainbow(15)

  seeds = sf::st_geometry(seeds)

  if (display)
    terra::plot(map, col = viridis::inferno(25), smooth = TRUE)

  if (is.null(network))
  {
    network = seeds[0]
    netseeds = seeds[0]
  }


  while (!is.null(seeds))
  {
    ans <- vector("list", length(seeds))
    for (i in seq_along(seeds))
    {
      tryCatch({
        res <- track_line(seeds[i], map, network = network, ...)
      },
      error = function(e) {
        seed = seeds[i]
        f1 = tempfile(fileext = ".gpkg")
        f2 = tempfile(fileext = ".gpkg")
        sf::st_write(seed, f1)
        sf::st_write(network, f2)
        message(e)
        stop("Internal error. \n seed logged at: ", f1, "\n network logged: ", f2, call. = FALSE)
      })

      len <- as.numeric(sf::st_length(res$road))
      sin <- sinuosity(res$road)
      den <- res$dintersection

      if (len > min_length & sin <= max_sinuosity)
      {
        if (display) plot(sf::st_geometry(res$road), add = TRUE, lwd = 2, col = col[level])
        if(!is.null(res$seeds) && display) plot(sf::st_geometry(res$seeds), add = TRUE, lwd = 1, col = "lightblue")
        ans[[i]] <- res
        network = c(network, res$road)
        netseeds = c(netseeds, seeds[i])
      }
      else
      {
        if (len <= min_length && display)
          plot(sf::st_geometry(res$road), add = T, lwd = 2, col = "gray")

        if (sin > max_sinuosity && display)
          plot(sf::st_geometry(res$road), add = T, lwd = 2, col = "pink")
      }
    }

    ans <- Filter(Negate(is.null), ans)
    if (length(ans) == 0) break;

    seeds <- lapply(ans, function(x) x$seeds)
    seeds <- Filter(Negate(is.null), seeds)
    seeds <- do.call(c, seeds)

    level <- level +1
  }

  network <- sf::st_as_sf(network)
  #return(list(network, netseeds))

  #if (display)
  #{
  #  plot(network, add = F, col = sample(rainbow(nrow(network))), lwd = 2)
  #  starts = lwgeom::st_startpoint(network)
  #  ends = lwgeom::st_endpoint(network)
  #  plot(starts, add = TRUE, col = "black", pch = 19, cex = 0.5)
  #  plot(ends, add = TRUE, col = "black", pch = 19, cex = 0.5)
  #}

  network <- sf::st_geometry(network)
  network <- sfnetworks::as_sfnetwork(network)
  #plot(network)

  network <- tidygraph::convert(network, sfnetworks::to_spatial_smooth)
  #plot(sf::st_geometry(network, "edges"), add = F, col = sample(rainbow(length(sf::st_geometry(network, "edges")))), lwd = 2)
  #plot(sf::st_geometry(network, "nodes"), add = TRUE, col = "black", pch = 19, cex = 0.5)

  network <- tidygraph::convert(network, sfnetworks::to_spatial_subdivision)
  #plot(sf::st_geometry(network, "edges"), add = F, col = sample(rainbow(length(sf::st_geometry(network, "edges")))), lwd = 2)
  #plot(sf::st_geometry(network, "nodes"), add = TRUE, col = "black", pch = 19, cex = 0.5)

  network <- sf::st_geometry(network, "edges")
  network <- sf::st_simplify(network, dTolerance = 2)

  return(network)
}

