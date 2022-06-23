#' Vectorize a raster representing a segmented roads network
#'
#' Vectorize a raster representing a segmented roads network. The raster can be either a binary
#' raster or a non binary raster (see examples). The algorithm is expected to be robust to gaps, and
#' false positive (see references)
#'
#' @param map
#' @param seeds
#' @param network
#' @param max_sinusoity
#' @param min_length
#' @param verbose
#' @param ...
#'
#' @examples
#' library(terra)
#' library(sf)
#'
#' # =======================
#' # On a perfect binary map
#' # =======================
#'
#' network <- system.file("extdata", "network.shp", package = "vecnet")
#' network = buffer(vect(network), 1)
#' map = rast(network, resolution = 2)
#' map = rasterize(network, map)
#' map[is.na(map)] = 0.1
#' seeds <- system.file("extdata", "seeds.shp", package = "vecnet")
#' seeds = st_read(seeds) |> st_geometry()
#'
#' res <- vectorize_network(map, seeds, min_length = 0)
#'
#' plot(map, col = viridis::inferno(25), smooth = TRUE)
#' col = sample(rainbow(nrow(res)))
#' plot(res, add = TRUE, col = col, lwd = 2)
#' starts = lwgeom::st_startpoint(res)
#' ends = lwgeom::st_endpoint(res)
#' plot(starts, add = TRUE, col = "black", pch = 19, cex = 0.5)
#' plot(ends, add = TRUE, col = "black", pch = 19, cex = 0.5)
#'
#' # =============================
#' # On a imperfect non binary map
#' # =============================
#'
#' map <- system.file("extdata", "network.tif", package = "vecnet")
#' seeds <- system.file("extdata", "seeds.shp", package = "vecnet")
#' map = rast(map)
#' seeds = st_read(seeds) |> st_geometry()
#' plot(map, col = viridis::inferno(25), smooth = TRUE)
#' plot(seeds, add = T, col = "cyan", lwd = 4)
#'
#' res <- vectorize_network(map, seeds, min_conductivity = 0.55)
#'
#' plot(map, col = viridis::inferno(25), smooth = TRUE, alpha = 0.5)
#' col = sample(rainbow(nrow(res)))
#' plot(res, add = TRUE, col = col, lwd = 2)
#' starts = lwgeom::st_startpoint(res)
#' ends = lwgeom::st_endpoint(res)
#' plot(starts, add = TRUE, col = "black", pch = 19, cex = 0.5)
#' plot(ends, add = TRUE, col = "black", pch = 19, cex = 0.5)
#' @export
vectorize_network = function(map, seeds, network = NULL, max_sinusoity = 2.6, min_length = 400, verbose = FALSE, ...)
{
  level <- 1
  col <- rainbow(8)

  if (is.null(network))
    network = seeds[0]

  while (!is.null(seeds))
  {
    ans <- vector("list", length(seeds))
    for (i in seq_along(seeds))
    {
      tryCatch({
        res <- track_line(seeds[i], map, network = network, min_conductivity = 0.6)#, ...)
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
      sin <- ALSroads:::sinuosity.sfc_LINESTRING(res$road)
      den <- res$dintersection

      if (len > min_length & sin <= max_sinusoity)
      {
        plot(st_geometry(res$road), add = TRUE, lwd = 2, col = col[level])
        if(!is.null(res$seeds)) plot(st_geometry(res$seeds), add = TRUE, lwd = 1, col = "lightblue")
        ans[[i]] <- res
        network = c(network, res$road)
      }
      else
      {
        if (len <= min_length)
          plot(st_geometry(res$road), add = T, lwd = 2, col = "gray")

        if (sin > max_sinusoity)
          plot(st_geometry(res$road), add = T, lwd = 2, col = "pink")
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
  #plot(network, add = F, col = sample(rainbow(nrow(network))), lwd = 2)
  #starts = lwgeom::st_startpoint(network)
  #ends = lwgeom::st_endpoint(network)
  #plot(starts, add = TRUE, col = "black", pch = 19, cex = 0.5)
  #plot(ends, add = TRUE, col = "black", pch = 19, cex = 0.5)

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
  #network <- sf::st_simplify(network, dTolerance = 2)
  return(network)
}

