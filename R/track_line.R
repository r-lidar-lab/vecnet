#' Drive and vectorize an unknown road
#'
#' Drive along an unknown road in a probability map starting from small piece
#' of road segment pointing in the right direction and vectorize one road by returning a spatial
#' line. This function is used by \link{vectorize_network} and is not intended to be used
#' by regular users except for trial-and-error tests.
#'
#' @param seed \code{sfc_LINESTRING}. A seed to start driving a road
#' @param conductivity \code{SpatRaster}. A raster with value ranging in [0,1] where roads pixels are
#' close to 1 and non road pixels are close to 0 (see references).
#' @param network \code{sf/sfc} with LINSTRING. An already existing network such as the algorithm
#' can stop when vectorizing an already vectorized part of the network. Not mandatory. It is used
#' automatically by \link{vectorize_network}.
#' @param fov  numeric. Field of view (degrees) ahead of the search vector (see references).
#' @param sightline  numeric (distance unit). Search distance used to find the next most probable
#' point on the road (see references).
#' @param min_conductivity numeric. Between 0 and 1. Corresponds to the sensitivity of the method. A
#' value close to 1 indicates that the algorithm follows only the pixels with a very high
#' conductivity and stops easily when the average conductivity is lower than this value. With an high
#' value the algorithm may miss roads. A low value indicates that the algorithm follows the pixel even
#' with low conductivity and is likely to vectorize roads that do not exist but is less likely to miss
#' existing roads.
#' @param th_conductivity numeric. Between 0 and 1. A conductivity map with too low values is not allowed
#' (see reference). All the values lower than this value are clamped and replaced by this value. The default
#' is 0.1 meaning that every values from \code{conductivity} lower than 0.1 is set to 0.1
#' @param partial_gap_size numeric. The algorithm can drive on low conductivity segments even if the
#' conductivity is too low and should trigger a stop signal if it looks like a sharp line and an
#' high conductivity road is reach after the gap (see references). The allowed distance is a multiple
#' of the sightline. Default is 2.5 which roughly means that it can follow a too low conductivity
#' track for two iterations and at the third one it must be back on a high conductivity track.
#' @param ... Unused
#' @param disp bool. Display in realtime the progress on images. For debugging purposes.
#'
#' @return list.
#' \enumerate{
#' \item{\code{road} contains the road found in a vectorial format (\code{sfc}) }
#' \item{\code{seed} conntains all the potential intersection and seeds found along
#'  the road in a vectorial format (\code{sfc}) }}
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
#' map <- system.file("extdata", "network.tif", package = "vecnet")
#' seeds <- system.file("extdata", "seeds.shp", package = "vecnet")
#' map <- rast(map)
#' seeds <- st_read(seeds, quiet = TRUE)
#' seed <- st_geometry(seeds)[1]
#'
#' res <- track_line(seed, map, min_conductivity = 0.5)
#'
#' plot(map, col = viridis::inferno(25))
#' plot(seed, add = TRUE, col = "green", lwd = 3)
#' plot(res$road, add = TRUE, col = "red", lwd = 2)
#' plot(res$seeds, add = TRUE, col = "green", lwd = 3)
#' }
track_line <- function(seed,
                       conductivity,
                       network = NULL,
                       fov = 160,
                       sightline = 100,
                       min_conductivity = 0.6,
                       th_conductivity = 0.1,
                       partial_gap_size = 2.5,
                       ...,
                       disp = FALSE)
{
  # network = NULL ; fov = 160 ; sightline = 100 ; min_conductivity = 0.6 ; th_conductivity = 0.1 ; partial_gap_size = 2.5 ;
  t0 <- Sys.time()

  if (methods::is(seed, "sf")) seed = sf::st_geometry(seed)
  if (!methods::is(seed, "sfc_LINESTRING")) stop("seed must be sfc_LINESTRING")
  if (length(seed) != 1L) stop("seed must be of length 1")
  if (!methods::is(conductivity, "SpatRaster")) stop("conductivity must be a SpatRaster")

  if (!is.null(network)) network = sf::st_geometry(network)

  # We need an orthogonal raster
  resolution <- terra::res(conductivity)
  if (resolution[1] != resolution[2]) stop("'conductivity' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  # Retains a trace of already explored areas
  trace <- sf::st_geometry(seed)

  # Because we want an "overall" direction and we want to drop
  # smooth little imperfections in the seed if any.
  seed <- sf::st_simplify(seed, dTolerance = 5)

  # Threshold of cost to stop the search.
  cost_max <- sightline * 1/min_conductivity

  # Initialization of list of coordinates with starting line
  start_pts <- sf::st_cast(sf::st_geometry(seed), "POINT")
  n <- length(start_pts)
  start <- start_pts[n-1][[1]]
  end <- start_pts[n][[1]]
  #list_coords <- list(start, end)
  list_lines <- list(sf::st_geometry(seed)[[1]])
  list_intersections <- list()
  heading <- get_heading(start, end)

  sub_aoi_conductivity <- query_conductivity_aoi(conductivity, seed, smooth = T)
  sub_aoi_conductivity[sub_aoi_conductivity < th_conductivity] = th_conductivity

  #plot(conductivity, col = viridis::inferno(50))
  #plot(terra::ext(sub_aoi_conductivity), col = "red", add = T, alpha = 0.3)
  #plot(seed, add = T, lwd = 3, col = "blue")
  #plot(sub_aoi_conductivity, col = viridis::inferno(50))

  # Mask the existing road to avoid driving a known road
  dots = list(...)
  mask_value = 0
  buffer = 10
  find_seek_mask = NULL
  if (!is.null(dots$find_seed_mode))
  {
    mask_value = 2
    buffer = 2
    find_seek_mask = terra::vect(dots$find_seed_mode)
  }

  sub_aoi_conductivity <- mask_existing_network(sub_aoi_conductivity, network, mask_value, buffer)

  if (!is.null(find_seek_mask))
    sub_aoi_conductivity = terra::mask(sub_aoi_conductivity, find_seek_mask, inverse = TRUE, updatevalue = th_conductivity)

  # Init view angles as a function of the resolution of the raster and the sightline
  angles_rad <- generate_angles(resolution, sightline, fov)

  # Loop initialization
  current_cost <- 0
  k <- 2
  n <- 1
  overcost = 0
  novercost = 0
  dovercost = 0

  while (dovercost <= partial_gap_size*sightline & !is.infinite(cost_max))
  {
    if (k %% 5 == 0)
    {
      dist = (k-1)*sightline*0.8
      cat(dist, " m (", get_speed(dist, t0), " km/h)\r", sep = "")
    }
    utils::flush.console()

    # Compute heading from the previous segment
    l <- list_lines[[k-1]]
    l <- lwgeom::st_linesubstring(l, 0.5, 1)
    p1 <- lwgeom::st_startpoint(l)[[1]]
    p2 <- lwgeom::st_endpoint(l)[[1]]
    heading <- get_heading(p1, p2)

    # Create all possible ends ahead of the previous points
    start <- sf::st_sfc(p2)
    ends <- generate_ends(p2, angles_rad, sightline, heading)
    ends$angle = angles_rad

    # Generate a very small area of interest corresponding to what can be seen in the line of sight
    search_zone <- sf::st_bbox(c(start, sf::st_geometry(ends)))
    search_zone <- terra::ext(c(search_zone[1], search_zone[3],search_zone[2],search_zone[4])) + 10
    aoi <- terra::crop(sub_aoi_conductivity, search_zone)
    #aoi <- terra::stretch(aoi, maxv = 1)
    p <- sf::st_polygon(list(rbind(sf::st_coordinates(start), sf::st_coordinates(ends), sf::st_coordinates(start))))
    p <- sf::st_sfc(p)
    p <- sf::st_buffer(p, 2)
    p <- terra::vect(p)
    terra::crs(p) <- terra::crs(aoi)
    aoi <- terra::mask(aoi, p)
    #plot(aoi, col = viridis::inferno(50), range = c(0,1))
    #plot(start,add =T, col = "red", pch = 19)
    #plot(ends$geometry,add =T, col = "red", pch = 19)

    # Compute the transition matrix for this AOI
    trans <- transition(aoi, geocorrection = TRUE)

    # Check if some ends fall outside of the AOI
    # If any we  reached the loaded part of the conductivity.
    # Reload another raster part further ahead.
    val <- terra::extract(aoi, sf::st_coordinates(ends))[[1]]
    if (any(is.na(val)))
    {
      # Make a newly sub aoi conductivity raster further ahead
      line <- sf::st_cast(c(p1, p2), "LINESTRING") |> sf::st_sfc()
      line <- sf::st_set_crs(line, sf::st_crs(seed))

      sub_aoi_conductivity <- query_conductivity_aoi(conductivity, line, smooth = T)
      sub_aoi_conductivity[sub_aoi_conductivity < th_conductivity] = th_conductivity
      n = n+1
      #terra::plot(conductivity, col = viridis::inferno(50))
      #terra::plot(terra::ext(sub_aoi_conductivity), col = "red", add = T)
      #plot(seed, add = T, lwd = 3, col = "red")

      sub_aoi_conductivity <- mask_existing_network(sub_aoi_conductivity, network, mask_value, buffer)
      sub_aoi_conductivity <- mask_passage(sub_aoi_conductivity, trace, 0, 5, sf::st_crs(seed))

      if (!is.null(find_seek_mask))
        sub_aoi_conductivity = terra::mask(sub_aoi_conductivity, find_seek_mask, inverse = TRUE, updatevalue = th_conductivity)

      # Check again if some ends fall outside of the newly cropped conductivity raster
      # If yes we are close to the edge of the raster. Try again with half the sightline
      val <- terra::extract(sub_aoi_conductivity, sf::st_coordinates(ends))[[1]]
      if (any(is.na(val)))
      {
        tmp_angles_rad <- generate_angles(resolution, sightline/2, fov)
        ends <- generate_ends(p2, tmp_angles_rad, sightline/2, heading)
        ends$angle <- tmp_angles_rad
      }

      # If we still have NA we reached the border of the conductivity map
      val <- terra::extract(sub_aoi_conductivity, sf::st_coordinates(ends))[[1]]
      if (any(is.na(val)))
      {
        message("Drive stopped early. Edge of conductivity raster has been reached.")
        cost_max <- -Inf
        break
      }

      # Re-extract the AOI and recompute the transition
      aoi <- terra::crop(sub_aoi_conductivity, search_zone)
      #aoi <- terra::stretch(aoi, maxv = 1)
      trans <- transition(aoi, geocorrection = TRUE)
    }

    # Estimates the cost to reach each end points and returns the main direction and
    # other possible directions (i.e. intersections)
    ans <- find_reachable(start, ends, trans, cost_max)

    # The only way to get infinite values is to reach an existing road
    # We need to reduce the sightline to reduce down the approach speed
    j = 0
    while (any(ans$cost >= 9999) & sightline/(2^j) > 25)
    {
      j = j+1
      tmp_angles_rad <- generate_angles(resolution, sightline/(2^j), min(fov, 90))
      ends <- generate_ends(p2, tmp_angles_rad, sightline/(2^j), heading)
      ends$angle <- tmp_angles_rad
      ans <- find_reachable(start, ends, trans, cost_max/(2^j))
    }

    if (any(ans$cost >= 9999))
    {
      message("Driving stopped because it reached another road")
      cost_max = -Inf
      break
    }

    if (is.null(ans))
    {
      warning("Driving stopped because not reachable point have been found", call. = FALSE)
      cost_max = -Inf
      break
    }

    idx_main <- ans$minima$idx[1]
    depth_main = ans$minima$depth[1]
    idx_other <- ans$minima$idx[-1]
    cost <- ans$cost
    releative_cost <- ans$minima$rcost
    current_cost <- releative_cost[1]

    end <- ends[idx_main,]
    ends_other <- ends[idx_other,]
    cost_other <- releative_cost[-1]

    delta_angle = abs(end$angle - ends_other$angle)*180/pi
    ends_other = ends_other[delta_angle > 15,]
    cost_other = cost_other[delta_angle > 15]

    L <- gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(end), output = "SpatialLines")
    L <- sf::st_geometry(sf::st_as_sf(L))
    C <- sf::st_coordinates(L)
    C <- C[1:(nrow(C)*0.8),-3]
    L <- sf::st_sfc(sf::st_linestring(C))
    end <- lwgeom::st_endpoint(L)

    #list_coords[[k+1]] <- sf::st_geometry(end)[[1]]
    list_lines[[k]] <- L[[1]]
    trace <- st_join_linestring(trace, L)

    # Protect against infinite loops
    sub_aoi_conductivity <- mask_passage(sub_aoi_conductivity, trace, 0, 5,sf::st_crs(seed))
    #plot(sub_aoi_conductivity, col = viridis::inferno(50))

    if (current_cost > 1)
    {
      novercost = novercost + 1
      dovercost = dovercost + sightline
      overcost = overcost + (current_cost - cost_max)
    }
    else
    {
      novercost = 0
      dovercost = 0
      overcost = 0
    }

    I <- NULL
    if (length(sf::st_geometry(ends_other)) > 0)
    {
      I <- gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(ends_other), output = "SpatialLines") |> suppressWarnings()
      I <- sf::st_geometry(sf::st_as_sf(I))
      I <- sf::st_difference(I,L)

      I = Filter(function(x){
        if (sf::st_geometry_type(x) != "LINESTRING")
          return(FALSE)

        if (lwgeom::st_endpoint(L) ==  lwgeom::st_startpoint(x))
          return(FALSE)

        return(TRUE)
      }, I)

      if (length(I) > 0)
      {
        list_intersections <- c(list_intersections, list(I))
        sub_aoi_conductivity <- mask_passage(sub_aoi_conductivity, I, 5, 0, sf::st_crs(seed))
      }
    }

    k <- k + 1

    if (disp)
    {

      if (length(angles_rad) == length(cost))
      {
        #svglite::svglite("/home/jr/Documents/Ulaval/2022 Post-doc/Articles/vectnet/img/svg/followroadcostprofile.svg", height = 8, width = 8)
        graphics::par(family = "serif")
        plot(angles_rad*180/pi, cost, type = "l", lwd = 2, ylim = c(sightline, max(cost) *1.1), col = "red", xlab = "Angle (\u00B0)", ylab = "Cost (1)")
        smooth <- 9
        scost = ma(cost, n = smooth)
        scost = ma(scost, n = 5)
        graphics::lines(angles_rad*180/pi, scost, col = "darkgreen", lwd = 2)
        graphics::abline(v = angles_rad[c(idx_main, idx_other)]*180/pi, col = "black",lwd = 2)
        graphics::abline(h = cost_max, lty = 3)
        graphics::abline(h = 2.5*cost_max, lty = 3)
        graphics::legend(80, 990, legend=c("Cost profile", "Smoothed profile"),
               col=c("red", "darkgreen"), lty=1:2, cex=0.8)
        #abline(v = angles_rad[idx_main]*180/pi + c(-0.2612, 0.2612)*180/pi, lty = 3)
        #dev.off()
      }
      terra::plot(terra::crop(sub_aoi_conductivity, terra::ext(aoi) + 80),
                  col = viridis::inferno(50),
                  main = paste0(k,"/", n),
                  range = c(0,1),
                  axes = FALSE)
      terra::plot(aoi, add = T, col = viridis::inferno(50))
      terra::plot(terra::ext(aoi), add = T)
      ends$cost = cost
      tryCatch({
        plot(ends["cost"], pal = viridis::viridis, pch = 19, add = T, cex = 0.25, breaks = "quantile")
      }, error = function(x) {
        plot(ends$geometry, col = "red", add = T, pch = 19, cex = 0.5)
      })

      tryCatch({
        paths = gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(ends), output = "SpatialLines") |> suppressWarnings()
        paths = sf::st_as_sf(paths)
        paths$cost = ends$cost
        plot(paths, add = T,  pal = viridis::viridis)
        plot(p1, col = "red", pch = 19, add = T)
        plot(p2, col = "red", pch = 19, add = T)
      }, error = function(x) print(x))

      plot(L, add = T, col = "red", lwd = 4)
      if(!is.null(I)) plot(I, add = T, col = "blue", lwd = 3)
      #if(!is.null(I2)) plot(I2, add = T, col = "cornflowerblue", lwd = 3)
      plot(end, add = T, col = "green", pch = 19, cex = 2)
      cat("step = ", k-1, "=====\n")
      cat("cost =", current_cost, " | ")
      cat("overcost =", overcost, " | ")
      cat("dovercost =", dovercost, "\n\n")

      In = I
      #tm = paper_figure()
    }
  }

  list_lines <- list_lines[1:(length(list_lines)-novercost)]
  first = list_lines[[1]]
  list_lines <- lapply(list_lines, function(x) { x <- x[] ; return(x[-1,])})
  list_lines[[1]] = first

  newline <- do.call(rbind, list_lines) |>
    sf::st_linestring() |>
    sf::st_sfc() |>
    sf::st_set_crs(sf::st_crs(seed))


  if (!is.null(network) && length(network) >= 1)
  {
    tail_point = lwgeom::st_endpoint(newline)
    distance_to_network = min(sf::st_distance(network, tail_point))
    if (as.numeric(distance_to_network) < 75)
    {
      u = sf::st_simplify(newline, dTolerance = 2)
      u = lwgeom::st_linesubstring(u, 0.5, 1)

      #plot(u, xlim = st_bbox(u)[c(1, 3)] + c(-20,20),  ylim = st_bbox(u)[c(2, 4)] + c(-20,20))
      #plot(network, add = T)
      u = st_extend_line(u, 80, end = "TAIL" )
      #plot(u, xlim = st_bbox(u)[c(1, 3)] + c(-20,20),  ylim = st_bbox(u)[c(2, 4)] + c(-20,20))
      #plot(network, add = T)
      p = sf::st_intersection(u, network)
      if (length(p) == 1 && sf::st_geometry_type(p) == "POINT")
      {
        M = rbind(sf::st_coordinates(newline)[,1:2], sf::st_coordinates(p)[,1:2])
        newline = sf::st_linestring(M) |> sf::st_sfc() |> sf::st_set_crs(sf::st_crs(seed))
        cat("End of the road connected to the existing network.\n")
      }

      #plot(newline, add = T, col = "red", lwd = 2)
      #plot(network, add = T)
    }
  }

  if (length(list_intersections) >= 1)
  {
    intersections <- do.call(c, list_intersections) |>
      sf::st_set_crs(sf::st_crs(seed))

    p = lwgeom::st_startpoint(intersections)
    d = as.numeric(sf::st_distance(p,newline))
    intersections <- intersections[d < 1]
  }
  else
  {
    intersections <- NULL
  }


  nintersection <- length(intersections)
  len = as.numeric(sf::st_length(newline))
  dintersection = nintersection/(len/1000)

  h <- km <- NULL
  tf <- Sys.time()
  dt <- tf-t0
  dt <- round(units::as_units(dt),1)
  dist <- sf::st_length(newline)
  hdt = dt
  units(hdt) <- units::make_units(h)
  units(dist) <- units::make_units(km)
  dist = round(dist, 2)
  speed = round(dist/hdt)
  cat("Processed ended in", dt, units::deparse_unit(dt),
      ": road of", dist, units::deparse_unit(dist),
      "driven at", speed, units::deparse_unit(speed),
      "\n")

  return(list(road = newline, seeds = intersections, dintersection = dintersection))
}

# Very similar to st_ends_heading()
get_heading <- function(from, to)
{
  from <- sf::st_coordinates(from)
  to   <- sf::st_coordinates(to)
  M <- rbind(from, to)
  Ax <- M[2,1]
  Ay <- M[2,2]
  Bx <- M[1,1]
  By <- M[1,2]
  heading <- atan2(Ay-By, Ax-Bx)
}

query_conductivity_aoi <- function(conductivity, seed, smooth = TRUE)
{
  resolution <- terra::res(conductivity)[1]

  buf_dist <- c("ahead" = 900, "behind" = 100, "side" = 400)  # Could be set as a parameter

  aoi <- seed |>
    st_extend_line(c(buf_dist["ahead"], buf_dist["behind"])) |>
    sf::st_buffer(buf_dist["side"], endCapStyle = "FLAT")
  #plot(aoi, add = T, col= "red")
  #plot(seed, add = T, lwd = 3)

  conductivity_crop <- terra::crop(conductivity, terra::vect(aoi))
  #plot(conductivity_crop, col = viridis::inferno(50))
  #plot(seed, add = T, lwd = 3, col = "red")

  if (smooth)
  {
    w = matrix(1,3,3)
    w[5] = 10
    w = w/sum(w)
    conductivity_crop = terra::focal(conductivity_crop, w = w, fun = "mean")
  }

  conductivity_crop[is.na(conductivity_crop)] <- 0

  return(conductivity_crop)
}

generate_angles <- function(resolution, sightline, fov)
{
  angular_resolution <- floor(1.5 * (resolution / sightline) * (180 / pi))
  angles <- seq(angular_resolution, fov, angular_resolution)
  angles <- c(rev(-angles), 0, angles)
  angles_rad <- angles * pi / 180
  angles_rad
}

generate_ends <- function(center, angles, sightline, direction)
{
  center <- sf::st_coordinates(center)
  X <- center[1] + sightline * cos(direction + angles)
  Y <- center[2] + sightline * sin(direction + angles)
  M <- data.frame(X,Y)
  sf::st_as_sf(M, coords = c("X", "Y"))
}

find_reachable <- function(start, ends, trans, cost_max)
{
  if (length(sf::st_geometry(ends)) < 6) return(NULL)

  cost <- gdistance::costDistance(trans, sf::st_coordinates(start), sf::st_coordinates(ends))
  cost <- as.numeric(cost)

  if (all(is.infinite(cost))) return(NULL)

  cost[is.infinite(cost)] <- max(9999, max(cost[!is.infinite(cost)]))

  # Select local minima of cost
  smooth <- 9
  scost = -ma(cost, n = smooth)
  scost = ma(scost, n = 5)
  offset = floor(smooth/2) + floor(5/2)
  scost = as.numeric(stats::na.omit(scost))
  minima = pracma::findpeaks(scost, zero = "+")

  if (is.null(minima)) return(NULL)
  if (nrow(minima) == 0) return(NULL)

  depth = minima[,1]
  center = minima[,2] + offset
  left = minima[,3] + offset
  right = minima[,4] + offset
  cost2 = depth

  #plot(seq_along(cost), cost, type = "l", col = "red", log = "y")
  #lines(-scost, type = "l", col = "blue")
  #abline(v = center)
  #abline(v = left, col = "red")
  #abline(v = right, col = "blue")

  for (i in seq_along(center))
  {
    lhs = cost[left[i]]
    rhs = cost[right[i]]
    center[i] = which.min(cost[left[i]:right[i]]) + left[i] - 1
  }

  for (i in seq_along(center))
  {
    lhs = cost[left[i]]
    rhs = cost[right[i]]
    cost2[i] = cost[center[i]]
    depth[i] = mean(c(lhs, rhs)) - cost2[i]
  }

  minima = data.frame(idx = center, depth = depth, cost = cost2)

  depth_thresolds = c(50,400)
  factor_thresold = c(1,2.5)
  slope = diff(factor_thresold)/diff(depth_thresolds)
  intercept = factor_thresold[1] - slope * depth_thresolds[1]
  cost_max_multiplier_for_depth = minima$depth*slope + intercept
  cost_max_multiplier_for_depth[cost_max_multiplier_for_depth > factor_thresold[2]] <- factor_thresold[2]
  minima$rcost = minima$cost/(cost_max_multiplier_for_depth*cost_max)
  keep = (minima$rcost <= 1)

  if (sum(keep) > 0)
    minima = minima[keep,]
  else
    minima = minima[which.min(minima$cost),]


  data.table::setorder(minima, cost)

  return(list(minima = minima, cost = cost))
}

mask_existing_network <- function(x, network, updatevalue = 0, buffer = 10)
{
  if (!is.null(network) && length(network) > 0)
  {
    bb <- sf::st_bbox(x)
    bb <- sf::st_set_crs(bb, sf::st_crs(network))
    mask <- sf::st_crop(network, bb)
    if (length(mask) > 0)
    {
      mask <- sf::st_buffer(mask, dist = buffer)
      x <- terra::mask(x, terra::vect(mask), inverse = TRUE, updatevalue = updatevalue)
    }
  }

  return(x)
}

mask_passage <- function(raster, lines, start_cut, end_cut, crs)
{
  len = as.numeric(sf::st_length(lines))
  from = 1 - (len-start_cut)/len
  to = (len-end_cut)/len
  if (length(from) > 1)
  {
    from = mean(from)
    to = mean(to)
  }
  if (methods::is(lines, "sfg")) lines <- sf::st_sfc(lines)
  mask <- lwgeom::st_linesubstring(lines, from, to)
  mask <- sf::st_buffer(mask, dist = 8, endCapStyle = "FLAT")
  mask <- sf::st_set_crs(mask, crs)
  raster <- terra::mask(raster, terra::vect(mask), inverse = TRUE, updatevalue = 0.05)
  return(raster)
}

get_speed = function(dist, t0)
{
  h <- m <- km <- NULL
  tf <- Sys.time()
  dt <- tf-t0
  dt <- round(units::as_units(dt),1)
  hdt = dt
  units(hdt) <- units::make_units(h)
  units(dist) <- units::make_units(m)
  dist = round(dist, 2)
  speed = round(dist/hdt)
  units(speed) <- units::make_units(km/h)
  speed = round(speed, 0)
  return(speed)
}

# paper_figure = function()
# {
#   library(tmap)
#   crs = st_crs(sub_aoi_conductivity)
#   A = ends["cost"] |> st_set_crs(crs)
#   P = gdistance::shortestPath(trans, sf::st_coordinates(start), sf::st_coordinates(ends), output = "SpatialLines") |> suppressWarnings()
#   P = sf::st_as_sf(P)|> st_set_crs(crs)
#   P$cost = ends$cost
#   a = st_sfc(p1) |> st_set_crs(crs)
#   b = st_sfc(p2) |> st_set_crs(crs)
#   L = L |> st_set_crs(crs)
#   end = end |> st_set_crs(crs)
#
#   tm = tm_shape(terra::crop(sub_aoi_conductivity, terra::ext(aoi) + 80)) +
#     tm_raster(palette = viridis::inferno(25), n = 15, style = "cont", title = "Probability",  breaks = c(0,0.25,0.5,0.75,1)) +
#     tm_shape(aoi) +
#     tm_raster(palette = viridis::inferno(25), n = 15, style = "cont", legend.show = F) +
#     tm_shape(A) + tm_dots("cost", pal = viridis::viridis(25), size = 0.05, style = "cont") +
#     tm_shape(P) + tm_lines("cost", pal = viridis::viridis(25), legend.col.show = FALSE) +
#     #tm_shape(a) + tm_dots(col = "red", shape = 19, size = 1) +
#     tm_shape(b) + tm_dots(col = "red", shape = 19, size = 1) +
#     tm_shape(L) + tm_lines(col = "red", lwd = 5) +
#     tm_shape(end) + tm_dots(col = "green", size = 1, shape = 17) +
#     tm_layout(legend.bg.color = "white", fontfamily = "serif") +
#     tm_add_legend("line", label = c("Road found", "Intersection found"), col = c("red", "blue"), lwd = 5) +
#     tm_scale_bar(text.color = "white", text.size = 1, breaks = c(0,0.1), width = 0.1)
#
#   if (!is.null(In))
#   {
#     In = In |> st_set_crs(crs)
#     tm = tm + tm_shape(In) + tm_lines(col = "blue", lwd = 4)
#   }
#
#   tm
# }
