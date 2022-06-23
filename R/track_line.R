#' Drive along an unknown road
#'
#' Drive along an unknown road in a conductivity raster starting from small piece
#' of road segment pointing in the right direction.
#'
#' @param seed  \code{sfc_LINESTRING} a seed to start driving a road
#' @param conductivity  raster (\code{raster} format)
#' @param fov  numeric. Field of view (degrees) ahead of the search vector.
#' @param sightline  numeric (distance unit). Search distance used to find the next most probable
#' point on the road.
#'
#' @return list, \code{line} being the road found (as \code{sfc}) and \code{cost} a numeric vector
#' of cost at each invidual vertices of the line.
#' @export
#' @examples
#' library(terra)
#' library(sf)
#'
#' map <- system.file("extdata", "network.tif", package = "vecnet")
#' seeds <- system.file("extdata", "seeds.shp", package = "vecnet")
#' map = rast(map)
#' seeds = st_read(seeds)
#' seed = st_geometry(seeds)[1]
#'
#' res <- track_line(seed, map, min_conductivity = 0.5)
#'
#' plot(map, col = viridis::inferno(25))
#' plot(seed, add = TRUE, col = "green", lwd = 3)
#' plot(res$road, add = TRUE, col = "red", lwd = 2)
#' plot(res$seeds, add = TRUE, col = "green", lwd = 3)
track_line <- function(seed, conductivity, network = NULL, fov = 160, sightline = 100, min_conductivity = 0.4, ..., disp = FALSE)
{
  t0 <- Sys.time()

  if (is(seed, "sf")) seed = sf::st_geometry(seed)
  if (!is(seed, "sfc_LINESTRING")) stop("seed must be sfc_LINESTRING")
  if (length(seed) != 1L) stop("seed must be of length 1")
  if (!is(conductivity, "SpatRaster")) stop("conductivity must be a SpatRaster")

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
  list_coords <- list(start, end)
  list_lines <- list(sf::st_geometry(seed)[[1]])
  list_intersections <- list()
  heading <- get_heading(start, end)

  sub_aoi_conductivity <- query_conductivity_aoi(conductivity, seed, smooth = T)

  #plot(conductivity, col = viridis::inferno(50))
  #plot(terra::ext(sub_aoi_conductivity), col = "red", add = T, alpha = 0.3)
  #plot(seed, add = T, lwd = 3, col = "blue")
  #plot(sub_aoi_conductivity, col = viridis::inferno(50))

  # Mask the existing road to avoid driving a known road
  sub_aoi_conductivity <- mask_existing_network(sub_aoi_conductivity, network)

  # Init view angles as a function of the resolution of the raster and the sightline
  angles_rad <- generate_angles(resolution, sightline, fov)

  # Loop initialization
  current_cost <- 0
  k <- 2
  n <- 1
  overcost = 0
  novercost = 0
  dovercost = 0

  while (dovercost <= 250 & !is.infinite(cost_max))
  {
    if (k %% 5 == 0)
    {
      dist = (k-1)*sightline*0.8
      cat(dist, " m (", get_speed(dist, t0), " km/h)\r", sep = "")
    }
    flush.console()

    # Compute heading from the two previous points
    p1 <- list_coords[[k-1]]
    p2 <- list_coords[[k]]
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
      n = n+1
      #plot(conductivity, col = viridis::inferno(50))
      #plot(raster::extent(sub_aoi_conductivity), col = "red", add = T)
      #plot(seed, add = T, lwd = 3, col = "red")

      sub_aoi_conductivity <- mask_existing_network(sub_aoi_conductivity, network)
      sub_aoi_conductivity <- mask_passage(sub_aoi_conductivity, trace, 0, 5, sf::st_crs(seed))

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
        warning("Drive stopped early. Edge of conductivity raster has been reached.", call. = FALSE)
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

    # The only way to get infinity values is to reach an existing road
    # We need to reduce the sightline to reduce down the approach speed
    if (any(ans$cost >= 9999))
    {
      tmp_angles_rad <- generate_angles(resolution, sightline/2, min(fov, 90))
      ends <- generate_ends(p2, tmp_angles_rad, sightline/2, heading)
      ends$angle <- tmp_angles_rad
      ans <- find_reachable(start, ends, trans, cost_max/2)
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

    list_coords[[k+1]] <- sf::st_geometry(end)[[1]]
    list_lines[[k]] <- L[[1]][-1,]
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
        plot(angles_rad, cost, type = "l", ylim = c(sightline, max(cost) *1.1), col = "green", log = "y")
        lines(angles_rad, ma(cost, n = 7), col = "darkgreen")
        abline(v = angles_rad[c(idx_main, idx_other)], col = "darkgreen",lwd = 3)

        abline(h = cost_max, lty = 3)
        abline(h = 1.5*cost_max, lty = 3)
        abline(v = angles_rad[idx_main] + c(-0.2612, 0.2612), lty = 3)
      }
      terra::plot(terra::crop(sub_aoi_conductivity, terra::ext(aoi) + 80), col = viridis::inferno(50), main = paste0(k,"/", n), range = c(0,1))
      terra::plot(aoi, add = T, col = viridis::inferno(50))
      plot(terra::ext(aoi), add = T)
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
      cat("cost =", current_cost, "\n")
      cat("overcost =", overcost, "\n")
      cat("dovercost =", dovercost, "\n")
    }
  }

  list_lines = list_lines[1:(length(list_lines)-novercost)]

  newline <- do.call(rbind, list_lines) |>
    sf::st_linestring() |>
    sf::st_sfc() |>
    sf::st_set_crs(sf::st_crs(seed))


  if (!is.null(network) && length(network) > 1)
  {
    tail_point = lwgeom::st_endpoint(newline)
    distance_to_network = min(sf::st_distance(network, tail_point))
    if (as.numeric(distance_to_network) < 75)
    {
      u = sf::st_simplify(newline, dTolerance = 2)

      #plot(u, xlim = st_bbox(u)[c(1, 3)] + c(-20,20),  ylim = st_bbox(u)[c(2, 4)] + c(-20,20))
      #plot(network, add = T)
      u = st_extend_line(u, 75, end = "TAIL" )
      #plot(u, xlim = st_bbox(u)[c(1, 3)] + c(-20,20),  ylim = st_bbox(u)[c(2, 4)] + c(-20,20))
      #plot(network, add = T)
      p = sf::st_intersection(u, network)
      if (length(p) == 1 && st_geometry_type(p) == "POINT")
      {
        M = rbind(sf::st_coordinates(newline)[,1:2], sf::st_coordinates(p)[,1:2])
        newline = sf::st_linestring(M) |> sf::st_sfc() |> sf::st_set_crs(sf::st_crs(seed))
      }

      #plot(newline)
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
  scost = as.numeric(na.omit(scost))
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

mask_existing_network <- function(x, network)
{
  if (!is.null(network) && length(network) > 0)
  {
    bb <- sf::st_bbox(x)
    bb <- sf::st_set_crs(bb, sf::st_crs(network))
    mask <- sf::st_crop(network, bb)
    if (length(mask) > 0)
    {
      mask <- sf::st_buffer(mask, dist = 10)
      x <- terra::mask(x, terra::vect(mask), inverse = TRUE, updatevalue = 0)
    }
  }

  return(x)
}

mask_passage <- function(raster, lines, start_cut, end_cut, crs)
{
  len =as.numeric(st_length(lines))
  from = 1 - (len-start_cut)/len
  to = (len-end_cut)/len
  if (length(from) > 1)
  {
    from = mean(from)
    to = mean(to)
  }
  if (is(lines, "sfg")) lines <- sf::st_sfc(lines)
  mask <- lwgeom::st_linesubstring(lines, from, to)
  mask <- sf::st_buffer(mask, dist = 5, endCapStyle = "FLAT")
  mask <- sf::st_set_crs(mask, crs)
  raster <- terra::mask(raster, terra::vect(mask), inverse = TRUE, updatevalue = 0)
  return(raster)
}

get_speed = function(dist, t0)
{
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
