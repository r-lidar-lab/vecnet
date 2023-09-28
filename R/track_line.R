#' Drive and vectorize an unknown road
#'
#' Drive along an unknown road in a probability map starting from small piece
#' of road segment pointing in the right direction and vectorize one road by returning a spatial
#' line. This function is used by \link{vectorize_network} and is not intended to be used
#' by regular users except for trial-and-error tests.
#'
#' @param seed \code{sfc_LINESTRING}. A seed to start driving a road
#' @param map \code{SpatRaster}. A raster with value ranging in [0,1] where roads pixels are
#' close to 1 and non road pixels are close to 0 (see references). Also called conductivity map.
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
#' @param ... Unused
#' @param disp bool. Display in realtime the progress on images. For debugging purposes.
#'
#' @return list.
#' \enumerate{
#' \item{\code{road} contains the road found in a vectorial format (\code{sfc}) }
#' \item{\code{seed} contains all the potential intersection and seeds found along
#'  the road in a vectorial format (\code{sfc}) }}
#'
#' @references
#' Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
#' (2023) Vectorial and topologically valid segmentation of forestry road networks from ALS data.
#' Journal of Applied Earth Observation and Geoinformation.\cr\cr
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
#' seeds <- system.file("extdata", "seeds2.gpkg", package = "vecnet")
#' map <- rast(map)
#' seeds <- st_read(seeds, quiet = TRUE)
#' seed <- st_geometry(seeds)[1]
#'
#' res <- track_line(seed, map)
#'
#' plot(map, col = viridis::inferno(25), smooth =T)
#' plot(seed, add = TRUE, col = "green", lwd = 3)
#' plot(res$road, add = TRUE, col = "red", lwd = 2)
#' plot(res$seeds, add = TRUE, col = "green", lwd = 3)
#' }
track_line <- function(seed,
                       map,
                       network = NULL,
                       fov = 160,
                       sightline = 100,
                       min_conductivity = 0.6,
                       th_conductivity = 0.1,
                       ...,
                       disp = FALSE)
{
  # network = NULL ; fov = 160 ; sightline = 100 ; min_conductivity = 0.6 ; th_conductivity = 0.1 ; partial_gap_size = 2.5 ;
  debug = FALSE

  # Defensive programming
  if (methods::is(seed, "sf")) seed = sf::st_geometry(seed)
  if (!methods::is(seed, "sfc_LINESTRING")) stop("seed must be sfc_LINESTRING")
  if (length(seed) != 1L) stop("seed must be of length 1")
  if (!methods::is(map, "SpatRaster")) stop("map must be a SpatRaster")
  if (!is.null(network)) network = sf::st_geometry(network)
  resolution <- terra::res(map)
  if (resolution[1] != resolution[2]) stop("'map' raster must have the same resolution in both X and Y axis.", call. = FALSE)
  resolution <- resolution[1]

  t0 <- Sys.time()

  # We are working with large rasters not loaded in memory. We must load in memory some small parts
  # iteratively wile driving. This object handle the logic of querying subsets in a larger map
  mapper <- MapManager$new(map, th_conductivity = th_conductivity, smooth = TRUE)
  mapper$set_search_seed_mode(list(...)$seed_mode)
  mapper$query_aoi(seed)
  mapper$mask_network(network)

  # This object handles the logic of least cost based driving
  driver = DrivingAgent$new(sightline = sightline, fov = fov, resolution = resolution, min_conductivity = min_conductivity)
  driver$set_seed(seed)

  # Loop while the driver can move forward
  while (!driver$reached_end())
  {
    driver$print_speed(t0)
    driver$query_sightview(mapper$aoi)

    if (driver$is_on_edge())
    {
      # If the driving agent is on the edge of the loaded aoi we must load a new aoi ahead
      mapper$query_aoi(driver$get_seed())
      mapper$mask_network(network)
      mapper$mask_trace(driver$get_trace())
      driver$query_sightview(mapper$aoi)

      # If the driving agent is still on the edge of the loaded aoi it means that we are really on
      # the edge of the full raster. Divide the sightline to reach the edge doing smaller steps
      if (driver$is_on_edge())
        driver$divide_sightline()

      if (driver$is_on_edge())
        driver$divide_sightline()

      if (driver$is_on_edge())
        break
    }

    # We are not on an edge we can compute the costs
    tryCatch({
      driver$compute_cost()
    },
    error = function(e) { warning(paste0("Internal errror while computing the cost: stop driving the road. (", e, ")")); })

    # If there is no path it means no point was reachable
    if (!driver$has_path())
    {
      break
    }

    if (debug) driver$plot()

    # We reached another road from the existing network
    if (driver$has_reached_a_road())
    {
      driver$network_connection(network)
      break
    }

    # Leave a trace of already seen roads and protect against infinite loops /!\
    mapper$mask_trace(driver$get_seed())
    if (driver$has_new_intersections())
      mapper$mask_trace(driver$get_intersections())

    if (disp)
    {
      bb = sf::st_buffer(driver$get_start(), 1.25*sightline, endCapStyle = "SQUARE") |> terra::vect()
      plot(terra::crop(mapper$aoi, terra::ext(bb)), col = viridis::inferno(25))
      driver$plot(add = T)
    }

    driver$move()
  }

  driver$print_speed(t0, done = TRUE)

  return(driver$get_road())
}

# ==== MapMananer =====

# We are working with large rasters not loaded in memory. We must load in memory some small part
# iteratively wile driving. This object handle the logic of querying subsets in a larger map

MapManager <- R6::R6Class("MapManager",
cloneable = FALSE,
public = list(
  map = NULL,
  aoi = NULL,
  trace = NULL,
  initialize = function(map, th_conductivity = 0.1, smooth = FALSE)
  {
    self$map = map
    private$smooth = smooth
    private$th_conductivity = th_conductivity
    private$mask_value = 0
    private$buffer = 10
    self$set_search_seed_mode(FALSE)
  },

  query_aoi = function(seed)
  {
    resolution <- terra::res(self$map)[1]
    ahead = 900
    behind = 100
    side = 400
    smooth = TRUE

    seed <- sf::st_simplify(seed, FALSE, dTolerance = 5)
    seed <- sf::st_sfc(seed)
    sf::st_crs(seed) <- sf::st_crs(self$map)
    private$seed <- seed

    line <- seed |> st_extend_line(c(ahead, behind)) |>  sf::st_buffer(side, endCapStyle = "FLAT")
    self$aoi <- terra::crop(self$map, terra::vect(line))

    if (private$smooth)
    {
      w = matrix(1,3,3)
      w[5] = 10
      w = w/sum(w)
      self$aoi  = terra::focal(self$aoi, w = w, fun = "mean")
    }

    self$aoi[is.na(self$aoi)] <- 0
    self$aoi[self$aoi < private$th_conductivity] = private$th_conductivity
  },

  mask_network = function(network)
  {
    self$aoi <- mask_existing_network(self$aoi, network, private$mask_value, private$buffer)

    if (private$seed_mode)
    {
      bb = sf::st_as_sfc(sf::st_bbox(self$aoi))
      polygon = sf::st_cast(network, "POLYGON")
      polygon = sf::st_difference(bb, polygon)
      self$aoi = terra::mask(self$aoi, terra::vect(polygon), inverse = TRUE, updatevalue = private$th_conductivity)
    }
  },

  mask_trace = function(trace)
  {
    #start_cut = 0
    #end_cut = 5

    #len = as.numeric(sf::st_length(trace))
    #from = 1 - (len-start_cut)/len
    #to = (len-end_cut)/len
    #if (length(from) > 1)
    #{
    #  from = mean(from)
    #  to = mean(to)
    #}

    if (methods::is(trace, "sfg")) trace <- sf::st_sfc(trace)
    #mask <- lwgeom::st_linesubstring(trace, from, to)
    mask <- sf::st_buffer(trace, dist = private$buffer, endCapStyle = "FLAT")
    mask <- sf::st_set_crs(mask,  sf::st_crs(private$seed))
    self$aoi <- terra::mask(self$aoi, terra::vect(mask), inverse = TRUE, updatevalue = 0.001)
  },

  set_search_seed_mode = function(val)
  {
    if (isTRUE(val))
    {
      private$mask_value = 2
      private$buffer = 2
      private$seed_mode = TRUE
      private$smooth = FALSE
    }
    else
    {
      private$mask_value = 0
      private$buffer = 10
      private$seed_mode = FALSE
    }
  },

  plot = function(aoi = FALSE)
  {
    if (isFALSE(aoi))
    {
      terra::plot(self$map, col = viridis::inferno(50))
      terra::plot(terra::ext(self$aoi), col = "red", add = T, alpha = 0.3)
      plot(private$seed, add = T, lwd = 3, col = "blue")
    }
    else
    {
      terra::plot(self$aoi, col = viridis::inferno(50))
    }
  }
),
private = list(
  smooth = NULL,
  th_conductivity = NULL,
  mask_value = NULL,
  seed = NULL,
  buffer = NULL,
  seed_mode = NULL)
)

# ==== Driving Agent ====

DrivingAgent <- R6::R6Class("DrivingAgent",
cloneable = FALSE,
private = list(
  sightline = NULL,          # Sightline of the driving agent input by the user
  current_sightline = NULL,  # The sightline maybe reduced locally. We use current_sightline and sightline hold the original input
  resolution = NULL,         # Resolution of the map on which we are suppose to drive
  rcost = NULL,              # Relative cost compared to the maximal allowed value (which is not necessarily cost_max)
  cost_max = NULL,           # Maximal cost allowed (user input). This value may be locally changed.
  fov = NULL,                # Field of view (in degrees) of the driving agent
  angle = NULL,              # The angle of the 0 degree for the field of view
  map = NULL,                # [SpatRaster] on which the driving agent drives on step
  trans = NULL,              # [TransitionMatrix] the transition matrix associated with map
  start = NULL,              # [sfc] the position of the driving agent
  ends = NULL,               # [sfc] the potentially reachable points
  p1 = NULL,                 # [sfc] The first point of the seed used to compute the angle
  p2 = NULL,                 # [sfc] The second point of the seed used to compute the angle (same than start actually)
  k = NULL,                  # A counter of move
  cost_profile = NULL,       # [vect] The angle-cost profile
  cost_peaks = NULL,         # [data.frame] The location of the peaks in the cost profile
  lcp = NULL,                # [sfc] line of the least cost path
  intersections = NULL,      # [sfc] lines of the intersection seeds founds
  novercost = NULL,          # Number of time we moved forward despite a cost higher than the maximum allowed cost
  dovercost = NULL,          # Distance travelled  despite a cost higher than the maximum allowed cost
  list_lines = NULL,         # [list<sfg>] list of sfg linestring of each step driven
  list_intersections = NULL, # [list<sfg>] list of sfg linestring of each intersection seed
  crs = sf::NA_crs_,         # [crs]
  sightline_min = 15,        # Miminum allowed sightline
  pahead = 0.8,              # When moving forward we do not move of 100% of the sightline but rather 80% (this prevent missing intersection)
  partial_gap_size = 2.5
),

public = list(
  initialize = function(sightline = 100, fov = 160, start = c(0,0), resolution = 2, min_conductivity = 0.6, crs = sf::NA_crs_)
  {
    private$sightline <- sightline
    private$fov <- fov
    private$start <- start
    private$angle <- 0
    private$resolution <- resolution
    private$cost_max <- sightline * 1/min_conductivity
    private$current_sightline  <- sightline
    private$k <- 2L
    private$novercost = 0L
    private$dovercost = 0
    private$list_lines <- list()
    private$list_intersections <- list()
  },

  query_sightview = function(val)
  {
    bb <- sf::st_bbox(c(private$start, sf::st_geometry(private$ends)))
    bb <- terra::ext(bb) + 10
    aoi <- terra::crop(val, bb)

    p <- sf::st_polygon(list(rbind(sf::st_coordinates(private$start), sf::st_coordinates(private$ends), sf::st_coordinates(private$start))))
    p <- sf::st_sfc(p)
    p <- sf::st_buffer(p, 2)
    p <- terra::vect(p)
    terra::crs(p) <- terra::crs(val)

    private$map <- terra::mask(aoi, p)
    private$trans <- transition(private$map, geocorrection = TRUE)
  },

  get_novercost = function()
  {
    return(private$novercost)
  },

  set_start = function(val)
  {
    private$start <- val

    angles_rad = generate_angles(private$resolution, private$current_sightline,  private$fov)
    private$ends <- generate_ends(private$p2, angles_rad, private$current_sightline, private$angle)
    private$ends$angle <- angles_rad
  },

  get_start = function()
  {
    return(private$start)
  },

  set_seed = function(line)
  {
    if (is.na(private$crs))
      private$crs <- sf::st_crs(line)

    if (length(private$list_lines) == 0)
      private$list_lines[[1]] = sf::st_geometry(line)[[1]]

    l  <- lwgeom::st_linesubstring(line, 0.5, 1)
    p1 <- lwgeom::st_startpoint(l)
    p2 <- lwgeom::st_endpoint(l)

    private$p1 = p1
    private$p2 = p2

    from <- sf::st_coordinates(p1)
    to   <- sf::st_coordinates(p2)
    M <- rbind(from, to)
    Ax <- M[2,1]
    Ay <- M[2,2]
    Bx <- M[1,1]
    By <- M[1,2]

    private$angle <- atan2(Ay-By, Ax-Bx)
    self$set_start(p2)
  },

  get_seed = function()
  {
    n = length(private$list_lines)
    return(private$list_lines[[n]])
  },

  set_fov = function(val)
  {
    private$fov <- val
  },

  set_sightline = function(val)
  {
    private$current_sightline <- val
    private$set_sighline <- val
  },

  get_intersections = function()
  {
    return(private$intersections)
  },

  get_road = function()
  {
    roads <- private$list_lines
    roads <- roads[1:(length(roads) - private$novercost)]
    first = private$list_lines[[1]]
    roads <- lapply(roads, function(x) { x <- x[] ; return(x[-1,])})
    roads[[1]] = first

    newline <- do.call(rbind, roads) |>
      sf::st_linestring() |>
      sf::st_sfc() |>
      sf::st_set_crs(private$crs)

    intersections <- NULL
    if (length(private$list_intersections) >= 1)
    {
      intersections <- do.call(c, private$list_intersections)
      intersections <- sf::st_set_crs(intersections, private$crs)

      p = lwgeom::st_startpoint(intersections)
      d = as.numeric(sf::st_distance(p, newline))
      intersections <- intersections[d < 1]
    }


    return(list(road = newline, seeds = intersections))
  },

  get_trace = function()
  {
    u = self$get_road()
    road = u$road
    seed = u$seeds
    return(c(road, seed))
  },

  has_new_intersections = function()
  {
    return(!is.null(private$intersections))
  },

  has_path = function()
  {
    return(!is.null(self$least_cost_path()))
  },

  network_connection = function(network)
  {
    # Create a prolongation line
    angle <- private$angle
    xstart <- sf::st_coordinates(private$start)[1]
    ystart <- sf::st_coordinates(private$start)[2]
    xend <- xstart+50*cos(angle)
    yend <- ystart+50*sin(angle)
    m = matrix(c(xstart, xend, ystart, yend), ncol = 2)
    prolongation <- sf::st_sfc(sf::st_linestring(m), crs = sf::st_crs(network))

    # Intersection between the prolongation and the existing network
    p = sf::st_intersection(prolongation, network)

    # No intersection: double the length of the prolongation
    if (length(p) == 0)
    {
      xend <- xstart+100*cos(angle)
      yend <- ystart+100*sin(angle)
      m = matrix(c(xstart, xend, ystart, yend), ncol = 2)
      prolongation <- sf::st_sfc(sf::st_linestring(m), crs = sf::st_crs(network))

      # Intersection between the prolongation and the existing network
      p = sf::st_intersection(prolongation, network)
    }

    if (length(p) > 1)
    {
      p = p[1]
      warning("Connection problem. More than one line intersection")
    }

    # Still no intersection we may be in the wrong direction e.g. if we are trying to
    # connect to a parallel line and we have 0 tolerance. Connect to the NN
    if (length(p) == 0)
    {
      bbb = sf::st_buffer(private$start, dist = 50)
      sf::st_crs(private$start) = sf::st_crs(network)
      sf::st_crs(bbb) = sf::st_crs(network)
      subnet = sf::st_intersection(network, bbb)
      prolongation = sf::st_nearest_points(private$start, subnet)
    }
    # we have a connection but we must round the value to match a node of the line
    else
    {
      bbb = sf::st_buffer(p, dist = 5)
      sf::st_crs(p) = sf::st_crs(network)
      sf::st_crs(bbb) = sf::st_crs(network)
      subnet = sf::st_intersection(network, bbb)
      subnet = sf::st_cast(subnet, "POINT")
      i = sf::st_nearest_feature(p, subnet)
      p = subnet[i]
      prolongation = sf::st_coordinates(prolongation)[,1:2]
      prolongation[2,] = sf::st_coordinates(p)
      prolongation = sf::st_sfc(sf::st_linestring(prolongation), crs = sf::st_crs(network))
    }

    private$novercost = 0L
    private$dovercost = 0

    private$list_lines[[private$k]] = prolongation[[1]]
  },

  least_cost_path = function()
  {
    return(private$lcp)
  },

  has_reached_a_road = function()
  {
    any_infinite_cost(private$cost_profile)
  },

  divide_sightline = function(div = 2)
  {
    private$current_sightline <- private$current_sightline/2
    angles_rad <- generate_angles(private$resolution, private$current_sightline, private$fov)
    private$ends <- generate_ends(private$p2, angles_rad, private$current_sightline, private$angle)
    private$ends$angle <- angles_rad
    self$query_sightview(private$map)
  },

  is_on_edge = function()
  {
    val <- terra::extract(private$map, sf::st_coordinates(private$ends))[[1]]
    return(any(is.na(val)))
  },

  move = function()
  {
    private$k = private$k+1
    private$current_sightline = private$sightline
    private$cost_profile = NULL
    private$cost_peaks = NULL
    private$intersections = NULL

    self$set_seed(private$lcp)

    private$lcp = NULL
  },

  reached_end = function()
  {
    if (is.infinite(private$cost_max))
      return(TRUE)

    # @param partial_gap_size numeric. The algorithm can drive on low conductivity segments even if the
    # conductivity is too low and should trigger a stop signal if it looks like a sharp line and an
    # high conductivity road is reach after the gap (see references). The allowed distance is a multiple
    # of the sightline. Default is 2.5 which roughly means that it can follow a too low conductivity
    # track for two iterations and at the third one it must be back on a high conductivity track.
    if (private$dovercost > private$partial_gap_size*private$sightline)
      return(TRUE)

    return(FALSE)
  },

  compute_cost = function()
  {
    ans <- find_reachable(private$start, private$ends, private$trans, private$cost_max)

    if (is.null(ans)) return(NULL)

    # The only way to get infinite values is to reach an existing road
    # We need to reduce the sightline to reduce down the approach speed and connect accurately
    current_cost = private$cost_max
    while (any_infinite_cost(ans$cost) & private$current_sightline > private$sightline_min)
    {
      self$divide_sightline()
      current_cost = current_cost/2
      ans <- find_reachable(private$start, private$ends, private$trans, current_cost)
    }

    private$cost_profile = ans$cost
    private$cost_peaks = ans$minima

    # We now have the cost to reach each point of the sightline and the peaks in the cost
    # profile indicating which point we must reach
    cost <- private$cost_profile

    # The main direction
    idx_main <- private$cost_peaks$idx[1]
    depth_main = private$cost_peaks$depth[1]
    end <- private$ends[idx_main,]

    # Compute the shortest path between the starting points and the valid end points
    L <- gdistance::shortestPath(private$trans, sf::st_coordinates(private$start), sf::st_coordinates(end), output = "SpatialLines")
    L <- sf::st_geometry(sf::st_as_sf(L))
    C <- sf::st_coordinates(L)
    C <- C[1:(nrow(C)*private$pahead),-3]
    L <- sf::st_sfc(sf::st_linestring(C))
    private$lcp <- L
    private$list_lines[[private$k]] <- L[[1]]

    # The other direction
    idx_other <- private$cost_peaks$idx[-1]
    private$rcost <- private$cost_peaks$rcost[1]
    ends_other <- private$ends[idx_other,]
    cost_other <- private$cost_peaks$rcost[-1]

    # Check the angle between the main road and the other road. If the angle is < 15 degrees there
    # are false positive.
    # (This is probably not useful anymore but was need to fix bugs when using very poor map)
    delta_angle = abs(end$angle - ends_other$angle)*180/pi
    ends_other  = ends_other[delta_angle > 15,]
    cost_other  = cost_other[delta_angle > 15]

    # Now we are processing the intersection founds (if any)
    I <- NULL
    if (length(sf::st_geometry(ends_other)) > 0)
    {
      I <- gdistance::shortestPath(private$trans, sf::st_coordinates(private$start), sf::st_coordinates(ends_other), output = "SpatialLines") |> suppressWarnings()
      I <- sf::st_geometry(sf::st_as_sf(I))
      I <- sf::st_difference(I,L)

      I = Filter(function(x){
        if (sf::st_geometry_type(x) != "LINESTRING")
          return(FALSE)

        if (lwgeom::st_endpoint(L) ==  lwgeom::st_startpoint(x))
          return(FALSE)

        return(TRUE)
      }, I)

      private$intersections <- I
      private$list_intersections <- c(private$list_intersections, list(I))
    }

    # If the current cost is > 1 it means that we drove on a track that is partially or totally
    # not labelled as road. We record this state
    if (private$rcost > 1)
    {
      private$novercost = private$novercost + 1L
      private$dovercost = private$dovercost + private$current_sightline
    }
    else
    {
      private$novercost = 0L
      private$dovercost = 0
    }
  },

  plot = function(add = FALSE)
  {
    # Plot the background map
    if (!is.null(private$map))
    {
      terra::plot(private$map, col = viridis::inferno(50), range = c(0,1), add = add)
      add = TRUE
    }

    # If we already computed the lcp, then recompute all the possible paths (they were not recorded)
    # (For display only)
    if (!is.null(private$lcp))
    {
      tryCatch(
      {
        paths = gdistance::shortestPath(private$trans, sf::st_coordinates(private$start), sf::st_coordinates(private$ends), output = "SpatialLines") |> suppressWarnings()
        paths = sf::st_as_sf(paths)
        paths$cost = private$cost_profile
        plot(paths["cost"], add = add, pal = viridis::viridis, breaks = "quantile")
      },
      error = function(x) print(x))

      end <- lwgeom::st_endpoint(private$lcp)
      plot(private$lcp, add = T, col = "red", lwd = 4)
      plot(end, add = add, col = "green", pch = 19, cex = 2)

      graphics::title(paste0("rcost: ", round(private$rcost,2), " | radius: ", private$current_sightline, "\n k: ", private$k, " | nover: ", private$novercost, " (+", private$dovercost, "m)"))
    }

    # Plot the intersection
    if (!is.null(private$intersections))
    {
      plot(private$intersections, add = add, col = "blue", lwd = 3)
      add = TRUE
    }


    if (!is.null(private$ends))
    {
      if (is.null(private$cost_profile))
      {
        plot(private$ends$geometry, add = add, col = "red", pch = 19, cex = 0.3)
      }
      else
      {
        geom = private$ends
        geom$cost = private$cost_profile
        plot(geom["cost"], add = add, pal = viridis::viridis, pch = 19, cex = 0.25, breaks = "quantile")
      }

      add = TRUE
    }

    if (!is.null(private$start))
    {
      plot(private$start, add = add, col = "red", pch = 19)
      add = TRUE
    }

    if (!is.null(private$p1))
    {
      plot(private$p1, add = add, col = "blue", pch = 19)
      add = TRUE
    }
  },

  print_speed = function(t0, done = FALSE)
  {
    if (isFALSE(done))
    {
      if (private$k %% 5 == 0)
      {
        dist = (private$k-1)*private$sightline*private$pahead
        cat(dist, " m (", get_speed(dist, t0), " km/h)\r", sep = "")
        utils::flush.console()
      }
    }
    else
    {
      h <- km <- NULL
      tf <- Sys.time()
      dt <- tf-t0
      dt <- round(units::as_units(dt),1)
      dist <- sf::st_length(self$get_road()$road)
      hdt = dt
      units(hdt) <- units::make_units(h)
      units(dist) <- units::make_units(km)
      dist = round(dist, 2)
      speed = round(dist/hdt)
      cat("Processed ended in", dt, units::deparse_unit(dt),
          ": road of", dist, units::deparse_unit(dist),
          "driven at", speed, units::deparse_unit(speed),
          "\n")
    }
  })
)

# ==== Tools ====

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
  crs <- sf::st_crs(center)
  center <- sf::st_coordinates(center)
  X <- center[1] + sightline * cos(direction + angles)
  Y <- center[2] + sightline * sin(direction + angles)
  M <- data.frame(X,Y)
  M <- sf::st_as_sf(M, coords = c("X", "Y"))
  sf::st_crs(M) = crs
  return(M)
}

find_reachable <- function(start, ends, trans, cost_max)
{
  if (length(sf::st_geometry(ends)) < 6) return(NULL)

  cost <- gdistance::costDistance(trans, sf::st_coordinates(start), sf::st_coordinates(ends))
  cost <- as.numeric(cost)

  if (all(is.infinite(cost))) return(NULL)

  cost[is.infinite(cost)] <- max(99999, max(cost[!is.infinite(cost)]))

  # Select local minima of cost
  if (length(cost) > 36)
    smooth <- 9
  else if (length(cost > 18))
    smooth = 2
  else
    smooth <- 1

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

has_na_in_sightline = function(map, ends)
{
  val <- terra::extract(map, sf::st_coordinates(ends))[[1]]
  return(any(is.na(val)))
}

any_infinite_cost = function(cost)
{
  return(any(cost >= 99999))
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



# if (length(angles_rad) == length(cost))
# {
#
#   #svglite::svglite("/home/jr/Documents/Ulaval/2022 Post-doc/Articles/vectnet/img/svg/followroadcostprofile.svg", height = 8, width = 8)
#   graphics::par(family = "serif")
#   plot(angles_rad*180/pi, cost, type = "l", lwd = 2, ylim = c(sightline, max(cost) *1.1), col = "red", xlab = "Angle (\u00B0)", ylab = "Cost (1)")
#   smooth <- 9
#   scost = ma(cost, n = smooth)
#   scost = ma(scost, n = 5)
#   graphics::lines(angles_rad*180/pi, scost, col = "darkgreen", lwd = 2)
#   graphics::abline(v = angles_rad[c(idx_main, idx_other)]*180/pi, col = "black",lwd = 2)
#   graphics::abline(h = cost_max, lty = 3)
#   graphics::abline(h = 2.5*cost_max, lty = 3)
#   graphics::legend(80, 990, legend=c("Cost profile", "Smoothed profile"),
#          col=c("red", "darkgreen"), lty=1:2, cex=0.8)
#   #abline(v = angles_rad[idx_main]*180/pi + c(-0.2612, 0.2612)*180/pi, lty = 3)
#   #dev.off()
# }
