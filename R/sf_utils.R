#' Get heading of both ends of a line
#'
#' Retrieve heading of both ends of a line, pointing away from it.
#'
#' @param line  line (\code{sfc} or \code{sfg} format)
#'
#' @return numeric of length 2 expressing respectively the head and tail heading in degrees
#' (range [-180°, 180°] using \code{atan2()}.
#' @noRd
st_ends_heading <- function(line)
{
  M <- sf::st_coordinates(line)
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)

  headings <- mapply(i, j, FUN = function(i, j) {
    Ax = M[i-j,1]
    Ay = M[i-j,2]
    Bx = M[i,1]
    By = M[i,2]
    atan2(Ay-By, Ax-Bx)*180/pi
  })
  names(headings) <- c("head", "tail")
  return(headings)
}


st_join_linestring = function(x,y)
{
  X = sf::st_coordinates(x)
  Y = sf::st_coordinates(y)
  M = rbind(X,Y)[,1:2]
  L = sf::st_sfc(sf::st_linestring(M))
  return(L)
}

#' Extend line by given distance
#'
#' Extend one or both ends of a line by a given distance.
#' No new vertices are added, instead, the first/last is moved
#' to its new location following the same heading as before.
#'
#' @param line  line (\code{sfc} or \code{sfg} format)
#' @param distance  numeric; distance by which the line will be extended. In case of \code{end = "BOTH"},
#' a vector of length 2 can be provided to extend respectively the head and tail end by different values.
#' @param end  character; (\code{BOTH}, \code{HEAD} or \code{TAIL}; define which end will be extended.
#'
#' @return Same line as input but now extended.
#' @noRd
st_extend_line <- function(line, distance, end = "BOTH")
{
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")

  M <- sf::st_coordinates(line)[,-3]
  keep <- !(end == c("TAIL", "HEAD"))

  ends <- c(1, nrow(M))[keep]
  headings <- st_ends_heading(line)[keep] / 180 * pi
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])

  M[ends, 1:2] <- M[ends, 1:2] + distances[keep] * c(cos(headings), sin(headings))
  newline <- sf::st_linestring(M)

  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- sf::st_sfc(newline, crs = sf::st_crs(line))

  return(newline)
}

#' @examples
#' f <- system.file("extdata", "sinuosity", package="vectnet")
#' m = sf::st_read(f)
#' m = sf::st_transform(m, 27572)
sinuosity <- function(x)
{
  lines <- x
  n <- length(lines)
  S <- numeric(n)
  for (i in 1:n)
  {
    line <- lines[i]

    if (sf::st_geometry_type(line) != "LINESTRING")
      stop(paste0("Geometry ", i, " is not a LINESTRING"))

    spline <- adjust_spline(line)
    S[i] = sf::st_length(line)/sf::st_length(spline)
  }

  round(S,2)
}

adjust_spline = function(points)
{
  # Adjust a spline to create a smooth line from points
  xroad <- sf::st_coordinates(points)[,1]
  yroad <- sf::st_coordinates(points)[,2]

  if (length(xroad) <= 4)
    return(points)

  n = length(xroad)
  troad <- cumsum(c(0,sqrt((xroad[-1]-xroad[-n])^2 + (yroad[-1]-yroad[-n])^2)))

  ux <- stats::smooth.spline(troad, xroad, spar = 0.5, all.knots = TRUE)
  uy <- stats::smooth.spline(troad, yroad, spar = 0.5, all.knots = TRUE)

  # Sometime one point is missing in the spline for an unknown reason. If the output
  # does not have the same length than the input we resize the input to fit the output
  # and prevent failures
  rm <- FALSE
  if (length(ux$x) < length(troad))
  {
    rm <- troad %in% ux$x
    xroad = xroad[rm]
    yroad = yroad[rm]
    troad = troad[rm]
  }

  # Maybe it is possible to have an extra point. This case is handled here but never observed
  if (length(ux$x) > length(troad))
  {
    print(dput(points))
    stop("Different length")
  }

  spline = cbind(ux$y, uy$y)
  spline = sf::st_sfc(sf::st_linestring(spline))

  return(spline)
}
