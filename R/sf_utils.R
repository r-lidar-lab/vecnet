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
