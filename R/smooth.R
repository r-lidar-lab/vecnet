smooth_network = function(net)
{
    smoothed = lapply(net, adjust_spline)
    smoothed = do.call(c, smoothed)
    sf::st_crs(smoothed) = sf::st_crs(net)
    return(smoothed)
}

adjust_spline = function(points)
{
  # Adjust a spline to create a smooth line from points
  xroad <- sf::st_coordinates(points)[,1]
  yroad <- sf::st_coordinates(points)[,2]

  # Retrieve the "temporal" position along the line so we can order the points later
  xend <- xroad[length(xroad)]
  yend <- yroad[length(yroad)]
  xstart <- xroad[1]
  ystart <- yroad[1]
  dx  <- c(xroad, xend) - c(xstart, xroad)
  dy  <- c(yroad, yend) - c(ystart, yroad)
  dd  <- sqrt(dx*dx+dy*dy)
  dist <- cumsum(dd[-length(dd)])
  troad <- dist/dist[length(dist)]

  # Weigth the first and last node
  wroad <- rep(1, length(xroad))
  wroad[0] = 1e10
  wroad[length(wroad)] = 1e10

  ux <- stats::smooth.spline(troad, xroad, wroad, spar = 0.4, all.knots = TRUE)
  uy <- stats::smooth.spline(troad, yroad, wroad, spar = 0.4, all.knots = TRUE)

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
    wroad = wroad[rm]
  }

  # Maybe it is possible to have an extra point. This case is handled here but never observed
  if (length(ux$x) > length(troad))
  {
    print(dput(points))
    stop("Different length")
  }

  # The adjusted spline can be far from some outlier points. This is the role of the spline
  # to smooth the road. But if some important outliers are detected far from the spline (more than
  # 10 m) we remove those points and fit again the spline
  rm2 <- FALSE
  rmx <- abs(ux$y - xroad) < 10
  rmy <- abs(uy$y - yroad) < 10

  if (any(!rmx) | any(!rmy))
  {
    rm2   <- rmx & rmy
    xroad <- xroad[rm2]
    yroad <- yroad[rm2]
    troad <- troad[rm2]
    wroad <- wroad[rm2]

    if (length(xroad) > 3)
    {
      ux <- stats::smooth.spline(troad, xroad, w = wroad)
      uy <- stats::smooth.spline(troad, yroad, w = wroad)
    }
    else
    {
      rm2 <- FALSE
    }
  }

  if (!isFALSE(rm))  spline <- spline[rm,]
  if (!isFALSE(rm2)) spline <- spline[rm2,]

  spline <- cbind(ux$y, uy$y)
  spline[1,1] = xstart
  spline[1,2] = ystart
  spline[nrow(spline),1] = xend
  spline[nrow(spline),2] = yend
  spline <- sf::st_sfc(sf::st_linestring(spline))
  sf::st_crs(spline) <- st_crs(points)
  spline = sf::st_simplify(spline, dTolerance = 2)
  return(spline)
}
