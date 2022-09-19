#' @importClassesFrom  gdistance TransitionLayer
transition <- function(conductivity, directions = 8, geocorrection = TRUE)
{
  x = conductivity
  use_terra = FALSE

  if (methods::is(x, "RasterLayer"))
  {
    bb = raster::extent(x)
    ncells = raster::ncell(x)
    val = raster::values(x)
  }
  else if (methods::is(x, "SpatRaster"))
  {
    bb = terra::ext(x)
    bb = raster::extent(bb[1:4])
    ncells = terra::ncell(x)
    val = terra::values(x)[,1]
    use_terra = TRUE
  }
  else
  {
    stop("Only RasterLayer and SpatRaster are supported")
  }

  symm = TRUE

  tr <- methods::new("TransitionLayer",
            nrows=as.integer(nrow(x)),
            ncols=as.integer(ncol(x)),
            extent=bb,
            crs=sp::CRS(),
            transitionMatrix = Matrix::Matrix(0, ncells,ncells),
            transitionCells = 1:ncells)

  transitionMatr <- gdistance::transitionMatrix(tr)
  Cells <- which(!is.na(val))

  if (use_terra)
    adj <- terra::adjacent(x, cells=Cells, pairs=TRUE, directions=directions)
  else
    adj = raster::adjacent(x, cells=Cells, pairs=TRUE, directions=directions)

  if(symm)
    adj <- adj[adj[,1] < adj[,2],]

  dataVals <- cbind(val[adj[,1]], val[adj[,2]])
  transition.values <- rowMeans(dataVals)
  transition.values[is.na(transition.values)] <- 0

  if(!all(transition.values >= 0))
    warning("transition function gives negative values")

  transitionMatr[adj] <- as.vector(transition.values)
  if(symm)
    transitionMatr <- Matrix::forceSymmetric(transitionMatr)

  gdistance::transitionMatrix(tr) <- transitionMatr
  gdistance::matrixValues(tr) <- "conductance"


  #trans <- gdistance::transition(conductivity, transitionFunction = mean, directions = 8)
  trans = tr

  if (isTRUE(geocorrection))
  {
    trans <- gdistance::geoCorrection(trans)
  }

  return(trans)
}
