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
    adj <- raster::adjacent(x, cells=Cells, pairs=TRUE, directions=directions)

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
    trans <- geocorrection(trans)

  return(trans)
}

# Modification, simplification and use terra of the function gdistance::geoCorrection()
geocorrection = function(x)
{
  adj <- gdistance::adjacencyFromTransition(x)
  cell1 = terra::xyFromCell(x, adj[,1])
  cell2 = terra::xyFromCell(x,adj[,2])

  distance = (cell1-cell2)^2
  distance = sqrt(distance[,1] + distance[,2])

  if (gdistance::matrixValues(x) == "conductance")
    correctionValues <- 1/distance
  else if (gdistance::matrixValues(x) == "resistance")
    correctionValues <- distance
  else
    stop("Internal error")

  i <- as.integer(adj[,1] - 1)
  j <- as.integer(adj[,2] - 1)
  xv <- as.vector(correctionValues) #check for Inf values!
  dims <- terra::ncell(x)
  correctionMatrix <- methods::new("dgTMatrix", i = i, j = j, x = xv, Dim = as.integer(c(dims,dims)))
  correctionMatrix <- (methods::as(correctionMatrix,"sparseMatrix"))

  if (methods::is(gdistance::transitionMatrix(x), "dsCMatrix")) #isSymmetric?
    correctionMatrix <- Matrix::forceSymmetric(correctionMatrix)

  transitionCorrected <- correctionMatrix * gdistance::transitionMatrix(x)
  gdistance::transitionMatrix(x) <- transitionCorrected
  return(x)
}

