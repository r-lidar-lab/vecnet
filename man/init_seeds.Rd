% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_seeds.R
\name{init_seeds}
\alias{init_seeds}
\title{Find seeds for the vectorization}
\usage{
init_seeds(map, contour = sf::st_as_sfc(sf::st_bbox(map)), ...)
}
\arguments{
\item{map}{SpatRaster. Raster of probabilities/conductivities. Either a binary map or a probability
map the former being a special case of the later (see references)}

\item{contour}{an sfc_POLYGON that correspond to the contour of the raster. If the raster is
rectangular the default corresponds to the bounding box.}

\item{...}{propagated to \link{track_line}}
}
\description{
Drive along a contour line of the raster using \link{track_line} to find seeds that will be used
to vectorized the network
}
\examples{
\dontrun{
library(terra)
library(sf)

map <- system.file("extdata", "network2.tif", package = "vecnet")
map <- rast(map)

seeds = init_seeds(map, min_conductivity = 0.8)

plot(map, col = viridis::inferno(25))
plot(seeds, add = TRUE, col = "green", lwd = 3)
}
}
\references{
Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
(2022) Vectorial and topologically valid segmentation of forestry road networks from ALS data. In prep.\cr\cr
Jean-Romain Roussel, Jean-François Bourdon , Ilythia D. Morley , Nicholas C. Coops, Alexis Achim
(2022) Correction, update, and enhancement of vectorial forestry road maps using ALS data, a
pathfinder, and seven metrics. Journal of Applied Earth Observation and Geoinformation.
}
