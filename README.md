![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)

# Vecnet

`vecnet` provides tools to vectorize a road network either from a binary map or a probability map derived from machine learning or any other method.  The algorithm is expected to be robust to gaps (false negative) and false positive. This package is experimental and must be considered as a proof of concept. It contains the implementation of the algorithm described in:

> Roussel, J., Bourdon, J., Morley, I. D., Coops, N. C., & Achim, A. (2023). Vectorial and topologically valid segmentation of forestry road networks from ALS data. International Journal of Applied Earth Observation and Geoinformation, 118, 103267. https://doi.org/10.1016/j.jag.2023.103267

The version [0.1.0](https://github.com/Jean-Romain/vecnet/releases/tag/v0.1.0) corresponds to the state of the software when submitted to the journal. Next versions include enhancements and bug fixes.

## Installation

``` r
remotes::install_github("Jean-Romain/vecnet")
```
