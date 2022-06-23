#' Moving average
#' @noRd
ma <- function(x, n = 3)
{
  if (anyNA(x))x <- zoo::na.approx(x, na.rm = FALSE)
  as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))
}
