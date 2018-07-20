##' Random lognormal
##' @param n number of draws
##' @param mean
##' @param sd
##' @param lower lower bound
##' @param upper upper bound
##' @author Jason Cope
##' @export

rtlnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
{
  ret <- numeric()
  if (length(n) > 1)
    n <- length(n)
  while (length(ret) < n) {
    y <- rlnorm(n - length(ret), mean, sd)
    y <- y[y >= lower & y <= upper]
    ret <- c(ret, y)
  }
  stopifnot(length(ret) == n)
  ret
}