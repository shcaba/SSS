##' Beta function
##' @param n Number of iteration
##' @param m mean value
##' @param s standard deviaiton
##' @param a Lower bound
##' @param b Upper bound 
##' @author Jason Cope
##' @export

rbeta.ab <- function(n, m, s, a, b)
{
  # calculate mean of corresponding standard beta dist
  mu.std <- (m-a)/(b-a)

  # calculate parameters of std. beta with mean=mu.std and sd=s
  alpha <- (mu.std^2 - mu.std^3 - mu.std*s^2) / s^2
  beta  <- (mu.std - 2*mu.std^2 + mu.std^3 - s^2 + mu.std*s^2) / s^2

  # generate n draws from standard beta
  b.std <- rbeta(n, alpha, beta)

  # linear transformation from beta(0,1) to beta(a,b)
  b.out <- (b-a)*b.std + a

  return(b.out)
}