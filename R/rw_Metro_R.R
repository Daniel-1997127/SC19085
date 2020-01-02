#' @title  Random Walk Metropolis sampler using R
#' @description A Random Walk Metropolis sampler using R
#' @param sigma the standard deviation of the normal random increment with a zero mean 
#' @param x0 the initial iteration point
#' @param N the designated number of random numbers (including the initial point x0)
#' @param f the density of target distribution
#' @return a list of length 3, x the random numbers of size \code{N}, k the number of rejection times, acceptance_rate the acceptance rate of the candidate points 
#' @examples
#' \dontrun{
#' laplace<-function(x) return(1/2*exp(-abs(x)))
#' lapR <- rw_Metro_R(2,25,2000,laplace);
#' plot(1:2000,lapR$x,type='l')
#' abline(h=c(-3*sqrt(2),3*sqrt(2)))
#' }
#' @export

rw_Metro_R <- function(sigma, x0, N, f) {

  
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (f(y) / f(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  
  accept_rate<-(N-1-k)/(N-1)
  return(list(x=x, k=k, accept_rate=accept_rate))
}
