#include <Rcpp.h>
using namespace Rcpp;

//' @title Random Walk Metropolis sampler using Rcpp
//' @description A Random Walk Metropolis sampler using Rcpp
//' @useDynLib SC19085
//' @import Rcpp
//' @param sigma the standard deviation of the normal random increment with a zero mean 
//' @param x0 the initial iteration point
//' @param N the designated number of random numbers (including the initial point x0)
//' @return a list of length 3, x the random numbers of size \code{N}, k the number of rejection times, acceptance_rate the acceptance rate of the candidate points 
//' @examples
//' \dontrun{
//' laplace<-function(x) return(1/2*exp(-abs(x)))
//' lapC <- rw_Metro_C(2,25,2000);
//' plot(1:2000,lapC$x,type='l')
//' abline(h=c(-3*sqrt(2),3*sqrt(2)))
//' }
//' @export
// [[Rcpp::export]]

List rw_Metro_C(double sigma, double x0, int N) {
  
  NumericVector x(N);
  
  NumericVector u = runif(N,0,1);
  
  x[0]=x0;
  
  int k=0;
  
  
  for(int i=1;i<N;i++){
    
    double y= as<double>(rnorm(1,x[i-1],sigma));
    
    if (u[i]<=exp(-abs(y)+abs(x[i-1]))){
      x[i]=y;}
    
    //double t1=as<double>(laplace(y));
    // double t2=as<double>(laplace(x[i-1]));
    
    // if(u[i]<= t1/t2){
    //x[i]=y;}
    else {
      x[i] = x[i-1];
      k++;
    }
  }
  
  // double acceptance_rate = (double) (N-1-k)/(N-1);
  return List::create(
    _["x"] = x, 
    _["k"] = k,
    _["acceptance_rate"]= (double)(N-1-k)/(N-1)
    
  );
}

