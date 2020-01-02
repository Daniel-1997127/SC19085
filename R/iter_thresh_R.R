
#' @title Iterative thresholding method for the top eigen decomposition of a matrix  
#' @description Iterative thresholding method for calculating sparse unit eigenvector corresponding to the eigenvalue with the largest absolute value
#' @param S the matrix whose sparse unit eigenvector to be calculated
#' @param u0 initial unit vector
#' @param theta the sparese parameter which is also the truncation threshold
#' @param tol the tolerance of error 
#' @return a list of length 2, one is the eigenvector, the other is the estimation for the  (absolute) largest eigenvalue
#' @examples
#' \dontrun{
#' S<-matrix(runif(100),nrow=10)
#' u0<-c(rep(0,9),1)
#' iter_thresh_R(S,u=u0)
#' }
#' @export
iter_thresh_R<-function(S,u0,theta=0.1,tol=1e-4){
  
  u.old<-u0
  t<-S%*%u.old
  t<-t*(abs(t)>=theta)
  u.new<-t/sqrt(sum(t^2))
  k<-1
  while(sum((u.new-u.old)^2)>tol){
    cat('iteration:',k,'loss is:',sum((u.new-u.old)^2),'\n')
    u.old<-u.new
    t<-S%*%u.old
    t<-t*(abs(t)>=theta)
    u.new<-t/sqrt(sum(t^2))
    k<-k+1
  }

  eigen_value<-as.numeric(t(u.new)%*%S%*%u.new)
  return(list(eigen_vector=u.new, eigen_value= eigen_value))
  
}
