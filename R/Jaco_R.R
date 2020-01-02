#' @title  Jacobian method for calculating the eigenvalues and eigenvectors of a symmetric matrix  
#' @description Jacobian method is for calculating the eigenvalues and orthonormal eigenvectors for a symmetric matrix 
#' @param A a symmetric matrix
#' @param x0 the initial iteration point
#' @param tol the tolerance of error
#' @return a list of length 2: values a vector of the estimated eigenvalues, vectors a matrix with each column as the corresponding unit eigenvector
#' @examples
#' \dontrun{
#' S<-matrix(runif(400),nrow=20)
#' for(i in 1:19){
#' for(j in (i+1):20) S[i,j]<-S[j,i]
#' }
#' Jaco_R(S,tol=.Machine$double.eps^0.25)
#' }
#' @export

Jaco_R<-function(A,tol=.Machine$double.eps^0.25){
  
  stopifnot(isSymmetric(A))
  
  n<-dim(A)[1]
  V<-D<-diag(rep(1,n))
  maxpq<-0
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(abs(A[i,j])>abs(maxpq)){
        maxpq<-A[i,j]
        p<-i
        q<-j
      }
    }
  }
  
  while (abs(maxpq)>tol){
    maxpq<-0
    
    phi<-1/2*atan2(2*A[p,q],A[p,p]-A[q,q])
    
    U<-diag(rep(1,n))
    U[p,p]<-cos(phi)
    U[q,q]<-cos(phi)
    U[p,q]<--sin(phi)
    U[q,p]<-sin(phi)
    D<-t(U)%*%A%*%U
    V<-V%*%U
    A<-D
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if(abs(A[i,j])>abs(maxpq)){
          maxpq<-A[i,j]
          p<-i
          q<-j
        }
      }
    }
  }
  
  return(list(values=diag(A),vectors=V))
}

