## ---- eval=FALSE--------------------------------------------------------------
#  
#    List rw_Metro_C(double sigma, double x0, int N, Function f) {
#  
#    NumericVector x(N);
#  
#    NumericVector u = runif(N,0,1);
#  
#    x[0]=x0;
#  
#    int k=0;
#  
#  
#    for(int i=1;i<N;i++){
#  
#      double y= as<double>(rnorm(1,x[i-1],sigma));
#      double t1=as<double>(f(y));
#      double t2=as<double>(f(x[i-1]));
#  
#      if(u[i]<= t1/t2){
#        x[i]=y;}
#      else {
#        x[i] = x[i-1];
#        k++;
#      }
#    }
#  
#    return List::create(
#      _["x"] = x,
#      _["k"] = k,
#      _["acceptance_rate"]= (double)(N-1-k)/(N-1)
#    );
#  }
#  

## ----eval=TRUE----------------------------------------------------------------
library(SC19085)
library(microbenchmark)
sigma<-2
x0<-5
N<-1000
f<-function(x) return(0.5*exp(-abs(x)))

 tm1 <- microbenchmark::microbenchmark(
   rwR = rw_Metro_R(sigma,x0,N,f),
   rwC = rw_Metro_C(sigma,x0,N)
   )
 
 print(summary(tm1)[,c(1,3,5,6)])

## ---- eval=FALSE--------------------------------------------------------------
#  
#  List Jaco_C(MatrixXd A,double tol){
#    //Environment base("package:base");
#    //Function atan=base["atan"];
#    //Function isSymmetric=global["isSymmetric"];
#    //Function stop=global["stop"];
#  
#    //stopifnot(isSymmetric(A));
#  
#    int n = A.rows();
#    MatrixXd V = MatrixXd::Identity(n,n);
#    MatrixXd D = MatrixXd::Identity(n,n);
#    MatrixXd U;
#  
#    double maxpq=0.0,phi;
#    int p,q;
#  
#    for(int i=0;i<n-1;i++){
#      for( int j=i+1;j<n;j++){
#        if(abs(A(i,j))>abs(maxpq)){
#          maxpq=A(i,j);
#           p=i;
#           q=j;
#        }
#      }
#    }
#  
#    while (abs(maxpq)>tol){
#  
#      maxpq=0;
#      phi=0.5*atan2(2*A(p,q),A(p,p)-A(q,q));
#  
#      U = MatrixXd::Identity(n,n);
#      U(p,p)=cos(phi);
#      U(q,q)=cos(phi);
#      U(p,q)=-sin(phi);
#      U(q,p)=sin(phi);
#      D=U.transpose()*A*U;
#      V=V*U;
#      A=D;
#  
#      for(int i=0;i<n-1;i++){
#        for(int j=i+1;j<n;j++){
#          if(abs(A(i,j))>abs(maxpq)){
#            maxpq=A(i,j);
#            p=i;
#            q=j;
#          }
#        }
#      }
#  
#  
#    }
#    VectorXd values = A.diagonal();
#    return List::create(
#      _["values"] = values,
#      _["vectors"] = V
#    );
#  
#  }
#  

## ----eval=TRUE----------------------------------------------------------------

A<-matrix(runif(100),nrow=10)
for(i in 1:9){
for(j in (i+1):10) A[i,j]<-A[j,i]
 }
JR <- Jaco_R(A,1e-5)
t(JR$vectors)%*%JR$vectors
JR$values
eigen(A)$values


## ---- eval=FALSE--------------------------------------------------------------
#  
#  List iter_thresh_C(MatrixXd S,MatrixXd u0,double theta=0.1,
#                     double tol=1e-4) {
#  
#    MatrixXd u_old = u0;
#  
#    MatrixXd t = S*u_old;
#  
#    int N = S.rows();
#  
#    for(int i = 0; i<N ;i++){
#      if(abs(t(i,0))<theta) t(i,0) = 0;
#    }
#  
#    double norm = t.norm();
#  
#    MatrixXd u_new = t/norm;
#  
#    int k=1;
#    double loss = (u_new-u_old).norm();
#  
#    while(loss>tol){
#      cout<<"iteration:"<< k <<" "<<"loss is:"<< loss <<endl;
#      u_old=u_new;
#      t=S*u_old;
#      for(int i=0; i<N ;i++){
#        if(abs(t(i,0))<theta) t(i,0)=0;
#      }
#      norm=t.norm();
#      u_new=t/norm;
#      k=k+1;
#      loss = (u_new-u_old).norm();
#    }
#  
#   VectorXd eigenvalue = u_new.transpose()*S*u_new;
#  
#    return List::create(
#      _["eigen_vector"] = u_new,
#      _["eigen_value"] = eigenvalue
#    );
#  
#  }

## ----eval=TRUE----------------------------------------------------------------

S<-matrix(runif(400),nrow=20)
S<-t(S)%*%S
u0<-matrix(c(rep(0,19),1))
ITR <- iter_thresh_R(S,u=u0)
ITR
eigen(S)$values[1]


