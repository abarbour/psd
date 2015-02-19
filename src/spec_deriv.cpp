#include <Rcpp.h>
using namespace Rcpp;

// this is actually much easier to imlement in pure R, so for now it remains as a placeholder

// [[Rcpp::export]]
NumericVector spec_deriv(NumericVector S, IntegerVector K) {
   
   int i, ns = S.size(), nt = K.size();
   NumericVector kopt(ns), dY(ns), d2Y(ns), Y(clone(S));
   
   if (nt==1){
     kopt = rep(K[0], ns);
   }
   
   Y = log(S);
   
   for (i = 0; i < ns; i++){
     dY[i] = Y[i] + kopt[i];
   }
   
   return dY;
}

/*** R
  
 set.seed(1235)
 n <- 20
 s <- runif(n, 1, 10)
 k <- 5
 spec_deriv(s, k)
 
 */
