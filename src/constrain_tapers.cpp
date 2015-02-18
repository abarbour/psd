#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
IntegerVector constrain_tapers_rcpp(IntegerVector tapvec, int maxslope = 1) {
  //
  // Simple FIR filter to constrain variations in tapers
  //
  
  IntegerVector tapvec_c = clone(tapvec);
  
  bool state = true;
  int i, im, valc, valo, slope, newval;
  int ssize = tapvec_c.size();
  
  // indices (f - forward, r - reverse)
  int i_min_f = 0,         i_max_f = ssize - 1;
  int i_min_r = ssize - 1, i_max_r = 0;
  
  Rcout << "\nsize " << ssize << " max-slp " << maxslope << " " << i_max_f << " " << i_min_f << " " << i_max_r << " " << i_min_r << "\n\n";
   
  if (maxslope <= 0) Rf_error( "max slope must greater than zero" );
  
  // APPLY CONSTRAINTS
  //   FORWARDS:
  for (i = 1; i < ssize; i++){
    
    im = i - 1;
    valc = tapvec_c[i];
    valo = tapvec_c[im];
    slope = valc - valo;
    newval = valo + maxslope;
    
    Rcout << "F -- i " << i << "\tim " << im << "\tval-c: " << valc << "\tval-o: " << valo << "\tslope: " << slope << "\tnew " << newval << "\tste " << state << "\n";
    
    if (state){
      if (slope > maxslope) {
        tapvec_c[i] = newval;
        state = false;
      }
    } else {
      if (valc >= newval) {
        tapvec_c[i] = newval;
      } else {
        state = true;
      }
    }
    
  }
  
  Rcout << "\n";
  
  // and BACKWARDS:
  for (i = ssize; i > 0; i--){
    
    im = i - 1;
    valc = tapvec_c[i];
    valo = tapvec_c[im];
    slope = valc - valo;
    newval = valc + maxslope;
    
    Rcout << "R -- i " << i << "\tim " << im << "\tval-c: " << valc << "\tval-o: " << valo << "\tslope: " << slope << "\tnew " << newval << "\tste " << state << "\n";
    
    if (state){
      if (slope > maxslope) {
        tapvec_c[im] = newval;
        state = false;
      }
    } else {
      if (valc > newval) {
        tapvec_c[im] = newval;
      } else {
        state = true;
      }
    }
    
  }
  
  return tapvec_c;
}

/*** R

set.seed(1237)
n <- 11
x <- seq_len(n)
x2 <- x*2
xn <- round(runif(n,1,n))

xc <- constrain_tapers_rcpp(x)
xc2 <- constrain_tapers_rcpp(x, 2)

x2c <- constrain_tapers_rcpp(x2)
x2c2 <- constrain_tapers_rcpp(x2, 2)

xnc <- constrain_tapers_rcpp(xn)
xnc2 <- constrain_tapers_rcpp(xn, 2)

#plot(x, type='b', pch=16)
#abline(a=0,b=1, col='red', lty=3); abline(a=0,b=2, col='blue', lty=3)
#lines(xc, type='b', col='red')
#lines(xc2, type='b', col='blue')
#lines(0.2+as.vector(ctap_simple(as.tapers(x))), type='b', pch=".", col='salmon')

plot(xn, type='b', pch=16, ylim=c(0,12))
grid()
abline(a=0,b=1, col='red', lty=3); abline(a=0,b=2, col='blue', lty=3)
lines(xnc, type='b', col='red')
lines(xnc2, type='b', col='blue')
lines(0.2+as.vector(ctap_simple(as.tapers(xn))), type='b', pch=".", col='salmon')

#plot(x2, type='b', pch=16, ylim=c(0,25))
#abline(a=0,b=1, col='red', lty=3); abline(a=1,b=2, col='blue', lty=3)
#lines(x2c, type='b', col='red')
#lines(x2c2, type='b', col='blue')
#lines(0.2+as.vector(ctap_simple(as.tapers(x2))), type='b', pch=".", col='salmon')

*/
