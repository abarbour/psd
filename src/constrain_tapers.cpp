// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
//   Taper constraints based on spectral derivatives
//
//    -- initial tests indicate upwards of a factor of 10 speed improvement vs. pure R implementations
//
//   Copyright (C) 2015  Andrew J. Barbour *
//
//   * Robert L. Parker authored the original algorithm.
//
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <Rcpp.h>
using namespace Rcpp;

// http://gallery.rcpp.org/articles/reversing-a-vector/
// [[Rcpp::export]]
IntegerVector irev(IntegerVector x) {
   IntegerVector revX = clone<IntegerVector>(x);
   std::reverse(revX.begin(), revX.end());
   ::Rf_copyMostAttrib(x, revX); 
   return revX;
}

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
  int ssize = tapvec_c.size();
  int i, im, valc, valo, slope, newval, dsize = ssize - 1;

  IntegerVector fdiff = diff( tapvec_c[seq(0,dsize)] );
  IntegerVector rdiff = diff( irev( tapvec_c[seq(1,ssize)] ) );

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
  
  return fdiff; //tapvec_c;
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
