// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
//   Apply constraints on tapers using simple derivatives
//
//     Copyright (C) 2017  Andrew J. Barbour *
//
//     * Robert L. Parker authored the original matlab algorithm
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

//' c++ implementation of the RLP constraint filter
//' 
//' @inheritParams tapers-refinement
//' @export
// [[Rcpp::export]]
IntegerVector rcpp_ctap_simple(IntegerVector tapvec, const int maxslope = 1) {
  //
  // Simple FIR filter to constrain variations in tapers
  // based on spectral derivatives so that we curb run-away growth of tapvec 
  // due to zeros in the psd''; limits slopes to be <= maxslope in magnitude, 
  // preventing greedy averaging.
  //
  
  if (maxslope < 0) stop( "max slope cannot be less than zero" );
  
  bool state = true;
  IntegerVector koptc(clone(tapvec));
  const int ssize = tapvec.size();
  int lasti = ssize - 1;
  int i, im, k_prev, k_curr, k_repl, slope;
  
  //  Scan FORWARD and create new series where slopes <= 1
  state = true;
  // orig: 2:nf
  for (i = 1; i <= lasti; i++){
    im = i - 1;
    k_prev = koptc[ im ];
    k_curr = koptc[ i ];
    k_repl = k_prev + maxslope;
    slope = k_curr - k_prev;
    if (state){
      if (slope >= maxslope){
        koptc[ i ] = k_repl;
        state = false;
      }
    } else {
      if (k_curr >= k_repl){
        koptc[ i ] = k_repl;
      } else {
        state = true;
      }
    }
  }
  // scan BACKWARDS to bound slopes >= -1
  state = true;
  //  orig: nf:-1:2
  for (i = lasti; i >= 1; i--){
    im = i - 1;
    k_prev = koptc[ i ];
    k_curr = koptc[ im ];
    k_repl = k_prev + maxslope;
    slope = k_curr - k_prev;
    if (state){
      if (slope >= maxslope){
        koptc[ im ] = k_repl;
        state = false;
      }
    } else {
      if (k_curr >= k_repl){
        koptc[ im ] = k_repl;
      } else {
        state = true;
      }
    }
  }
  
  // never let the the number of tapers be
  //  - less than one, or 
  //  - greater than the length of the series
  koptc = pmin(pmax(koptc, 1), ssize);
  
  return koptc;
}

//IntegerVector test_sample_indices( IntegerVector x ) {
//  int ssize = x.size();
//  IntegerVector outv(4);
//  outv[0] = x[0];
//  outv[1] = x[1];
//  outv[2] = x[ssize - 1];
//  outv[3] = x[ssize];
//  return outv;
//}
  
/*** R

ctap_simple_rcpp.default(1:10, 0)
ctap_simple_rcpp.default(1:10, -1)

#test_sample_indices(1:10)

taps <- c(100,rep(1,10),100,-1)
#matplot(t(rbind(taps, ctap_simple_rcpp(taps))), type='b')

taps <- c(100,rep(1,10),0,-1)
#matplot(t(rbind(taps, ctap_simple_rcpp(taps))), type='b')

taps <- c(100,round(runif(10,2,10)),0,-1)
#matplot(t(rbind(taps, ctap_simple_rcpp(taps))), type='b')

*/

