// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
//   fft resampling with parabolic weighting, in c++ with Rcpp/RcppArmadillo 
//
//	  -- initial tests indicate upwards of a factor of 10 speed improvement vs. pure R implementations
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

#include <RcppArmadillo.h>
// do not #include <Rcpp.h> once Arma has been
using namespace Rcpp; // otherwise add Rcpp::

// [[Rcpp::depends(RcppArmadillo)]]

//' @rdname psd-utilities
//' @export
//' @description \code{modulo_floor} returns the nearest \code{m}-length
//' value (downwards); this is different from \code{\link{nextn}}.
//' For example:
//' \code{1001} goes to \code{1000} for \code{m=2}
//' @param n integer; the number of terms
//' @param m integer; the modulo term
// [[Rcpp::export]]
IntegerVector modulo_floor_rcpp(IntegerVector n, int m = 2){

  int i, ni, ntrunc;
  int nn = n.size();
  IntegerVector ne(nn);
  
  for (i = 0; i < nn; i++){
  	ni = n[i];
  	ntrunc = (ni % m);
  	ne[i] = ni - ntrunc;
  }
  
  return ne;
}

//' @rdname parabolic_weights
//' @export
// [[Rcpp::export]]
List parabolic_weights_rcpp(int ntap) {
  //
  // return quadratic spectral weighting factors for a given number of tapers
  // Barbour and Parker (2014) Equation 7
  //
  
  //int K, K2;
  int K = pow(ntap, 1);
  
  NumericVector kseq(ntap), wgts(ntap);
  kseq = abs( seq_len( K ) - 1 );
  
  // orig: w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
  //   or: w = ( K2 - ksq ) * 6 / ( 4 * K3  +  3 * K2  -  K );
  //   or: w = ( K2 - ksq ) / ( K * (4 * K  -  1) * (K  +  1) );
  
  wgts = exp(log(1.5) + log( K * K - kseq * kseq ) - log( K * (K - 0.25) * (K + 1) ));
  
  List weights_out = List::create(
    Named("ntap")=K,
    Named("taper_seq")=kseq + 1,
    Named("taper_weights")=wgts
    );
    
    return weights_out;
}

//' @title Resample an fft using varying tapers
//' @param fftz complex; a vector representing the dual-length \code{\link[fftw]{FFT}}; see also \code{dbl}.
//' @param tapers integer; a vector of tapers
//' @param verbose logical; should messages be given?
//' @param dbl logical; should the code assume \code{fftz} is dual-length or singl-length?
//' @param tapcap integer; the maximum number of tapers which can be applied; note that the length is
//' automatically limited by the length of the series.
//' @export
// [[Rcpp::export]]
List resample_fft_rcpp( ComplexVector fftz, IntegerVector tapers, 
  bool verbose = true, bool dbl = true, const int tapcap=1000 ) {
  
  //
  // resample and reweight an fft estimates for a given number of tapers
  // using quadratic scaling in Barbour and Parker (2014)
  //
  
  // needs:
  //  - complex fft vector (assumes double-length)
  //  - integer tapers vector
  
  //Function warning("warning"); // use Rf_warning instead
  
  int sc, nf, nt, ne, ne2, nhalf, nfreq, m, m2, mleft1, mleft2, j1, j2, Kc, ki, ik;
  Rcomplex fdc;
  List bw;
  
  // double-length fft estimates assumed by default
  if (dbl){
  	sc = 2;
  } else {
	  sc = 1;
  }
  
  nf = fftz.size() / sc;
  nt = tapers.size();

  // even, double, and half lengths
  ne = nf - (nf % 2);
  
  if (verbose){
    Function msg("message");
    msg(std::string("fft resampling"));
  }
  
  if (ne < nf){
    Rf_warning("fft was not done on an even length series");
  }
  
  ne2 = 2 * ne;
  nhalf = ne / 2;

  if (nhalf < 1){
    stop("cannot operate on length-1 series");
  }
  
  if (nt == 1){
    Rf_warning("forced taper length");
    tapers = rep(tapers, nhalf);
  }
  
  // %  Select frequencies for PSD evaluation [0:nhalf]
  NumericVector Freqs = abs(seq_len(nhalf)) - 1; // add one since c++ indexes at zero
  
  //%  Calculate the psd by averaging over tapered estimates
  nfreq = Freqs.size();
  
  //printf ("NF: %i %i\n", nf, nfreq);
  
  NumericVector K(nfreq), absdiff(nfreq), psd(nfreq);
  
  for (int j = 0; j < nfreq; j++) {
    
    m = Freqs[j];
    m2 = 2*m;
    Kc = tapers[m + 1]; // number of tapers applied at a given frequency
    //^^^ is it ok that I removed the (m + 1) index? [no, it's not!]
    
    if (Kc > nhalf){
      Kc = nhalf;
      if (Kc > tapcap){
        Kc = tapcap;
      }
    }
    K[j] = Kc;
    
    NumericVector k(Kc);
    arma::rowvec sq_absdiff(Kc);
    
    // taper sequence and spectral weights
    bw = parabolic_weights_rcpp(Kc); // bw is Kc, k, w
    k = bw[1];
    arma::rowvec w = bw[2];
    
    // scan across ki to get vector of
    // spectral differences based on modulo indices
    for (ik = 0; ik < Kc; ik++) {
      ki = k[ik];
      mleft1 = m2 + ne2 - ki;
      mleft2 = m2 + ki;
      
      j1 = (mleft1 % ne2) + 1; // spectral indices
      j2 = (mleft2 % ne2) + 1;
      
      fdc = fftz[j1] - fftz[j2];
      
      // absolute value of fdc, squared:
      //    (sqrt(Re^2 + Im^2))^2 = Re^2 + Im^2
      sq_absdiff(ik) = pow(fdc.r,2) + pow(fdc.i,2);   
    }    
    // un-normalized psd is vector product of w and (sqrt(Re^2 + Im^2))^2
    psd(j) = arma::as_scalar( w * sq_absdiff.t() );
    
  }
  
  List psd_out = List::create(
    Named("freq.inds") = Freqs + 1,
    Named("k.capped") = K,
    Named("psd") = psd
    );
    return psd_out;
}

// end
