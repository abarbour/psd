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

//' @rdname psdcore
//' @export
// [[Rcpp::export]]
List resample_fft_rcpp(ComplexVector fftz, IntegerVector tapers) {
  
  //
  // resample and reweight an fft estimates for a given number of tapers
  // using quadratic scaling in Barbour and Parker (2014)
  //
  
  // needs:
  //  - complex fft vector (assumes double-length)
  //  - integer tapers vector
  
  //Function warning("warning"); // use Rf_warning instead
  
  int nf, ne, ne2, nhalf, nfreq, m, m2, mleft1, mleft2, j1, j2, Kc, ki, ik, tapcap=1000;
  Rcomplex fdc;
  List bw;
  
  // !careful: double-length fft estimates
  // perhaps Ill add a logical flag later [ ]
  nf = fftz.size() / 2;
  
  // even, double, and half lengths
  ne = nf - (nf % 2);
  if (ne < nf){
    // this hasnt been tested [ ]
    Rf_warning("fft was not done on an even length series");
  }
  ne2 = 2 * ne;
  nhalf = ne / 2;
  
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
    Named("freq.inds") = Freqs,
    Named("k.capped") = K,
    Named("psd") = psd
    );
    return psd_out;
}

// end
