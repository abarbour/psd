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

//' @title Nearest value below
//' @export
//' @description 
//' Returns the nearest \code{m}-length value (downwards from \code{n}).
//' @details
//' This function is different from \code{\link{nextn}} in that the value is floored.
//' For example:
//' \code{10} is the result for \code{n=11,m=2} whereas \code{\link{nextn}} would give \code{12}.
//' 
//' @param n integer; the number of terms (can be a vector)
//' @param m integer; the modulo term (cannot be zero)
//' 
//' @author A.J. Barbour
//' 
//' @seealso \code{\link{psd-utilities}}; \code{\link{psdcore}} uses this to 
//' truncate series to their nearest even length (i.e., \code{m=2}).
//' 
//' @examples
//' n <- 11
//' nextn(n) # 12
//' modulo_floor(n) # 10
//' 
//' # works on vectors too:
//' # defaults to m=2
//' modulo_floor(seq_len(n))
//' #[1]  0  2  2  4  4  6  6  8  8 10 10
//' 
//' # change the floor factor
//' modulo_floor(seq_len(n), 3)
//' #[1] 0 0 3 3 3 6 6 6 9 9 9
//' 
//' # zeros are not allowed for m
//' try(modulo_floor(n, 0))
//' 
// [[Rcpp::export]]
IntegerVector modulo_floor(IntegerVector n, int m = 2){

  int i, ni, ntrunc, nn = n.size();
  IntegerVector ne(clone(n));
  
  if (m == 0) stop("m = 0  is invalid");
  
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
List parabolic_weights_rcpp(int ntap = 1) {
  //
  // return quadratic spectral weighting factors for a given number of tapers
  // Barbour and Parker (2014) Equation 7
  //

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

//' @title Resample an fft using varying numbers of sine tapers
//' 
//' @description
//' Produce an un-normalized psd based on an fft and a vector of optimal sine tapers
//' 
//' @details
//' To produce a psd estimate with our adaptive spectrum estimation method, we need only make one 
//' fft calculation initially and then
//' apply the weighting factors given by \code{\link{parabolic_weights_rcpp}}, which this
//' function does.
//' 
//' @param fftz complex; a vector representing the dual-length \code{\link{fft}}; see also the \code{dbl} argument
//' @param tapers integer; a vector of tapers
//' @param verbose logical; should messages be given?
//' @param dbl logical; should the code assume \code{fftz} is dual-length or single-length?
//' @param tapcap integer; the maximum number of tapers which can be applied; note that the length is
//' automatically limited by the length of the series.
//' 
//' @seealso \code{\link{riedsid}}
//' 
//' @examples
//' fftz <- complex(real=1:8, imaginary = 1:8)
//' taps <- 1:4
//' try(resample_fft_rcpp(fftz, taps))
//' 
//' @export
// [[Rcpp::export]]
List resample_fft_rcpp( ComplexVector fftz, IntegerVector tapers, 
  bool verbose = true, const bool dbl = true, const int tapcap=1000 ) {
  
  //
  // resample and reweight an fft estimates for a given number of tapers
  // using quadratic scaling in Barbour and Parker (2014)
  //
  //
  // needs:
  //  - fftz: complex vector -- the FFT of the original signal
  //  - tapers: integer vector -- the number of tapers at each frequency
  //
  // optional args:
  //  - verbose: logical -- should warnings be given?
  //  - dbl: logical -- should the progam assume 'fftz' is for a double-length (padded) series? 
  //                    Otherwise a single-length series is assumed.
  //  - tapcap: integer -- the maximum number of tapers at any frequency
  //
  
  int sc, nf, nt, ne, ne2, nhalf, nfreq, m, m2, mleft1, mleft2, j1, j2, Kc, ki, ik;
  
  if (dbl){
    // double-length fft estimates assumed by default
  	sc = 2;
  } else {
    // but could be single-length
	  sc = 1;
  }
  
  // even, double, and half lengths
  nf = fftz.size() / sc; 
  nt = tapers.size();
  ne = nf - (nf % 2);

  if (verbose){
    Function msg("message");
    msg(std::string("\tfft resampling"));
  }
  
  if (ne < nf){
    warning("fft was not done on an even length series");
  }
  
  ne2 = 2 * ne;
  nhalf = ne / 2;

  if (nhalf < 1){
    stop("cannot operate on length-1 series");
  }
  
  if (nt == 1){
    warning("forced taper length");
    tapers = rep(tapers, nhalf);
  }

  // Select frequencies for PSD evaluation [0:nhalf]
  NumericVector Freqs = abs(seq_len(nhalf)) - 1; // add one since c++ indexes at zero
  nfreq = Freqs.size();
  
  //
  // Calculate the psd by averaging over tapered estimates
  //

  IntegerVector K(nfreq);
  NumericVector absdiff(nfreq), psd(nfreq);
  
  double wi, cpsd;
  Rcomplex fdc;
  
  for (int j = 0; j < nfreq; j++) {
    
    m = Freqs[j];
    m2 = 2*m;
    // number of tapers applied at a given frequency (do not remove m+1 index!)
    Kc = tapers[j];
    Kc = tapers[m]; // orig: m+1, but this was leading to an indexing bug
    
    // set the current number of tapers, limited by a few factors
    if (Kc > nhalf){
      Kc = nhalf;
      if (Kc > tapcap){
        Kc = tapcap;
      }
    } else if (Kc <= 0){
      Kc = 1;
    }
    K[j] = Kc;
    
    IntegerVector k(Kc);
    NumericVector w(Kc);
    arma::rowvec sq_absdiff(Kc), psdprod(Kc);
    
    // taper sequence and spectral weights
    List bw = parabolic_weights_rcpp(Kc); // bw is a list of Kc, k, w
    k = bw[1];
    w = bw[2];
    
    // scan across ki to get vector of
    // spectral differences based on modulo indices
    for (ik = 0; ik < Kc; ik++) {
      
      wi = w[ik];
      ki = k[ik];
      
      mleft1 = m2 + ne2 - ki;
      mleft2 = m2 + ki;
      
      j1 = (mleft1 % ne2); // orig: + 1; // spectral indices
      j2 = (mleft2 % ne2); // orig: + 1;
      
      fdc = fftz[j1] - fftz[j2];
      
      // absolute value of fdc, squared:
      //    (sqrt(Re^2 + Im^2))^2 = Re^2 + Im^2
      sq_absdiff(ik) = pow(fdc.r, 2) + pow(fdc.i, 2);
        
      // at the same time to a scalar product (pre-sum)
      psdprod(ik) = wi * sq_absdiff(ik);
      
    }    
    // un-normalized psd is vector product of w and (sqrt(Re^2 + Im^2))^2
    cpsd = sum(psdprod);    
    psd[j] = cpsd;
  }
  
  List psd_out = List::create(
    Named("freq.inds") = Freqs + 1,
    Named("k.capped") = K,
    Named("psd") = psd
    );
    return psd_out;
}

/*** R 
n. <- 100
set.seed(1234)
# random walk
x <- cumsum(sample(c(-1, 1), n., TRUE))
fftz <- fft(x)
taps <- ceiling(runif(n.,10,30))
try(rsz <- resample_fft_rcpp(fftz, taps))
str(rsz)
layout(matrix(1:3))
par(oma=c(0.2,0.2,0.2,0.2), mar=c(3,4,1,0.1))
plot(x)
with(rsz,{
  plot(freq.inds, psd)
  plot(freq.inds, k.capped, type='h')
})
*/

// end all
