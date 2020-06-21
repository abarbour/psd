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

#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR

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

//In C++ language we able to have set of overloaded functions under cmath such as:
//float       pow( float base, float exp ) ---- 2
//double      pow( double base, double exp ) ---- 3
//long double pow( long double base, long double exp ) ---- 4
//float       pow( float base, int iexp ) ---- 5
//double      pow( double base, int iexp ) ---- 6
//long double pow( long double base, int iexp ) ---- 7

//' @rdname parabolic_weights
//' @export
// [[Rcpp::export]]
List parabolic_weights_rcpp(const int ntap = 1) {
  //
  // return quadratic spectral weighting factors for a given number of tapers
  // Barbour and Parker (2014) Equation 7
  //
  
  NumericVector kseq(ntap), wgts(ntap);
  kseq = abs( seq_len( ntap ) - 1 );
  
  // orig: w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
  //   or: w = ( K2 - ksq ) * 6 / ( 4 * K3  +  3 * K2  -  K );
  //   or: w = ( K2 - ksq ) / ( K * (4 * K  -  1) * (K  +  1) );
  
  wgts = exp(log(1.5) + log( ntap * ntap - kseq * kseq ) - log( ntap * (ntap - 0.25) * (ntap + 1.0) ));
  
  List weights_out = List::create(
    Named("ntap")=ntap,
    Named("taper_seq")=kseq + 1.0,
    Named("taper_weights")=wgts
  );
  
  return weights_out;
}



// For all but short series this should be faster
//' @title parabolic_weights_field
//' @rdname parabolic_weights
//' 
//' @param ntap the maximum number of tapers
//' 
//' @export
//' 
// [[Rcpp::export]]
arma::field<arma::vec> parabolic_weights_field(const int ntap) {
  
  //
  // return quadratic spectral weighting factors for a given number of tapers
  // Barbour and Parker (2014) Equation 7
  //
  arma::field<arma::vec> f(ntap, 1);
  arma::vec kseq = arma::pow(arma::regspace<arma::vec>(0, ntap - 1), 2);
  
  double t3;
  for (int i = 1; i < ntap+1; i++) {
    t3 = log(i * (i - 0.25) * (i + 1.0));
    f(i-1,0) = arma::exp(log(1.5) + arma::log(i * i - kseq(arma::span(0, i-1))) - t3);
  }
  
  return(f);
  
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
                        bool verbose = true, const bool dbl = true, 
                        const int tapcap=1000 ) {
  
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
  NumericVector Freqs = abs( seq_len(nhalf) ) - 1; // add one since c++ indexes at zero
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
    } else if (Kc > tapcap){
      Kc = tapcap;
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
    Named("freq.inds") = Freqs + 1.0,
    Named("k.capped") = K,
    Named("psd") = psd
  );
  return psd_out;
}







arma::mat pad_data(const arma::mat& psd,
                   const int n,
                   const int n_max,
                   const double eps = 1e-78) {
  
  int nc = psd.n_cols;
  arma::mat ii(n - 1 + (n_max) * 2, nc);
  
  ii.head_rows(n_max) = arma::reverse(psd.head_rows(n_max));
  ii.tail_rows(n_max+1) = arma::reverse(psd.tail_rows(n_max+1));
  ii(arma::span(n_max-1, n_max + n-2), arma::span::all) = psd;
  
  return(ii + eps);
  
}

//' @title replaces time consuming portion of riedsid2
//' 
//' @param PSD vector or class \code{'amt'} or \code{'spec'}; the spectral values used to optimize taper numbers
//' @param ntaper scalar or vector; number of tapers to apply optimization
//' @param riedsid_column scalar integer; which column to use in multivariate optimization. If the value is 0 the maximum number of tapers for all columns is chosen. If the value is < 0 the minimum number of tapers for all columns is chosen. If the value is 1, 2, 3, etc. the number of tapers is based on the column selected.
//' 
//' @return kopt vector
//' @export
//' 
// [[Rcpp::export]]
arma::vec riedsid_rcpp(const arma::mat& PSD,
                       const arma::ivec& ntaper,
                       int riedsid_column = 0){
  
  double eps = 1e-78;
  double sc = 473.3736;
  double uzero, L, L2, CC;
  int nf = PSD.n_rows;
  int nc = PSD.n_cols;
  int nt = ntaper.n_elem;
  int j1, j2;
  
  arma::vec out(nf);
  arma::ivec taper_vec(nf);
  
  // check for single length taper
  if (nt == 1){
    taper_vec.fill(ntaper(0));
  } else {
    taper_vec = ntaper;
  }
  
  arma::ivec nspan(nf);
  
  for (int j = 0; j < nf; j++) {
    nspan[j] = std::ceil(std::min(7.0 * taper_vec(j) / 5.0, nf / 2.0));
  }
  
  int nadd = 1 + arma::max(nspan);
  arma::mat y = arma::log(pad_data(PSD, nf, nadd, eps = 1e-78));
  
  double dy, d2y;
  arma::mat yders(nf, nc);
  
// determine columns to use for riedsid
  int i_start, i_end;
  if (riedsid_column <= 0) {
    i_start = 0;
    i_end = nc;
  } else if (riedsid_column > nc) {
    Rcpp::warning("riedsid_column is greater than the input number of columns. Setting to the maximum column number");
    i_start = nc-1;
    i_end = nc;
  } else {
    i_start = riedsid_column - 1;
    i_end = riedsid_column;
  }
  
  
  for (int j = 0; j < nf; j++) {
    
    j1 = j - nspan[j] + nadd - 1;
    j2 = j + nspan[j] + nadd - 1;
    arma::vec u = arma::regspace<arma::vec>(j1+1, j2+1) - (j1+1 + j2+1) / 2;
    
    L = j2 - j1 + 1;
    L2 = L * L;
    CC = 12;
    
    uzero = (L2 - 1) / CC;
    
    
    for (int i = i_start; i < i_end; i++) {
      
      // first deriv
      dy = as_scalar(u.t() * y(arma::span(j1, j2), i) * CC / (L*(L2 - 1)));
      
      // second deriv
      d2y = as_scalar((u % u - uzero).t() *
        y(arma::span(j1, j2), i) * 360 / (L * (L2 - 1) * (L2 - 4)));
      
      yders(j, i) = fabs(dy*dy + d2y + eps);
    }
    
    if (riedsid_column == 0) {
      out(j) = arma::max(yders.row(j));
    } else if (riedsid_column < 0) {
      out(j) = arma::min(yders.row(j));
    } else {
      out(j) = yders(j, i_start);
    }
    
  }
  
  return(arma::round(pow(sc, 0.2) / arma::pow(out, 0.4)));
  
}





// This is for multivariate series  --------------------------------------------



//' @title Resample an fft using varying numbers of sine tapers
//' 
//' @description
//' Produce an un-normalized psd based on an fft and a vector of optimal sine 
//' tapers.
//' 
//' @details
//' To produce a psd estimate with our adaptive spectrum estimation method,
//' we need only make one fft calculation initially and then apply the weighting
//' factors given by \code{\link{parabolic_weights}}, which this function 
//' does.
//' 
//' @param fftz complex; a matrix representing the dual-length \code{\link{fft}}; see also the \code{dbl} argument
//' @param tapers integer; a vector of tapers
//' @param verbose logical; should messages be given?
//' @param dbl logical; should the code assume \code{fftz} is dual-length or single-length?
//' @param tapcap integer; the maximum number of tapers which can be applied; note that the length is
//' automatically limited by the length of the series.
//' 
//' @return list that includes the auto and cross-spectral density, and the
//' number of tapers
//' 
//' @seealso \code{\link{riedsid}}
//' 
//' @examples
//' fftz <- complex(real=1:8, imaginary = 1:8)
//' taps <- 1:4
//' try(resample_mvfft(fftz, taps))
//' 
//' @export
// [[Rcpp::export]]
List resample_mvfft( const arma::cx_mat& fftz, 
                     const arma::ivec& tapers,
                     bool verbose = true, 
                     const bool dbl = true, 
                     const int tapcap = 10000 ) {
  
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
  
  int sc, nf, nr, nc, nt, ne, ne2, nhalf, m2, mleft1, mleft2, Kc, ki, j1, j2;
  double wi;
  
  
  arma::cx_double fdc1, fdc2;
  
  
  if (dbl){
    // double-length fft estimates assumed by default
    sc = 2;
  } else {
    // but could be single-length
    sc = 1;
  }
  
  // even, double, and half lengths
  nr = fftz.n_rows;
  nf = fftz.n_rows / sc; 
  nc = fftz.n_cols;
  nt = tapers.n_elem;
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
  
  arma::ivec taper_vec(nhalf);
  
  if (nhalf < 1){
    stop("cannot operate on length-1 series");
  }
  
  
  if (nt == 1){
    warning("forced taper length");
    taper_vec.fill(tapers[0]);
  } else {
    taper_vec = tapers;
  }
  
  
  // set the current number of tapers, limited by a few factors
  arma::uvec wh = arma::find(taper_vec > nhalf);
  taper_vec(wh).fill(nhalf);
  wh = arma::find(taper_vec > tapcap);
  taper_vec(wh).fill(tapcap);
  wh = arma::find(taper_vec <= 0);
  taper_vec(wh).ones();
  
  
  int mm = taper_vec.max();
  
  arma::vec w(mm);
  arma::field<arma::vec> para = parabolic_weights_field(mm);
  
  
  //
  // Calculate the psd by averaging over tapered estimates
  //
  
  
  arma::cx_cube psd(nhalf, nc, nc);
  psd.fill(arma::cx_double(0.0, 0.0));
  
  // each frequency
  for (int jj = 0; jj < nhalf; jj++) {
    
    m2 = 2 * jj;
    // number of tapers applied at a given frequency (do not remove m+1 index!)
    
    Kc = taper_vec[jj]; // orig: m+1, but this was leading to an indexing bug
    
    // taper sequence and spectral weights
    w = para(Kc-1, 0);
    
    // scan across ki to get vector of
    // spectral differences based on modulo indices
    for (int ik = 0; ik < Kc; ik++) {
      
      wi = w[ik];
      ki = ik + 1;
      
      mleft1 = m2 + ne2 - ki;
      mleft2 = m2 + ki;
      
      j1 = mleft1 % ne2;
      j2 = mleft2 % ne2;
      
      // sum as loop progresses, auto and cross-spectra
      
      for (int ii = 0; ii < nc; ii++) {
        for (int kk = ii; kk < nc; kk++) {
            
            fdc1 = fftz(j1, ii) - fftz(j2, ii);
            fdc2 = fftz(j1, kk) - fftz(j2, kk);
          
          psd(jj, ii, kk) += (fdc1 * std::conj(fdc2)) * wi;
        }
      }
    }
    
  }
  
  // Add lower triangle for use in transfer function calculation
  for (int ii = 0; ii < nc; ii++) {
    for (int kk = 0; kk < nc; kk++) {
      if (kk < ii){
        psd(arma::span(), arma::span(ii), arma::span(kk)) = 
          arma::conj(psd(arma::span(), arma::span(kk), arma::span(ii)));
      }
    }
  }
  
  // return list to match previous function definition - jrk
  return Rcpp::List::create(
    Named("freq.inds") = arma::regspace<arma::ivec>(1, nhalf),
    Named("k.capped") = taper_vec(arma::span(0, nhalf-1)),
    Named("psd") = psd
  );
}



// solve for transfer function


//' @title
//' det_vector
//'
//' @description
//' Determinant for an array
//'
//' @param x \code{numeric array} values to evaluate
//'
//' @return vector of determinants
//'
//'
//' @export
// [[Rcpp::export]]
arma::cx_vec det_vector(const arma::cx_cube& x) {
  
  std::size_t n = x.n_rows;
  std::size_t nc = x.n_cols;
  
  arma::cx_vec output(n);
  arma::cx_mat mat(nc, nc);
  
  for (std::size_t i = 0; i < n; i++){
    mat = x(arma::span(i), arma::span::all, arma::span::all);
    output(i) = arma::det(mat);
  }
  
  return(output);
  
}



// [[Rcpp::export]]
arma::cx_mat solve_tf(arma::cx_cube x) {
  
  unsigned int n  = x.n_rows;
  unsigned int nr = x.n_cols -1;
  
  arma::cx_mat output(n, nr);
  arma::cx_cube numer_base = x(arma::span::all,
                               arma::span(1, nr),
                               arma::span(1, nr));;
  arma::cx_cube numer(n, nr, nr);
  
  
  arma::cx_cube numer_sol = x(arma::span::all,
                              arma::span(1, nr),
                              arma::span(0));
  
  arma::cx_vec denom = det_vector(numer_base);
  
  for (arma::uword i=0; i < nr; i++){
    
    // fill new column
    numer = numer_base;
    numer(arma::span::all, arma::span::all, arma::span(i)) = numer_sol;
    
    output.col(i) = det_vector(numer) / denom;
    
  }
  
  return(output);
  
}




/*** R 

library(microbenchmark)
library(psd)

n. <- 1000000
set.seed(1234)
nc <- 2
x <- matrix(cumsum(sample(c(-1, 1), n., TRUE)), ncol=nc)
fftz <- mvfft(x)
psd <- Re(fftz * Conj(fftz))
taps <- ceiling(runif(n./nc,10,300))

microbenchmark(
  a1 <- riedsid_rcpp(as.matrix(psd)[, 1, drop = FALSE], taps),
  a2 <- riedsid_rcpp(psd[, 2, drop = FALSE], taps),
  b  <- riedsid_rcpp(as.matrix(psd), taps),
  times = 2
)

all.equal(apply(cbind(a1,a2), 1, min), b)


microbenchmark(
  rsz1 <- resample_fft_rcpp(fftz[,2], taps, verbose = FALSE),
  # rsz2 <- resample_fft_rcpp2(fftz[,2], taps, verbose = FALSE),
  rsz3 <- resample_mvfft(fftz, taps, verbose = FALSE),
  times = 2
)

all.equal(rsz1, rsz2)
all.equal(rsz1$psd, Re(rsz3$psd[,2,2]))


# 
# microbenchmark(
#   tf <- solve_tf(rsz3$psd),
#   times = 5
# )


# 
# microbenchmark(
#   kopt1 <- riedsid2(PSD = rsz1$psd, ntaper = taps, constrained = FALSE, verbose = FALSE),
#   kopt2 <- riedsid2(PSD = rsz1$psd, ntaper = taps, constrained = FALSE, verbose = FALSE, fast = TRUE),
#   times = 10
# )


# 
# all.equal(kopt1, kopt2)
# 
# 
# microbenchmark(
#   psd1 <- psdcore(x, verbose = FALSE),
#   psd2 <- psdcore(x, verbose = FALSE, fast = TRUE),
#   times = 10
# )
# 
# 
# all.equal(psd1, psd2)
# 
# 
# 
# microbenchmark(
#   pspec1 <- pspectrum(x, verbose = FALSE),
#   pspec2 <- pspectrum(x, verbose = FALSE, fast = TRUE),
#   times = 10
# )
# 
# all.equal(pspec1, pspec2)
# 
# 
# microbenchmark(
#   pilot1 <- pilot_spec(x, verbose = FALSE),
#   pilot2 <- pilot_spec(x, verbose = FALSE, fast = TRUE),
#   times = 10
# )
# 
# all.equal(pilot1, pilot2)


*/

