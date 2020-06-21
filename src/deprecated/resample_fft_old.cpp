

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
//' try(resample_fft_rcpp2(fftz, taps))
//' 
//' @export
// [[Rcpp::export]]
List resample_fft_rcpp2( const arma::cx_vec& fftz, 
                         const arma::ivec& tapers,
                         bool verbose = true, 
                         const bool dbl = true, 
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
  
  int sc, nf, nt, ne, ne2, nhalf, m2, mleft1, mleft2, Kc, ki, j1, j2;
  double wi;
  
  arma::cx_double fdc;
  
  
  if (dbl){
    // double-length fft estimates assumed by default
    sc = 2;
  } else {
    // but could be single-length
    sc = 1;
  }
  
  // even, double, and half lengths
  nf = fftz.n_elem / sc; 
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
  arma::vec psd(nhalf);
  psd.zeros();
  
  
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
  
  
  //
  // Calculate the psd by averaging over tapered estimates
  //
  
  for (int j = 0; j < nhalf; j++) {
    
    m2 = 2*j;
    // number of tapers applied at a given frequency (do not remove m+1 index!)
    
    Kc = taper_vec[j]; // orig: m+1, but this was leading to an indexing bug
    
    // taper sequence and spectral weights
    arma::vec w = parabolic_weights_rcpp2(Kc);
    
    // scan across ki to get vector of
    // spectral differences based on modulo indices
    for (arma::uword ik = 0; ik < Kc; ik++) {
      
      wi = w[ik];
      ki = ik + 1;
      
      mleft1 = m2 + ne2 - ki;
      mleft2 = m2 + ki;
      
      j1 = mleft1 % ne2;
      j2 = mleft2 % ne2;
      
      fdc = fftz[j1] - fftz[j2];
      
      // un-normalized psd is vector product of w and (sqrt(Re^2 + Im^2))^2
      // sum as loop progresses
      psd[j] += (pow(std::real(fdc), 2) +
        pow(std::imag(fdc), 2)) * wi;
      
    }
    
  }
  
  // return list to match previous function definition - jrk
  return Rcpp::List::create(
    Named("freq.inds") = arma::regspace<arma::vec>(1, nhalf),
    Named("k.capped") = taper_vec(arma::span(0, nhalf-1)),
    Named("psd") = psd
  );
}




//' @rdname parabolic_weights
//' @export
// [[Rcpp::export]]
arma::vec parabolic_weights_rcpp2(const int ntap = 1) {
  
  //
  // return quadratic spectral weighting factors for a given number of tapers
  // Barbour and Parker (2014) Equation 7
  //
  if(ntap == 0) {
    return(arma::zeros(1));
  }
  
  arma::vec kseq = arma::pow(arma::regspace<arma::vec>(0, ntap - 1), 2);
  
  
  double t3 = log(ntap * (ntap - 0.25) * (ntap + 1.0));
  
  arma::vec wgts = arma::exp(log(1.5) + arma::log(ntap * ntap - kseq) - t3);
  
  // return just the weights
  return wgts;
}