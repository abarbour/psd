#include <Rcpp.h>
using namespace Rcpp; // otherwise add Rcpp::

//' @rdname parabolic_weights
//' @export
// [[Rcpp::export]]
List parabolic_weights_rcpp(int ntap) {
    //
    // return quadratic spectral weighting factors for a given number of tapers
    // Barbour and Parker (2014) Equation 7
    //
    NumericVector kseq(ntap), w(ntap);
    int K = pow(ntap, 1);
    kseq = abs( seq_len( K ) - 1 );
    
    // orig: w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
    //   or: w = ( K2 - ksq ) * 6 / ( 4 * K3  +  3 * K2  -  K );
    //   or: w = ( K2 - ksq ) / ( K * (4 * K  -  1) * (K  +  1) );
    
    w = ( K * K  -  kseq * kseq ) * ( 1.5 / ( K * (K  -  0.25) * (K  +  1) ) );
    
    List Weights = List::create(
                                Named("ntap") = K,
                                Named("taper_seq") = kseq + 1,
                                Named("taper_weights") = w
                                );
    
    return Weights;
}
