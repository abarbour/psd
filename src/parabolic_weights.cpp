#include <Rcpp.h>
using namespace Rcpp; // otherwise add Rcpp::

//' @rdname parabolic_weights
//' @export
// [[Rcpp::export]]
List parabolic_weights_rcpp(int ntap) {
    
    NumericVector kseq = abs( seq_len( ntap ) - 1 );
    
    NumericVector ksq = kseq * kseq;
    int K   = pow(ntap, 1);
    int K2  = pow(ntap, 2);
    //int K3  = pow(ntap, 3);
    
    // orig: w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
    //NumericVector w = ( K2 - ksq ) * 6 / ( 4 * K3  +  3 * K2  -  K );
    //NumericVector w = ( K2 - ksq ) / ( K * (4 * K  -  1) * (K  +  1) );
    NumericVector w = ( K2 - ksq ) * ( 1.5 / ( K * (K  -  0.25) * (K  +  1) ) );
    
    List Weights = List::create(
                              Named("ntap") = K,
                              Named("taper_seq") = kseq+1,
                              Named("taper_weights") = w
                              );
    
    return Weights;
}
