// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_ctap_simple
IntegerVector rcpp_ctap_simple(IntegerVector tapvec, const int maxslope);
RcppExport SEXP _psd_rcpp_ctap_simple(SEXP tapvecSEXP, SEXP maxslopeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type tapvec(tapvecSEXP);
    Rcpp::traits::input_parameter< const int >::type maxslope(maxslopeSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_ctap_simple(tapvec, maxslope));
    return rcpp_result_gen;
END_RCPP
}
// modulo_floor
IntegerVector modulo_floor(IntegerVector n, int m);
RcppExport SEXP _psd_modulo_floor(SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(modulo_floor(n, m));
    return rcpp_result_gen;
END_RCPP
}
// parabolic_weights_rcpp
List parabolic_weights_rcpp(const int ntap);
RcppExport SEXP _psd_parabolic_weights_rcpp(SEXP ntapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type ntap(ntapSEXP);
    rcpp_result_gen = Rcpp::wrap(parabolic_weights_rcpp(ntap));
    return rcpp_result_gen;
END_RCPP
}
// resample_fft_rcpp
List resample_fft_rcpp(ComplexVector fftz, IntegerVector tapers, bool verbose, const bool dbl, const int tapcap);
RcppExport SEXP _psd_resample_fft_rcpp(SEXP fftzSEXP, SEXP tapersSEXP, SEXP verboseSEXP, SEXP dblSEXP, SEXP tapcapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ComplexVector >::type fftz(fftzSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tapers(tapersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type dbl(dblSEXP);
    Rcpp::traits::input_parameter< const int >::type tapcap(tapcapSEXP);
    rcpp_result_gen = Rcpp::wrap(resample_fft_rcpp(fftz, tapers, verbose, dbl, tapcap));
    return rcpp_result_gen;
END_RCPP
}
// parabolic_weights_rcpp2
arma::vec parabolic_weights_rcpp2(const int ntap);
RcppExport SEXP _psd_parabolic_weights_rcpp2(SEXP ntapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type ntap(ntapSEXP);
    rcpp_result_gen = Rcpp::wrap(parabolic_weights_rcpp2(ntap));
    return rcpp_result_gen;
END_RCPP
}
// resample_fft_rcpp2
List resample_fft_rcpp2(const arma::cx_vec& fftz, const arma::ivec& tapers, bool verbose, const bool dbl, const int tapcap);
RcppExport SEXP _psd_resample_fft_rcpp2(SEXP fftzSEXP, SEXP tapersSEXP, SEXP verboseSEXP, SEXP dblSEXP, SEXP tapcapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cx_vec& >::type fftz(fftzSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type tapers(tapersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type dbl(dblSEXP);
    Rcpp::traits::input_parameter< const int >::type tapcap(tapcapSEXP);
    rcpp_result_gen = Rcpp::wrap(resample_fft_rcpp2(fftz, tapers, verbose, dbl, tapcap));
    return rcpp_result_gen;
END_RCPP
}
// riedsid_rcpp
arma::vec riedsid_rcpp(const arma::mat& PSD, const arma::ivec& ntaper);
RcppExport SEXP _psd_riedsid_rcpp(SEXP PSDSEXP, SEXP ntaperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type PSD(PSDSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ntaper(ntaperSEXP);
    rcpp_result_gen = Rcpp::wrap(riedsid_rcpp(PSD, ntaper));
    return rcpp_result_gen;
END_RCPP
}
// parabolic_weights_field
arma::field<arma::vec> parabolic_weights_field(const int ntap);
RcppExport SEXP _psd_parabolic_weights_field(SEXP ntapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type ntap(ntapSEXP);
    rcpp_result_gen = Rcpp::wrap(parabolic_weights_field(ntap));
    return rcpp_result_gen;
END_RCPP
}
// resample_mvfft
List resample_mvfft(const arma::cx_mat& fftz, const arma::ivec& tapers, bool verbose, const bool dbl, const int tapcap);
RcppExport SEXP _psd_resample_mvfft(SEXP fftzSEXP, SEXP tapersSEXP, SEXP verboseSEXP, SEXP dblSEXP, SEXP tapcapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cx_mat& >::type fftz(fftzSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type tapers(tapersSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type dbl(dblSEXP);
    Rcpp::traits::input_parameter< const int >::type tapcap(tapcapSEXP);
    rcpp_result_gen = Rcpp::wrap(resample_mvfft(fftz, tapers, verbose, dbl, tapcap));
    return rcpp_result_gen;
END_RCPP
}
// det_vector
arma::cx_vec det_vector(const arma::cx_cube& x);
RcppExport SEXP _psd_det_vector(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cx_cube& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(det_vector(x));
    return rcpp_result_gen;
END_RCPP
}
// solve_tf
arma::cx_mat solve_tf(arma::cx_cube x);
RcppExport SEXP _psd_solve_tf(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_cube >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_tf(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_psd_rcpp_ctap_simple", (DL_FUNC) &_psd_rcpp_ctap_simple, 2},
    {"_psd_modulo_floor", (DL_FUNC) &_psd_modulo_floor, 2},
    {"_psd_parabolic_weights_rcpp", (DL_FUNC) &_psd_parabolic_weights_rcpp, 1},
    {"_psd_resample_fft_rcpp", (DL_FUNC) &_psd_resample_fft_rcpp, 5},
    {"_psd_parabolic_weights_rcpp2", (DL_FUNC) &_psd_parabolic_weights_rcpp2, 1},
    {"_psd_resample_fft_rcpp2", (DL_FUNC) &_psd_resample_fft_rcpp2, 5},
    {"_psd_riedsid_rcpp", (DL_FUNC) &_psd_riedsid_rcpp, 2},
    {"_psd_parabolic_weights_field", (DL_FUNC) &_psd_parabolic_weights_field, 1},
    {"_psd_resample_mvfft", (DL_FUNC) &_psd_resample_mvfft, 5},
    {"_psd_det_vector", (DL_FUNC) &_psd_det_vector, 1},
    {"_psd_solve_tf", (DL_FUNC) &_psd_solve_tf, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_psd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
