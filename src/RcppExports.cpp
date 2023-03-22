// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// straw
Rcpp::DataFrame straw(std::string norm, std::string fname, std::string chr1loc, std::string chr2loc, const std::string& unit, int32_t binsize, std::string matrix);
RcppExport SEXP _trackViewer_straw(SEXP normSEXP, SEXP fnameSEXP, SEXP chr1locSEXP, SEXP chr2locSEXP, SEXP unitSEXP, SEXP binsizeSEXP, SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type norm(normSEXP);
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type chr1loc(chr1locSEXP);
    Rcpp::traits::input_parameter< std::string >::type chr2loc(chr2locSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type unit(unitSEXP);
    Rcpp::traits::input_parameter< int32_t >::type binsize(binsizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type matrix(matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(straw(norm, fname, chr1loc, chr2loc, unit, binsize, matrix));
    return rcpp_result_gen;
END_RCPP
}
// readHicChroms
Rcpp::DataFrame readHicChroms(std::string fname);
RcppExport SEXP _trackViewer_readHicChroms(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(readHicChroms(fname));
    return rcpp_result_gen;
END_RCPP
}
// readHicBpResolutions
Rcpp::NumericVector readHicBpResolutions(std::string fname);
RcppExport SEXP _trackViewer_readHicBpResolutions(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(readHicBpResolutions(fname));
    return rcpp_result_gen;
END_RCPP
}
// readHicNormTypes
Rcpp::CharacterVector readHicNormTypes(std::string fname);
RcppExport SEXP _trackViewer_readHicNormTypes(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(readHicNormTypes(fname));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_trackViewer_straw", (DL_FUNC) &_trackViewer_straw, 7},
    {"_trackViewer_readHicChroms", (DL_FUNC) &_trackViewer_readHicChroms, 1},
    {"_trackViewer_readHicBpResolutions", (DL_FUNC) &_trackViewer_readHicBpResolutions, 1},
    {"_trackViewer_readHicNormTypes", (DL_FUNC) &_trackViewer_readHicNormTypes, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_trackViewer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
