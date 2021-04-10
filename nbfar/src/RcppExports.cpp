// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// nb_dev
double nb_dev(arma::mat Y, arma::mat MU, arma::vec Phi, arma::mat naind);
RcppExport SEXP _nbfar_nb_dev(SEXP YSEXP, SEXP MUSEXP, SEXP PhiSEXP, SEXP naindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type MU(MUSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type naind(naindSEXP);
    rcpp_result_gen = Rcpp::wrap(nb_dev(Y, MU, Phi, naind));
    return rcpp_result_gen;
END_RCPP
}
// get_sc
double get_sc(arma::mat X, arma::mat Y);
RcppExport SEXP _nbfar_get_sc(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sc(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// grad_eta_nb
arma::mat grad_eta_nb(arma::mat Y, arma::mat eta, arma::vec Phi);
RcppExport SEXP _nbfar_grad_eta_nb(SEXP YSEXP, SEXP etaSEXP, SEXP PhiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Phi(PhiSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_eta_nb(Y, eta, Phi));
    return rcpp_result_gen;
END_RCPP
}
// grad_mu_nb
arma::mat grad_mu_nb(const arma::mat& Y, const arma::mat& mu, const arma::vec& Phi);
RcppExport SEXP _nbfar_grad_mu_nb(SEXP YSEXP, SEXP muSEXP, SEXP PhiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Phi(PhiSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_mu_nb(Y, mu, Phi));
    return rcpp_result_gen;
END_RCPP
}
// grad_mu_nb_uv
arma::mat grad_mu_nb_uv(const arma::mat& Y, const arma::mat& mu, const arma::vec& Phi, arma::mat& d2l);
RcppExport SEXP _nbfar_grad_mu_nb_uv(SEXP YSEXP, SEXP muSEXP, SEXP PhiSEXP, SEXP d2lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type d2l(d2lSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_mu_nb_uv(Y, mu, Phi, d2l));
    return rcpp_result_gen;
END_RCPP
}
// update_mu_phi
arma::vec update_mu_phi(arma::mat Y, arma::mat mu, arma::vec Phi);
RcppExport SEXP _nbfar_update_mu_phi(SEXP YSEXP, SEXP muSEXP, SEXP PhiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Phi(PhiSEXP);
    rcpp_result_gen = Rcpp::wrap(update_mu_phi(Y, mu, Phi));
    return rcpp_result_gen;
END_RCPP
}
// update_mu_alpha
arma::vec update_mu_alpha(const arma::mat& Y, const arma::mat& mu, const arma::vec& Phi, const arma::mat& naind);
RcppExport SEXP _nbfar_update_mu_alpha(SEXP YSEXP, SEXP muSEXP, SEXP PhiSEXP, SEXP naindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type naind(naindSEXP);
    rcpp_result_gen = Rcpp::wrap(update_mu_alpha(Y, mu, Phi, naind));
    return rcpp_result_gen;
END_RCPP
}
// nbrrr_likelihood
arma::vec nbrrr_likelihood(const arma::mat& Y, const arma::mat& MU, const arma::mat& ETA, const arma::vec& Phi, const arma::mat& naind);
RcppExport SEXP _nbfar_nbrrr_likelihood(SEXP YSEXP, SEXP MUSEXP, SEXP ETASEXP, SEXP PhiSEXP, SEXP naindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type MU(MUSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ETA(ETASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type naind(naindSEXP);
    rcpp_result_gen = Rcpp::wrap(nbrrr_likelihood(Y, MU, ETA, Phi, naind));
    return rcpp_result_gen;
END_RCPP
}
// mySdiff
arma::uvec mySdiff(arma::uvec x, arma::uvec y);
RcppExport SEXP _nbfar_mySdiff(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(mySdiff(x, y));
    return rcpp_result_gen;
END_RCPP
}
// nbrrr_cpp
Rcpp::List nbrrr_cpp(arma::mat Y, arma::mat X0, int rnk, arma::vec cindex, arma::mat ofset, arma::mat Zini, arma::vec PhiIni, arma::mat Cini, Rcpp::List control, int msind, arma::mat naind);
RcppExport SEXP _nbfar_nbrrr_cpp(SEXP YSEXP, SEXP X0SEXP, SEXP rnkSEXP, SEXP cindexSEXP, SEXP ofsetSEXP, SEXP ZiniSEXP, SEXP PhiIniSEXP, SEXP CiniSEXP, SEXP controlSEXP, SEXP msindSEXP, SEXP naindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< int >::type rnk(rnkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ofset(ofsetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zini(ZiniSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type PhiIni(PhiIniSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cini(CiniSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type msind(msindSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type naind(naindSEXP);
    rcpp_result_gen = Rcpp::wrap(nbrrr_cpp(Y, X0, rnk, cindex, ofset, Zini, PhiIni, Cini, control, msind, naind));
    return rcpp_result_gen;
END_RCPP
}
// get_sv1
double get_sv1(arma::cube xyx, arma::vec ue, int q);
RcppExport SEXP _nbfar_get_sv1(SEXP xyxSEXP, SEXP ueSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type xyx(xyxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ue(ueSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sv1(xyx, ue, q));
    return rcpp_result_gen;
END_RCPP
}
// get_sv2
double get_sv2(const arma::mat& xyx, const arma::mat& Y, int q);
RcppExport SEXP _nbfar_get_sv2(SEXP xyxSEXP, SEXP YSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xyx(xyxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sv2(xyx, Y, q));
    return rcpp_result_gen;
END_RCPP
}
// get_sv
double get_sv(const arma::cube& xyx, const arma::vec& ue, int q, arma::uvec tem_uvec);
RcppExport SEXP _nbfar_get_sv(SEXP xyxSEXP, SEXP ueSEXP, SEXP qSEXP, SEXP tem_uvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type xyx(xyxSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ue(ueSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tem_uvec(tem_uvecSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sv(xyx, ue, q, tem_uvec));
    return rcpp_result_gen;
END_RCPP
}
// softThres
double softThres(double x, double lambda);
RcppExport SEXP _nbfar_softThres(SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(softThres(x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// softT
arma::vec softT(arma::vec x, arma::vec lambda);
RcppExport SEXP _nbfar_softT(SEXP xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(softT(x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// nzcount
int nzcount(arma::vec x);
RcppExport SEXP _nbfar_nzcount(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(nzcount(x));
    return rcpp_result_gen;
END_RCPP
}
// nbzerosol_cpp
Rcpp::List nbzerosol_cpp(arma::mat Y, arma::mat X0, arma::mat ofset, Rcpp::List control, int msind, arma::mat naind);
RcppExport SEXP _nbfar_nbzerosol_cpp(SEXP YSEXP, SEXP X0SEXP, SEXP ofsetSEXP, SEXP controlSEXP, SEXP msindSEXP, SEXP naindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X0(X0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ofset(ofsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type msind(msindSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type naind(naindSEXP);
    rcpp_result_gen = Rcpp::wrap(nbzerosol_cpp(Y, X0, ofset, control, msind, naind));
    return rcpp_result_gen;
END_RCPP
}
// nbfar_cpp
Rcpp::List nbfar_cpp(arma::mat Y, arma::mat Xm, int nlam, arma::vec cindex, arma::mat ofset, Rcpp::List initw, double Dini, arma::mat Zini, arma::vec PhiIni, arma::mat Uini, arma::vec Vini, double lmax, Rcpp::List control, int msind, arma::mat naind, int maxit, double epsilon);
RcppExport SEXP _nbfar_nbfar_cpp(SEXP YSEXP, SEXP XmSEXP, SEXP nlamSEXP, SEXP cindexSEXP, SEXP ofsetSEXP, SEXP initwSEXP, SEXP DiniSEXP, SEXP ZiniSEXP, SEXP PhiIniSEXP, SEXP UiniSEXP, SEXP ViniSEXP, SEXP lmaxSEXP, SEXP controlSEXP, SEXP msindSEXP, SEXP naindSEXP, SEXP maxitSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xm(XmSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ofset(ofsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type initw(initwSEXP);
    Rcpp::traits::input_parameter< double >::type Dini(DiniSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zini(ZiniSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type PhiIni(PhiIniSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Uini(UiniSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Vini(ViniSEXP);
    Rcpp::traits::input_parameter< double >::type lmax(lmaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type msind(msindSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type naind(naindSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(nbfar_cpp(Y, Xm, nlam, cindex, ofset, initw, Dini, Zini, PhiIni, Uini, Vini, lmax, control, msind, naind, maxit, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// cv_nbfar_cpp
Rcpp::List cv_nbfar_cpp(arma::mat Y, arma::mat Xm, int nlam, arma::vec cindex, arma::mat ofset, Rcpp::List initw, arma::mat Zini, arma::vec PhiIni, Rcpp::List xx, double lmax, Rcpp::List control, int maxit, double epsilon, int nfold);
RcppExport SEXP _nbfar_cv_nbfar_cpp(SEXP YSEXP, SEXP XmSEXP, SEXP nlamSEXP, SEXP cindexSEXP, SEXP ofsetSEXP, SEXP initwSEXP, SEXP ZiniSEXP, SEXP PhiIniSEXP, SEXP xxSEXP, SEXP lmaxSEXP, SEXP controlSEXP, SEXP maxitSEXP, SEXP epsilonSEXP, SEXP nfoldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xm(XmSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ofset(ofsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type initw(initwSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zini(ZiniSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type PhiIni(PhiIniSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< double >::type lmax(lmaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type nfold(nfoldSEXP);
    rcpp_result_gen = Rcpp::wrap(cv_nbfar_cpp(Y, Xm, nlam, cindex, ofset, initw, Zini, PhiIni, xx, lmax, control, maxit, epsilon, nfold));
    return rcpp_result_gen;
END_RCPP
}
// cv_nbfar_par
Rcpp::List cv_nbfar_par(arma::mat Y, arma::mat Xm, int nlam, arma::vec cindex, arma::mat ofset, Rcpp::List initw, arma::mat Zini, arma::vec PhiIni, Rcpp::List xx, double lmax, Rcpp::List control, int maxit, double epsilon, int nfold);
RcppExport SEXP _nbfar_cv_nbfar_par(SEXP YSEXP, SEXP XmSEXP, SEXP nlamSEXP, SEXP cindexSEXP, SEXP ofsetSEXP, SEXP initwSEXP, SEXP ZiniSEXP, SEXP PhiIniSEXP, SEXP xxSEXP, SEXP lmaxSEXP, SEXP controlSEXP, SEXP maxitSEXP, SEXP epsilonSEXP, SEXP nfoldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xm(XmSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cindex(cindexSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ofset(ofsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type initw(initwSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zini(ZiniSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type PhiIni(PhiIniSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< double >::type lmax(lmaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type control(controlSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type nfold(nfoldSEXP);
    rcpp_result_gen = Rcpp::wrap(cv_nbfar_par(Y, Xm, nlam, cindex, ofset, initw, Zini, PhiIni, xx, lmax, control, maxit, epsilon, nfold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nbfar_nb_dev", (DL_FUNC) &_nbfar_nb_dev, 4},
    {"_nbfar_get_sc", (DL_FUNC) &_nbfar_get_sc, 2},
    {"_nbfar_grad_eta_nb", (DL_FUNC) &_nbfar_grad_eta_nb, 3},
    {"_nbfar_grad_mu_nb", (DL_FUNC) &_nbfar_grad_mu_nb, 3},
    {"_nbfar_grad_mu_nb_uv", (DL_FUNC) &_nbfar_grad_mu_nb_uv, 4},
    {"_nbfar_update_mu_phi", (DL_FUNC) &_nbfar_update_mu_phi, 3},
    {"_nbfar_update_mu_alpha", (DL_FUNC) &_nbfar_update_mu_alpha, 4},
    {"_nbfar_nbrrr_likelihood", (DL_FUNC) &_nbfar_nbrrr_likelihood, 5},
    {"_nbfar_mySdiff", (DL_FUNC) &_nbfar_mySdiff, 2},
    {"_nbfar_nbrrr_cpp", (DL_FUNC) &_nbfar_nbrrr_cpp, 11},
    {"_nbfar_get_sv1", (DL_FUNC) &_nbfar_get_sv1, 3},
    {"_nbfar_get_sv2", (DL_FUNC) &_nbfar_get_sv2, 3},
    {"_nbfar_get_sv", (DL_FUNC) &_nbfar_get_sv, 4},
    {"_nbfar_softThres", (DL_FUNC) &_nbfar_softThres, 2},
    {"_nbfar_softT", (DL_FUNC) &_nbfar_softT, 2},
    {"_nbfar_nzcount", (DL_FUNC) &_nbfar_nzcount, 1},
    {"_nbfar_nbzerosol_cpp", (DL_FUNC) &_nbfar_nbzerosol_cpp, 6},
    {"_nbfar_nbfar_cpp", (DL_FUNC) &_nbfar_nbfar_cpp, 17},
    {"_nbfar_cv_nbfar_cpp", (DL_FUNC) &_nbfar_cv_nbfar_cpp, 14},
    {"_nbfar_cv_nbfar_par", (DL_FUNC) &_nbfar_cv_nbfar_par, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_nbfar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
