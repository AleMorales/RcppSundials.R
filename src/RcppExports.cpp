// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/RcppSundials.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// wrap_cvodes
NumericMatrix wrap_cvodes(NumericVector times, NumericVector states_, NumericVector parameters_, List forcings_data_, List settings, SEXP model_, SEXP jacobian_);
static SEXP RcppSundials_wrap_cvodes_try(SEXP timesSEXP, SEXP states_SEXP, SEXP parameters_SEXP, SEXP forcings_data_SEXP, SEXP settingsSEXP, SEXP model_SEXP, SEXP jacobian_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type states_(states_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters_(parameters_SEXP);
    Rcpp::traits::input_parameter< List >::type forcings_data_(forcings_data_SEXP);
    Rcpp::traits::input_parameter< List >::type settings(settingsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type jacobian_(jacobian_SEXP);
    __result = Rcpp::wrap(wrap_cvodes(times, states_, parameters_, forcings_data_, settings, model_, jacobian_));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP RcppSundials_wrap_cvodes(SEXP timesSEXP, SEXP states_SEXP, SEXP parameters_SEXP, SEXP forcings_data_SEXP, SEXP settingsSEXP, SEXP model_SEXP, SEXP jacobian_SEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(RcppSundials_wrap_cvodes_try(timesSEXP, states_SEXP, parameters_SEXP, forcings_data_SEXP, settingsSEXP, model_SEXP, jacobian_SEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// cvode_calc_derivs
List cvode_calc_derivs(SEXP model_, NumericVector t, NumericVector states, NumericVector parameters, List forcings_data_);
static SEXP RcppSundials_cvode_calc_derivs_try(SEXP model_SEXP, SEXP tSEXP, SEXP statesSEXP, SEXP parametersSEXP, SEXP forcings_data_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type states(statesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< List >::type forcings_data_(forcings_data_SEXP);
    __result = Rcpp::wrap(cvode_calc_derivs(model_, t, states, parameters, forcings_data_));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP RcppSundials_cvode_calc_derivs(SEXP model_SEXP, SEXP tSEXP, SEXP statesSEXP, SEXP parametersSEXP, SEXP forcings_data_SEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(RcppSundials_cvode_calc_derivs_try(model_SEXP, tSEXP, statesSEXP, parametersSEXP, forcings_data_SEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// cvode_calc_jac
NumericMatrix cvode_calc_jac(SEXP jacobian_, NumericVector t, NumericVector states, NumericVector parameters, List forcings_data_);
static SEXP RcppSundials_cvode_calc_jac_try(SEXP jacobian_SEXP, SEXP tSEXP, SEXP statesSEXP, SEXP parametersSEXP, SEXP forcings_data_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< SEXP >::type jacobian_(jacobian_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type states(statesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< List >::type forcings_data_(forcings_data_SEXP);
    __result = Rcpp::wrap(cvode_calc_jac(jacobian_, t, states, parameters, forcings_data_));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP RcppSundials_cvode_calc_jac(SEXP jacobian_SEXP, SEXP tSEXP, SEXP statesSEXP, SEXP parametersSEXP, SEXP forcings_data_SEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(RcppSundials_cvode_calc_jac_try(jacobian_SEXP, tSEXP, statesSEXP, parametersSEXP, forcings_data_SEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
// ida_Cpp_stl
NumericMatrix ida_Cpp_stl(NumericVector times, NumericVector states_, NumericVector derivatives_, NumericVector parameters_, List forcings_data_, List settings, SEXP model_, SEXP jacobian_);
RcppExport SEXP RcppSundials_ida_Cpp_stl(SEXP timesSEXP, SEXP states_SEXP, SEXP derivatives_SEXP, SEXP parameters_SEXP, SEXP forcings_data_SEXP, SEXP settingsSEXP, SEXP model_SEXP, SEXP jacobian_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type states_(states_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type derivatives_(derivatives_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters_(parameters_SEXP);
    Rcpp::traits::input_parameter< List >::type forcings_data_(forcings_data_SEXP);
    Rcpp::traits::input_parameter< List >::type settings(settingsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type jacobian_(jacobian_SEXP);
    __result = Rcpp::wrap(ida_Cpp_stl(times, states_, derivatives_, parameters_, forcings_data_, settings, model_, jacobian_));
    return __result;
END_RCPP
}
// ida_calc_res
List ida_calc_res(SEXP model_, NumericVector t, NumericVector states, NumericVector derivatives, NumericVector parameters, List forcings_data_);
RcppExport SEXP RcppSundials_ida_calc_res(SEXP model_SEXP, SEXP tSEXP, SEXP statesSEXP, SEXP derivativesSEXP, SEXP parametersSEXP, SEXP forcings_data_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type states(statesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type derivatives(derivativesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< List >::type forcings_data_(forcings_data_SEXP);
    __result = Rcpp::wrap(ida_calc_res(model_, t, states, derivatives, parameters, forcings_data_));
    return __result;
END_RCPP
}
// ida_calc_jac
NumericMatrix ida_calc_jac(SEXP jacobian_, NumericVector t, NumericVector states, NumericVector parameters, List forcings_data_);
RcppExport SEXP RcppSundials_ida_calc_jac(SEXP jacobian_SEXP, SEXP tSEXP, SEXP statesSEXP, SEXP parametersSEXP, SEXP forcings_data_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type jacobian_(jacobian_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type states(statesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< List >::type forcings_data_(forcings_data_SEXP);
    __result = Rcpp::wrap(ida_calc_jac(jacobian_, t, states, parameters, forcings_data_));
    return __result;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int RcppSundials_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("NumericMatrix(*wrap_cvodes)(NumericVector,NumericVector,NumericVector,List,List,SEXP,SEXP)");
        signatures.insert("List(*cvode_calc_derivs)(SEXP,NumericVector,NumericVector,NumericVector,List)");
        signatures.insert("NumericMatrix(*cvode_calc_jac)(SEXP,NumericVector,NumericVector,NumericVector,List)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP RcppSundials_RcppExport_registerCCallable() { 
    R_RegisterCCallable("RcppSundials", "RcppSundials_wrap_cvodes", (DL_FUNC)RcppSundials_wrap_cvodes_try);
    R_RegisterCCallable("RcppSundials", "RcppSundials_cvode_calc_derivs", (DL_FUNC)RcppSundials_cvode_calc_derivs_try);
    R_RegisterCCallable("RcppSundials", "RcppSundials_cvode_calc_jac", (DL_FUNC)RcppSundials_cvode_calc_jac_try);
    R_RegisterCCallable("RcppSundials", "RcppSundials_RcppExport_validate", (DL_FUNC)RcppSundials_RcppExport_validate);
    return R_NilValue;
}
