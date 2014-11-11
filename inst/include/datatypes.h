#ifndef __RcppSundialsDataType_h__

#define __RcppSundialsDataType_h__
#include <Rcpp.h>

// Signature of the model written in C++
typedef SEXP ode_in_Cpp(SEXP time, SEXP states, SEXP parameters, SEXP forcings);

// Struct that contains the data to run C++ models
struct data_Cpp {
  Rcpp::NumericVector parameters;
  Rcpp::List forcings_data;
  int neq;
  ode_in_Cpp* model;
};

// Struct that contains the data to run R models
struct data_R {
  Rcpp::NumericVector parameters;
  Rcpp::List forcings_data;
  int neq;
  Rcpp::Function model;
};

#endif