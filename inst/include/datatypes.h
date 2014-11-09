#ifndef __RcppSundialsDataType_h__

#define __RcppSundialsDataType_h__
#include <Rcpp.h>

// Signature of the model written in C++
typedef Rcpp::List ode_in_Cpp(double time, Rcpp::NumericVector states, 
                        Rcpp::NumericVector parameters, Rcpp::NumericVector forcings);

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