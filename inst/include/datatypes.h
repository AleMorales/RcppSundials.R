#ifndef __RcppSundialsDataType_h__

#define __RcppSundialsDataType_h__
#include <Rcpp.h>
#include <vector>
#include <RcppArmadillo.h>

// Signature of the model written in C++
typedef SEXP ode_in_Cpp(SEXP time, SEXP states, SEXP parameters, SEXP forcings);
typedef std::array<std::vector<double>, 2> ode_in_Cpp_stl(const double& t, const std::vector<double>& states, 
            const std::vector<double>& parameters, const std::vector<double>& forcings);

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

// Struct that contains the data to run C++ (with stl) models
struct data_Cpp_stl {
  std::vector<double> parameters;
  const std::vector<arma::mat> forcings_data;
  const int neq;
  ode_in_Cpp_stl* model;
};


#endif