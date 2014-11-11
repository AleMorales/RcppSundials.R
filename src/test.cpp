#include <Rcpp.h>
using namespace Rcpp;

/*
 This file contains a series of tests 
 */

//' An example of model (unit test)
//' @export
// [[Rcpp::export]]
List example_model(double t, NumericVector states, NumericVector parameters, NumericVector forcings) {
   NumericVector derivatives = -states*parameters[0];
   return List::create(derivatives);
}
