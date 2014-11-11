#include <Rcpp.h>
#include <array>
#include <vector>
using namespace std;
using namespace Rcpp;

/*
 This file contains a series of tests 
 */

//' An example of model (unit test)
//' @export
// [[Rcpp::export]]
List example_model(double t, NumericVector states, NumericVector parameters, NumericVector forcings) {
   NumericVector derivatives = -states*parameters[0];
   NumericVector observed{forcings[0]};
   return List::create(derivatives, observed);
}

extern "C" {
  array<vector<double>, 2> example_model_stl(const double& t, const vector<double>& states, 
            const vector<double>& parameters, const vector<double>& forcings) {
      
    vector<double> derivatives(states.size());  
    for(int i = 0; i < states.size(); i++) {
      derivatives[i] = -states[i]*parameters[0];
    }
    vector<double> observed{forcings[0]};
    array<vector<double>, 2> output{derivatives, observed};
    return output;
  }
  
};
