#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <array>
#include <vector>
using namespace std;
using namespace Rcpp;

/*
 This file contains a series of tests 
 */

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
  
  arma::mat example_jacobian_stl(const double& t, const vector<double>& states, 
            const vector<double>& parameters, const vector<double>& forcings) {
    arma::mat output = arma::eye(states.size(), states.size());
    output = -parameters[0]*output;
    return output;
  }
  
  array<vector<double>, 2> example_dae_stl(const double& t, const vector<double>& states, 
            const vector<double>& derivatives, const vector<double>& parameters, 
            const vector<double>& forcings) {
    
    vector<double> residues(states.size());  
    residues[0] = -0.04*states[0] + 1e4*states[1]*states[2];
    residues[1] = -residues[0] - 3e7*states[1]*states[1] - derivatives[1];
    residues[0] = residues[0] - derivatives[0];
    residues[2] = states[0] + states[1] + states[2] - 1.0;
    vector<double> observed{forcings[0]};
    array<vector<double>, 2> output{residues, observed};
    //cout << t << " " << states[1] << " " << states[2] << derivatives[0] << " " << derivatives[1] <<  '\n';
    return output;
  }
  
};
