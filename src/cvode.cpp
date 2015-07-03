#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cvodes/cvodes.h>           // CVODES functions and constants
#include <nvector/nvector_serial.h>  // Serial N_Vector
#include <cvodes/cvodes_dense.h>     // CVDense
#include <datatypes.h>
#include <string> 
#include <limits> 
#include <array>
#include <vector>
#include <time.h>

using namespace Rcpp; 
using namespace std;
using arma::mat;
using arma::vec;

// [[Rcpp::interfaces(r, cpp)]]

/*
 *
 Functions when the model is written in C++ using the standard library and the Armadillo library
 *
 */

// Interface between cvode integrator and the model function written in Cpp
int cvode_to_Cpp_stl(double t, N_Vector y, N_Vector ydot, void* inputs) {
  // Cast the void pointer back to the correct data structure
  data_Cpp_stl* data = static_cast<data_Cpp_stl*>(inputs);
  // Interpolate the forcings
  vector<double> forcings(data->forcings_data.size());
  if(data->forcings_data.size() > 0) forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  vector<double> states(data->neq);
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i);
  // Run the model
  array<vector<double>, 2> output = data->model(t, states, data->parameters, forcings); 
  // Return the states to the NV_Ith_S
  vector<double> derivatives = output[0];
  for(auto i = 0; i < data->neq; i++)  NV_Ith_S(ydot,i) = derivatives[i];
  return 0; 
}

// Interface between cvode integrator and the jacobian function written in R
int cvode_to_Cpp_stl_jac(long int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, 
                   void *inputs, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  // Cast the void pointer back to the correct data structure
  data_Cpp_stl* data = static_cast<data_Cpp_stl*>(inputs);
  // Interpolate the forcings
  vector<double> forcings(data->forcings_data.size());
  if(data->forcings_data.size() > 0) forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  vector<double> states(data->neq);
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i);
  // Get the Jacobian
  mat output = data->jacobian(t, states, data->parameters, forcings);  
  // Return the DlsMat
  for(int j = 0; j < output.n_cols; j++) {
    for(int i= 0; i < output.n_rows; i++) {
      DENSE_ELEM(Jac, i, j) = output(i,j);
    }
  }
  return 0;
}

//' Returns the solution to the IVP problem at specified timepoints. 
//' Object returned is a matrix where nrow = number timepoints &
//' ncol = no. states + no. observed + [(no. states)x(no. parameters)]
//' the third term above is only returned when sensitivity equations are being
//' calculated
//' @export
// [[Rcpp::export]]
NumericMatrix wrap_cvodes(NumericVector times, NumericVector states_, 
                        NumericVector parameters_, List forcings_data_, 
                        List settings, SEXP model_, SEXP jacobian_) {
  // Wrap the pointer to the model function with the correct signature                        
  ode_in_Cpp_stl* model =  (ode_in_Cpp_stl *) R_ExternalPtrAddr(model_);
  // Wrap the pointer to the jacobian function with the correct signature                        
  jac_in_Cpp_stl* jacobian =  nullptr;
  if(as<int>(settings["jacobian"]) == 1) 
      jacobian = (jac_in_Cpp_stl *) R_ExternalPtrAddr(jacobian_);
  // Store all inputs in the data struct, prior conversion to stl and Armadillo classes
  int neq = states_.size();
  vector<double> parameters{as<vector<double>>(parameters_)};
  vector<double> states{as<vector<double>>(states_)};
  vector<mat> forcings_data(forcings_data_.size());
  if(forcings_data_.size() > 0) 
    for(int i = 0; i < forcings_data_.size(); i++)
      forcings_data[i] = as<mat>(forcings_data_[i]);
  data_Cpp_stl data_model{parameters, forcings_data, neq, model, jacobian};
  
  /*
   *
   Initialize CVODE and pass all initial inputs
   *
   */
  N_Vector y = nullptr;
  y = N_VNew_Serial(neq);
  for(int i = 0; i < neq; i++) {
    NV_Ith_S(y,i) = states[i];
  }
  void *cvode_mem = nullptr;
  if(as<std::string>(settings["method"]) == "bdf") {
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);     
  } else if(as<std::string>(settings["method"]) == "adams"){
    cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);     
  } else {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
    ::Rf_error("Please choose bdf or adams as method");
  }
  
  // Shut up Sundials (errors should not be printed to the screen)
  int flag = CVodeSetErrFile(cvode_mem, NULL);
  if(flag < 0.0) {
   if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
   if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}       
   ::Rf_error("Error in the CVodeSetErrFile function");  
  }
  
  // Initialize the Sundials solver. Here we pass initial N_Vector, the interface function and the initial time
  flag = CVodeInit(cvode_mem, cvode_to_Cpp_stl, times[0], y);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
    ::Rf_error("Error in the CVodeInit function"); 
  }
  
  // Tell Sundials the tolerance settings for error control
  Rcpp::NumericVector abstol = settings["atol"]; 
  if(abstol.size() > 1) {
    N_Vector Nabstol = nullptr;
    Nabstol = N_VNew_Serial(neq);
    for(int i = 0; i < neq; i++) {
      NV_Ith_S(Nabstol,i) = abstol[i];
    }
    flag = CVodeSVtolerances(cvode_mem, settings["rtol"], Nabstol);
  } else {
    flag = CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"]);    
  }

  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSStolerances function");
  }
  
  // Tell Sundials the number of state variables, so that I can allocate memory for the linear solver
  flag = CVDense(cvode_mem, neq);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVDense function");
  }

  // Give Sundials a pointer to the struct where all the user data is stored. It will be passed (untouched) to the interface as void pointer
  flag = CVodeSetUserData(cvode_mem, &data_model);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetUserData function");
  }
  
  // If we want to provide our own Jacobian, set the interface function to Sundials
  if(as<int>(settings["jacobian"]) == 1) {
    flag = CVDlsSetDenseJacFn(cvode_mem, cvode_to_Cpp_stl_jac);
    if(flag < 0.0) {
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
      ::Rf_error("Error in the CVDlsSetDenseJacFn function");
    }
  }

  // Set maximum number of steps
  flag = CVodeSetMaxNumSteps(cvode_mem, settings["maxsteps"]);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetUserData function");
  }
  
  // Set maximum order of the integration
  flag = CVodeSetMaxOrd(cvode_mem, settings["maxord"]); 
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetMaxOrd function");
  }
  
  // Set the initial step size
  flag = CVodeSetInitStep(cvode_mem, settings["hini"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetInitStep function");
  }
  
  // Set the minimum step size
  flag = CVodeSetMinStep(cvode_mem, settings["hmin"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetMinStep function");
  }
  
  // Set the maximum step size
  flag = CVodeSetMaxStep(cvode_mem, settings["hmax"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetMaxStep function");  
  }
  
  // Set the maximum number of error test fails
  flag = CVodeSetMaxErrTestFails(cvode_mem, settings["maxerr"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetMaxErrTestFails function");  
  }
  
  // Set the maximum number of nonlinear iterations per step
  flag = CVodeSetMaxNonlinIters(cvode_mem, settings["maxnonlin"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetMaxNonlinIters function");  
  }
  
  // Set the maximum number of convergence failures
  flag = CVodeSetMaxConvFails(cvode_mem, settings["maxconvfail"]);   
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetMaxConvFails function");
  }
  
  // Set stability limit detection
  flag = CVodeSetStabLimDet(cvode_mem, as<bool>(settings["stability"]));  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSetStabLimDet function");  
  }
  
  /*
   * 
   Make a first call to the model to check that everything is ok and retrieve the number of observed variables
   *
   */
  vector<double> forcings(forcings_data.size());
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[0]);
  std::array<vector<double>, 2> first_call;
  try {
    first_call = model(times[0], states, parameters, forcings); 
  } catch(std::exception &ex) {
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
      forward_exception_to_r(ex);
  } catch(...) { 
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  // Check length of time derivatives against the information passed through settings
  if(first_call[0].size() != neq) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
    ::Rf_error("Length of time derivatives returned by the model does not coincide with the number of state variables.");     
  }
  
  /*
   * 
   Fill up the output matrix with the values for the initial time
   *
   */
  vector<double> observed;
  int nder = 0;
  if(first_call.size() == 2) {
    vector<double> temp =  first_call[1];
    observed.resize(temp.size());
    observed = temp;
    nder = observed.size();
  }
  
  vector<int> extract_states = as<vector<int>>(settings["which_states"]);
  vector<int> extract_observed = as<vector<int>>(settings["which_observed"]);
  
  mat output(times.size(), extract_states.size() + extract_observed.size() + 1, arma::fill::zeros);
  output.at(0,0) = times[0];
  //for(auto i = 0; i < neq; i++) {
  //    output(0,i+1) = states[i];
  for(auto it = extract_states.begin(); it != extract_states.end(); it++) {
    if(*it > NV_LENGTH_S(y)) {
      Rcout << "The index " << *it << " exceeds the number of state variables of the model" << '\n';
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}        
      ::Rf_error("Simulation exited because of error in extracting state variables");
    }
    output.at(0,*it) = NV_Ith_S(y,*it - 1);
  }     

  if(extract_observed.size()  > 0) {
      //for(auto i = 0; i < nder; i++)  
      //    output(0,i + 1 + neq) = observed[i];
      for(auto it = extract_observed.begin(); it != extract_observed.end(); it++) {
        if(*it > observed.size()) {
          Rcout << "The index " << *it << " exceeds the number of observed variables returned by the model" << '\n';
          if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
          if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}            
          ::Rf_error("Simulation exited because of error in extracting observed variables");
        }
        output.at(0,*it + extract_states.size()) = observed[*it - 1];
      }     
  }
  
  
  
  /*
   * 
   Main time loop. Each timestep call cvode. Handle exceptions and fill up output
   *
   */

  double t = times[0];
  for(int i = 1; i < times.size(); i++) {
    try {
      flag = CVode(cvode_mem, times[i], y, &t, CV_NORMAL);
      if(as<bool>(settings["positive"])) {
        for(auto h = 0; h < neq; h++) {
         if(NV_Ith_S(y,h) < as<double>(settings["minimum"])) {
           Rcout << "The state variable at position " << h + 1 << " became smaller than minimum: " << NV_Ith_S(y,h) << " at time: " << times[i] << '\n';
           if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
           if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}  
           ::Rf_error("At least one of the states became smaller than minimum"); 
         }
        }
      }
      if(flag < 0.0) {
        switch(flag) {
          case -1:
            throw std::runtime_error("The solver took mxstep internal steps but could not reach tout."); break;
          case -2:
            throw std::runtime_error("The solver could not satisfy the accuracy demanded by the user for some internal step"); break;  
          case -3:
            throw std::runtime_error("Error test failures occured too many times during one internal time step or minimum step size was reached"); break;    
          case -4:
            throw std::runtime_error("Convergence test failures occurred too many times during one internal time step or minimum step size was reached."); break;  
          case -5:
            throw std::runtime_error("The linear solver’s initialization function failed."); break; 
          case -6:
            throw std::runtime_error("The linear solver’s setup function failed in an unrecoverable manner"); break; 
          case -7:
            throw std::runtime_error("The linear solver’s solve function failed in an unrecoverable manner"); break;  
          case -8:
            throw std::runtime_error("The right hand side function failed in an unrecoverable manner"); break;
          case -9:
            throw std::runtime_error("The right-hand side function failed at the first call."); break;    
          case -10:
            throw std::runtime_error("The right-hand side function had repeated recoverable errors."); break;   
          case -11:
            throw std::runtime_error("The right-hand side function had a recoverable errors but no recovery is possible."); break; 
          case -25:
            throw std::runtime_error("The time t is outside the last step taken."); break; 
          case -26:
            throw std::runtime_error("The output derivative vector is NULL."); break;   
          case -27:
            throw std::runtime_error("The output and initial times are too close to each other."); break;              
        }
      }
    } catch(std::exception &ex) {
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
      forward_exception_to_r(ex);
    } catch(...) { 
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
      ::Rf_error("c++ exception (unknown reason)"); 
    }

     // Write to the output matrix the new values of state variables and time
    output.at(i,0) = times[i];
    //for(auto h = 0; h < neq; h++) output(i,h + 1) = NV_Ith_S(y,h);
    for(auto it = extract_states.begin(); it != extract_states.end(); it++) {
      output.at(i,*it) = NV_Ith_S(y,*it - 1);
    }    
  }

  // If we have observed variables we call the model function again
  if(extract_observed.size() > 0 && flag >= 0.0) {
    for(unsigned i = 1; i < times.size(); i++) {
      // Get forcings values at time 0.
      if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[i]);
      // Get the simulate state variables
      for(auto j = 0; j < neq; j++) states[j] = output(i,j + 1);
      // Call the model function to retrieve total number of outputs and initial values for derived variables
      std::array<vector<double>, 2> model_call  = model(times[i], states, parameters, forcings); 
      observed =  model_call[1];
      // Derived variables already stored by the interface function
      //for(auto j = 0; j < nder; j++)  output(i,j + 1 + neq) = observed[j]; 
      for(auto it = extract_observed.begin(); it != extract_observed.end(); it++) {
        output.at(i,*it + extract_states.size()) = observed[*it - 1];
      }      
    } 
  }
              
  // De-allocate the N_Vector and the cvode_mem structures
  if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
  if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     

  return wrap(output);
}

//' Allows calling the model that calculates the time derivatives
//' @export
// [[Rcpp::export]]
List cvode_calc_derivs(SEXP model_, NumericVector t, NumericVector states, 
                          NumericVector parameters, List forcings_data_) {
   // Wrap the pointer to the model function with the correct signature                        
  ode_in_Cpp_stl* model =  (ode_in_Cpp_stl *) R_ExternalPtrAddr(model_); 
  // Interpolate the forcings
  vector<mat> forcings_data(forcings_data_.size());
  if(forcings_data_.size() > 0) 
    for(int i = 0; i < forcings_data_.size(); i++)
      forcings_data[i] = as<mat>(forcings_data_[i]);
  vector<double> forcings(forcings_data.size());
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, t[0]);
  // Call the model
  array<vector<double>, 2> output = model(t[0], as<vector<double>>(states),
                                          as<vector<double>>(parameters),
                                          forcings);
  // return the output as a list
  return List::create(_["Derivatives"] = wrap(output[0]),
                      _["Observed"] = wrap(output[1]));
}

//' Allows calling the function to calculate the Jacobian matrix of the model
//' @export
// [[Rcpp::export]]
NumericMatrix cvode_calc_jac(SEXP jacobian_, NumericVector t, NumericVector states, 
                          NumericVector parameters, List forcings_data_) {
   // Wrap the pointer to the model function with the correct signature                        
  jac_in_Cpp_stl* jacobian = (jac_in_Cpp_stl *) R_ExternalPtrAddr(jacobian_);
  // Interpolate the forcings
  vector<mat> forcings_data(forcings_data_.size());
  if(forcings_data_.size() > 0) 
    for(int i = 0; i < forcings_data_.size(); i++)
      forcings_data[i] = as<mat>(forcings_data_[i]);
  vector<double> forcings(forcings_data.size());
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, t[0]);
  // Call the model
  arma::mat output = jacobian(t[0], as<vector<double>>(states),
                                          as<vector<double>>(parameters),
                                          forcings);
  // return the output as a list
  return wrap(output);
} 
