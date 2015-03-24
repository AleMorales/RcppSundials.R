#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <ida/ida.h>           // CVODES functions and constants
#include <nvector/nvector_serial.h>  // Serial N_Vector
#include <ida/ida_dense.h>     // CVDense
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

// Interface between cvode integrator and the model function written in Cpp
int ida_to_Cpp_stl(double t, N_Vector y, N_Vector ydot, N_Vector rr, void* inputs) {
  // Cast the void pointer back to the correct data structure
  data_ida_Cpp_stl* data = static_cast<data_ida_Cpp_stl*>(inputs);
  // Interpolate the forcings
  vector<double> forcings(data->forcings_data.size());
  if(data->forcings_data.size() > 0) forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  vector<double> states(data->neq);
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i);
  // Extract the derivatives from the NV_Ith_S container
  vector<double> derivatives(data->neq);
  for(auto i = 0; i < data->neq ; i++) derivatives[i] = NV_Ith_S(ydot,i);  
  // Run the model
  array<vector<double>, 2> output = data->model(t, states, derivatives, data->parameters, forcings); 
  // Return the residuals to the NV_Ith_S
  vector<double> residuals = output[0];
  for(auto i = 0; i < data->neq; i++)  NV_Ith_S(rr,i) = residuals[i];
  return 0; 
}



//' Simulates the model when it is written as a C++ function using stl containers
//' @export
// [[Rcpp::export]]
NumericMatrix ida_Cpp_stl(NumericVector times, NumericVector states_, 
                          NumericVector derivatives_,
                        NumericVector parameters_, List forcings_data_, 
                        List settings, SEXP model_, SEXP jacobian_) {
  // Wrap the pointer to the model function with the correct signature                        
  dae_in_Cpp_stl* model =  (dae_in_Cpp_stl *) R_ExternalPtrAddr(model_);
  // Wrap the pointer to the jacobian function with the correct signature                        
  jac_in_Cpp_stl* jacobian =  nullptr;
  if(as<int>(settings["jacobian"]) == 1) 
      jacobian = (jac_in_Cpp_stl *) R_ExternalPtrAddr(jacobian_);
  // Store all inputs in the data struct, prior conversion to stl and Armadillo classes
  int neq = states_.size();
  vector<double> parameters{as<vector<double>>(parameters_)};
  vector<double> states{as<vector<double>>(states_)};
  vector<double> derivatives{as<vector<double>>(derivatives_)};
  vector<mat> forcings_data(forcings_data_.size());
  if(forcings_data_.size() > 0) 
    for(int i = 0; i < forcings_data_.size(); i++)
      forcings_data[i] = as<mat>(forcings_data_[i]);
  data_ida_Cpp_stl data_model{parameters, forcings_data, neq, model, jacobian};
  
  /*
   *
   Initialize IDA and pass all initial inputs
   *
   */
  N_Vector y = nullptr;
  y = N_VNew_Serial(neq);
  for(int i = 0; i < neq; i++) {
    NV_Ith_S(y,i) = states[i];
  }
  N_Vector yp = nullptr;
  yp = N_VNew_Serial(neq);
  for(int i = 0; i < neq; i++) {
    NV_Ith_S(yp,i) = derivatives[i];
  }  
  void *ida_mem = IDACreate();
  // Shut up Sundials (errors not printed to the screen)
  int flag = IDASetErrFile(ida_mem, NULL);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetErrFile function");  
  }
  
  // Initialize the Sundials solver. Here we pass initial N_Vector, the interface function and the initial time
  flag = IDAInit(ida_mem, ida_to_Cpp_stl, times[0], y, yp);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDAInit function"); 
  }
  
  // Tell Sundials the tolerance settings for error control
  flag = IDASStolerances(ida_mem, settings["rtol"], settings["atol"]);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);}     
    ::Rf_error("Error in the IDASStolerances function");
  }

  // Tell Sundials the number of state variables, so that I can allocate memory for the Jacobian
  flag = IDADense(ida_mem, neq);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDADense function");
  }
    
  // Give Sundials a pointer to the struct where all the user data is stored. It will be passed (untouched) to the interface as void pointer
  flag = IDASetUserData(ida_mem, &data_model);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetUserData function");
  }

  // Correct initial values
  flag = IDACalcIC(ida_mem, IDA_Y_INIT, times[1]);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDACalcIC function");  
  }

  // Set maximum number of steps
  flag = IDASetMaxNumSteps(ida_mem, settings["maxsteps"]);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetMaxNumSteps function");
  }
  
  // Set maximum order of the integration
  flag = IDASetMaxOrd(ida_mem, settings["maxord"]); 
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetMaxOrd function");
  }
  
  // Set the initial step size
  flag = IDASetInitStep(ida_mem, settings["hini"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetInitStep function");
  }
  
  // Set the maximum step size
  flag = IDASetMaxStep(ida_mem, settings["hmax"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetMaxStep function");  
  }
  
  // Set the maximum number of error test fails
  flag = IDASetMaxErrTestFails(ida_mem, settings["maxerr"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetMaxErrTestFails function");  
  }
  
  // Set the maximum number of nonlinear iterations per step
  flag = IDASetMaxNonlinIters(ida_mem, settings["maxnonlin"]);  
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetMaxNonlinIters function");  
  }
  
  // Set the maximum number of convergence failures
  flag = IDASetMaxConvFails(ida_mem, settings["maxconvfail"]);   
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("Error in the IDASetMaxConvFails function");
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
    first_call = model(times[0], states, derivatives, parameters, forcings); 
  } catch(std::exception &ex) {
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
      forward_exception_to_r(ex);
  } catch(...) { 
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
    ::Rf_error("c++ exception (unknown reason)"); 
  }

  // Check length of time derivatives against the information passed through settings
  if(first_call[0].size() != neq) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);}     
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
  
  mat output(times.size(), extract_states.size() + extract_observed.size() + 1);
  output(0,0) = times[0];
  //for(auto i = 0; i < neq; i++) {
  //    output(0,i+1) = states[i];
  for(auto it = extract_states.begin(); it != extract_states.end(); it++) {
    if(*it > NV_LENGTH_S(y)) {
      Rcout << "The index " << *it << " exceeds the number of state variables of the model" << '\n';
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
      ::Rf_error("Simulation exited because of error in extracting state variables");
    }
    output[0,*it] = NV_Ith_S(y,*it - 1);
  }     
  if(extract_observed.size()  > 0) {
      //for(auto i = 0; i < nder; i++)  
      //    output(0,i + 1 + neq) = observed[i];
      for(auto it = extract_observed.begin(); it != extract_observed.end(); it++) {
        if(*it > observed.size()) {
          Rcout << "The index " << *it << " exceeds the number of observed variables returned by the model" << '\n';
          if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
          if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
          ::Rf_error("Simulation exited because of error in extracting observed variables");
        }
        output[0,*it + extract_states.size()] = observed[*it - 1];
      }     
  }
  
  /*
   * 
   Main time loop. Each timestep call cvode. Handle exceptions and fill up output
   *
   */
  double t = times[0];
  for(unsigned i = 1; i < times.size(); i++) {
    try {
      flag = IDASolve(ida_mem, times[i], &t, y, yp, IDA_NORMAL);
      for(auto h = 0; h < neq; h++) {
        if(as<bool>(settings["positive"])) {
           if(NV_Ith_S(y,h) < as<double>(settings["minimum"])) {
             Rcout << "The state variable at positon " << h << " became smaller than minimum: " << NV_Ith_S(y,h) << '\n';
             if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
             if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);}  
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
      if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
      forward_exception_to_r(ex);
    } catch(...) { 
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
      ::Rf_error("c++ exception (unknown reason)"); 
    }

     // Write to the output matrix the new values of state variables and time
    output(i,0) = times[i];
    //for(auto h = 0; h < neq; h++) output(i,h + 1) = NV_Ith_S(y,h);
    for(auto it = extract_states.begin(); it != extract_states.end(); it++) {
      output[i,*it] = NV_Ith_S(y,*it - 1);
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
      std::array<vector<double>, 2> model_call  = model(times[i], states, derivatives, parameters, forcings); 
      observed =  model_call[1];
      // Derived variables already stored by the interface function
      //for(auto j = 0; j < nder; j++)  output(i,j + 1 + neq) = observed[j]; 
      for(auto it = extract_observed.begin(); it != extract_observed.end(); it++) {
        output[i,*it + extract_states.size()] = observed[*it - 1];
      }       
    } 
  }
              
  // De-allocate the N_Vector and the ida_mem structures
  if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
  if(ida_mem == nullptr) {free(ida_mem);} else {IDAFree(&ida_mem);} 
  
  return wrap(output);
}

//' Allows calling the model that calculates the time derivatives
//' @export
// [[Rcpp::export]]
List ida_calc_res(SEXP model_, NumericVector t, NumericVector states, 
                    NumericVector derivatives, NumericVector parameters, 
                    List forcings_data_) {
   // Wrap the pointer to the model function with the correct signature                        
  dae_in_Cpp_stl* model =  (dae_in_Cpp_stl *) R_ExternalPtrAddr(model_); 
  // Interpolate the forcings
  vector<mat> forcings_data(forcings_data_.size());
  if(forcings_data_.size() > 0) 
    for(int i = 0; i < forcings_data_.size(); i++)
      forcings_data[i] = as<mat>(forcings_data_[i]);
  vector<double> forcings(forcings_data.size());
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, t[0]);
  // Call the model
  array<vector<double>, 2> output = model(t[0], as<vector<double>>(states),
                                          as<vector<double>>(derivatives),
                                          as<vector<double>>(parameters),
                                          forcings);
  // return the output as a list
  return List::create(_["Residuals"] = wrap(output[0]),
                      _["Observed"] = wrap(output[1]));
}

//' Allows calling the function to calculate the Jacobian matrix of the model
//' @export
// [[Rcpp::export]]
NumericMatrix ida_calc_jac(SEXP jacobian_, NumericVector t, NumericVector states, 
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
