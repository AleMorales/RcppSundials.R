#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cvode/cvode.h>           // CVODES functions and constants
#include <nvector/nvector_serial.h>  // Serial N_Vector
#include <cvode/cvode_dense.h>     // CVDense
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
 Functions when the model is written in R
 *
 */
 
// Interface between cvode integrator and the model function written in R
int cvode_to_R(double t, N_Vector y, N_Vector ydot, void* inputs) {
  data_R* data = static_cast<data_R*>(inputs);
  // Interpolate the forcings
  NumericVector forcings(data->forcings_data.size());
  if(data->forcings_data.size() > 0) forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  NumericVector states(data->neq);
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i);
  // Run the model
  List output = data->model(wrap(t), states, data->parameters, forcings); 
  // Return the states to the NV_Ith_S
  NumericVector derivatives = output[0];
  for(auto i = 0; i < data->neq; i++)  NV_Ith_S(ydot,i) = derivatives[i];
  return 0; 
}

// Interface between cvode integrator and the jacobian function written in R
int cvode_to_R_jac(long int N, double t, N_Vector y, N_Vector fy, DlsMat Jac, 
                   void *inputs, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  data_R* data = static_cast<data_R*>(inputs);
  // Interpolate the forcings
  NumericVector forcings(data->forcings_data.size());
  if(data->forcings_data.size() > 0) forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  NumericVector states(data->neq);
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i); 
  // Get the Jacobian
  NumericMatrix output = data->jacobian(wrap(t), states, data->parameters, forcings);  
  // Return the DlsMat
  for(int j = 0; j < output.ncol(); j++) {
    for(int i= 0; i < output.nrow(); i++) {
      DENSE_ELEM(Jac, i, j) = output(i,j);
    }
  }
  return 0;
}

//' Simulates the model when it is written as a R function
//' @export
// [[Rcpp::export]]
NumericMatrix cvode_R(NumericVector times, NumericVector states, 
                        NumericVector parameters, List forcings_data, 
                        List settings, Function model, Function jacobian) {
  // Store all inputs in the data struct
  int neq = states.size();
  data_R data_model{parameters, forcings_data, neq, model, jacobian};
  
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
  
  // Shut up Sundials (errors not printed to the screen)
  int flag = CVodeSetErrFile(cvode_mem, NULL);
  if(flag < 0.0) {
   if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
   if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}       
   ::Rf_error("Error in the CVodeSetErrFile function");  
  }
  
  // Initialize the Sundials solver. Here we pass initial N_Vector, the interface function and the initial time
  flag = CVodeInit(cvode_mem, cvode_to_R, times[0], y);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
    ::Rf_error("Error in the CVodeInit function"); 
  }
  
  // Tell Sundials the tolerance settings for error control
  flag = CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"]);
  if(flag < 0.0) {
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 
    ::Rf_error("Error in the CVodeSStolerances function");
  }

  // Tell Sundials the number of state variables, so that I can allocate memory for the Jacobian
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
    flag = CVDlsSetDenseJacFn(cvode_mem, cvode_to_R_jac);
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
  
  /*
   * 
   Make a first call to the model to check that everything is ok and retrieve the number of observed variables
   *
   */
  NumericVector forcings(forcings_data.size());
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[0]);
  List first_call(0);
  try {
    first_call  = model(times[0], states, parameters, forcings); 
  } catch(std::exception &ex) {
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
      forward_exception_to_r(ex);
  } catch(...) { 
    if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
    if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}     
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  /*
   * 
   Fill up the output matrix with the values for the initial time
   *
   */
  
  NumericVector observed(0);
  int nder = 0;
  if(first_call.size() == 2) {
    NumericVector temp =  first_call[1];
    for(int i = 0; i < temp.size(); i++) {
      observed.push_back(temp[i]);
    }
    nder = observed.size();
  }
  
  // Indices that we actually want to extract
  vector<int> extract_states = as<vector<int>>(settings["which_states"]);
  vector<int> extract_observed = as<vector<int>>(settings["which_observed"]);
  
  NumericMatrix output(times.size(), extract_states.size() + extract_observed.size() + 1);
  output(0,0) = times[0];
  //for(auto i = 0; i < neq; i++) 
  //    output(0,i+1) = states[i];
  for(auto it = extract_states.begin(); it != extract_states.end(); it++) {
    if(*it > NV_LENGTH_S(y)) {
      Rcout << "The index " << *it << " exceeds the number of state variables of the model" << '\n';
      if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
      if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}        
      ::Rf_error("Simulation exited because of error in extracting state variables");
    }
    output(0,*it) = NV_Ith_S(y,*it - 1);
  }  
  if(extract_observed.size() > 0) {
    //for(auto i = 0; i < nder; i++)  
    //    output(0,i + 1 + neq) = observed[i];
    for(auto it = extract_observed.begin(); it != extract_observed.end(); it++) {
      if(*it > observed.size()) {
        Rcout << "The index " << *it << " exceeds the number of observed variables returned by the model" << '\n';
        if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
        if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);}  
        ::Rf_error("Simulation exited because of error in extracting observed variables");
      }
      output(0,*it + extract_states.size()) = observed[*it - 1];
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
      flag = CVode(cvode_mem, times[i], y, &t, CV_NORMAL);
      for(auto h = 0; h < neq; h++) {
       if(as<bool>(settings["positive"])) {
         if(NV_Ith_S(y,h) < as<double>(settings["minimum"])) {
           Rcout << "The state variable at positon " << h << " became smaller than minimum: " << NV_Ith_S(y,h) << '\n';
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
          case -28:
            throw std::runtime_error("The output became negative"); break;  
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
    output(i,0) = times[i];
    //for(auto h = 0; h < neq; h++) output(i,h + 1) = NV_Ith_S(y,h);
    for(auto it = extract_states.begin(); it != extract_states.end(); it++) {
      output(i,*it) = NV_Ith_S(y,*it - 1);
    }
  }

  // If we have observed variables we call the model function again
  if(extract_observed.size() > 0 && flag >= 0.0) {
    for(unsigned i = 1; i < times.size(); i++) {
      // Get forcings values at time 0.
      NumericVector forcings(forcings_data.size());
      if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[i]);
      // Get the state variables into the_data_states
      NumericVector simulated_states(neq);
      for(auto j = 0; j < neq; j++) simulated_states[j] = output(i,j + 1);
      // Call the model function to retrieve total number of outputs and initial values for derived variables
      List model_call  = model(wrap(times[i]), simulated_states, parameters, forcings); 
      observed  = model_call[1];
      // Derived variables already stored by the interface function
      //for(auto j = 0; j < nder; j++)  output(i,j + 1 + neq) = observed[j]; 
      for(auto it = extract_observed.begin(); it != extract_observed.end(); it++) {
        output(i,*it + extract_states.size()) = observed[*it - 1];
      }
    } 
  }
               
  // De-allocate the N_Vector and the cvode_mem structures
  if(y == nullptr) {free(y);} else {N_VDestroy_Serial(y);}
  if(cvode_mem == nullptr) {free(cvode_mem);} else {CVodeFree(&cvode_mem);} 

  return output;  
}
