#include <Rcpp.h>
#include <cvode/cvode.h>           // CVODES functions and constants
#include <nvector/nvector_serial.h>  // Serial N_Vector
#include <cvode/cvode_dense.h>     // CVDense
#include <datatypes.h>
#include <string> 
#include <limits> 
#include <time.h>

using namespace Rcpp; 
// [[Rcpp::interfaces(r, cpp)]]

/*
 *
 Functions to interpolate the forcings
 *
 */
 
// Interpolate the forcings when passed as a matrix
double interpolate(NumericMatrix data, double xout) {
  NumericVector x = data(_,0);
  NumericVector y = data(_,1);
  // Retrieve the range of x values where xout is located
  xout = xout + sqrt(std::numeric_limits<double>::epsilon());
  // use rule = 2 for xout > x[end]  
  if(xout > *(x.end() - 1)) {
      return(*(y.end() - 1));
  // use rule = 2 for xout < x[0]  
  } else if (xout < x[0]) {
      return(y[0]);
 // Otherwise do linear interpolation
 } else {
   // This is supposed to be fast if the vector is ordered...
   auto it = std::lower_bound(x.begin(),x.end(),xout);
   int index = (it - x.begin()) - 1;
   double inty = (y[index + 1] - y[index])/(x[index + 1] - x[index])*(xout -x[index]) + y[index];
   return(inty);
 }  
}

// Interpolate the list of forcings and return the values as a NumericVector
NumericVector interpolate_list(List data, double t) {
  NumericVector forcings(data.size());
  for(auto i = 0; i < data.size(); i++) {
    forcings[i] = interpolate(data[i], t);
  }
  return forcings;
}
  
/*
 *
 Functions when the model is written in C++
 *
 */

// Interface between cvode integrator and the model function written in Cpp
int cvode_to_Cpp(double t, N_Vector y, N_Vector ydot, void* inputs) {
  // Cast the void pointer back to the correct data structure
  data_Cpp* data = (data_Cpp*) inputs;
  data = static_cast<data_Cpp*>(inputs);
  // Interpolate the forcings
  NumericVector forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  NumericVector states;
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i);
  // Run the model
  List output = data->model(wrap(t), states, data->parameters, forcings); 
  // Return the states to the NV_Ith_S
  NumericVector derivatives = output[0];
  for(auto i = 0; i < data->neq; i++)  NV_Ith_S(ydot,i) = derivatives[i];
  return 0; 
}

//' Simulates the model when it is written as a C++ function
//' @export
// [[Rcpp::export]]
NumericMatrix cvode_Cpp(NumericVector times, NumericVector states, 
                        NumericVector parameters, List forcings_data, 
                        List settings, SEXP model_) {
  // Wrap the pointer to the model function with the correct signature                        
  ode_in_Cpp* model =  (ode_in_Cpp *) R_ExternalPtrAddr(model_);
  // Store all inputs in the data struct
  int neq = states.size();
  data_Cpp data_model = {parameters, forcings_data, states.size(), model};
  
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
    stop("Please choose bdf or adams as method");
  }
  // Shut up Sundials (errors not printed to the screen)
  int flag = CVodeSetErrFile(cvode_mem, NULL);
  if(flag < 0.0) stop("Error in the CVodeSetErrFile function");  
  
  // Initialize the Sundials solver. Here we pass initial N_Vector, the interface function and the initial time
  flag = CVodeInit(cvode_mem, cvode_to_Cpp, times[0], y);
  if(flag < 0.0) stop("Error in the CVodeInit function"); 

  // Tell Sundials the tolerance settings for error control
  flag = CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"]);
  if(flag < 0.0) stop("Error in the CVodeSStolerances function");

  // Tell Sundials the number of state variables, so that I can allocate memory for the Jacobian
  flag = CVDense(cvode_mem, neq);
  if(flag < 0.0) stop("Error in the CVDense function");

  // Give Sundials a pointer to the struct where all the user data is stored. It will be passed (untouched) to the interface as void pointer
  flag = CVodeSetUserData(cvode_mem, &data_model);
  if(flag < 0.0) stop("Error in the CVodeSetUserData function");

  // Set maximum number of steps
  flag = CVodeSetMaxNumSteps(cvode_mem, settings["maxsteps"]);
  if(flag < 0.0) stop("Error in the CVodeSetUserData function");
  
  // Set maximum order of the integration
  flag = CVodeSetMaxOrd(cvode_mem, settings["maxord"]); 
  if(flag < 0.0) stop("Error in the CVodeSetMaxOrd function");
  
  // Set the initial step size
  flag = CVodeSetInitStep(cvode_mem, settings["hini"]);  
  if(flag < 0.0) stop("Error in the CVodeSetInitStep function");
  
  // Set the minimum step size
  flag = CVodeSetMinStep(cvode_mem, settings["hmin"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMinStep function");
  
  // Set the maximum step size
  flag = CVodeSetMaxStep(cvode_mem, settings["hmax"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMaxStep function");  
  
  // Set the maximum number of error test fails
  flag = CVodeSetMaxErrTestFails(cvode_mem, settings["maxerr"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMaxErrTestFails function");  
  
  // Set the maximum number of nonlinear iterations per step
  flag = CVodeSetMaxNonlinIters(cvode_mem, settings["maxnonlin"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMaxNonlinIters function");  
  
  // Set the maximum number of convergence failures
  flag = CVodeSetMaxConvFails(cvode_mem, settings["maxconvfail"]);   
  if(flag < 0.0) stop("Error in the CVodeSetMaxConvFails function");
  
  /*
   * 
   Make a first call to the model to check that everything is ok and retrieve the number of observed variables
   *
   */
  NumericVector forcings;
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[0]);
  List first_call;
  try {
    first_call = model(wrap(times[0]), states, parameters, forcings); 
  } catch(std::exception &ex) {
      forward_exception_to_r(ex);
  } catch(...) { 
    ::Rf_error("c++ exception (unknown reason)"); 
  }
  
  /*
   * 
   Fill up the output matrix with the values for the initial time
   *
   */
  NumericVector observed;
  int nder = 0;
  if(first_call.size() == 2) {
    observed = first_call[1];
    nder = observed.size();
  }
  NumericMatrix output(times.size(), neq + nder + 1);
  output(0,0) = times[0];
  for(auto i = 0; i < neq; i++) 
      output(0,i+1) = states[i];
  if(nder  > 0) 
      for(auto i = 0; i < nder; i++)  
          output(0,i + 1 + neq) = observed[i];
  // Initialize counter to facilitate interruption of long simulations
  clock_t t_initial, t_varying; // initial and varying times
  t_initial = clock(); 
  double elapsed_time;
  
  
  /*
   * 
   Main time loop. Each timestep call cvode. Handle exceptions and fill up output
   *
   */
  double t = times[0];
  for(unsigned i = 1; i < times.size(); i++) {
    try {
      flag = CVode(cvode_mem, times[i], y, &t, CV_NORMAL);
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
    // Calculate the time that has elapsed since beginning simulation (in seconds)
    t_varying = clock();
    elapsed_time = ((float)t_varying - (float)t_initial)/CLOCKS_PER_SEC;
    // If the time exceeds maximum, throw an exception and finish the simulation
    if(elapsed_time >= as<double>(settings["maxtime"]))   
      throw std::runtime_error("The simulation time exceeded the maximum.");
    } catch(std::exception &ex) {
      forward_exception_to_r(ex);
    } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
    }

     // Write to the output matrix the new values of state variables and time
    output(i,0) = times[i];
    for(auto h = 0; h < neq; h++) output(i,h + 1) = NV_Ith_S(y,h);
  }
  
  // If we have observed variables we call the model function again
  if(nder > 0 && flag >= 0.0) {
    for(unsigned i = 1; i < times.size(); i++) {
      // Get forcings values at time 0.
      NumericVector forcings;
      if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[i]);
      // Get the state variables into the_data_states
      NumericVector simulated_states(neq);
      for(auto j = 0; j < neq; j++) simulated_states[j] = output(i,j + 1);
      // Call the model function to retrieve total number of outputs and initial values for derived variables
      List model_call  = model(wrap(times[i]), simulated_states, parameters, forcings); 
      observed  = model_call[1];
      // Derived variables already stored by the interface function
      for(auto j = 0; j < nder; j++)  output(i,j + 1 + neq) = observed[j]; 
    } 
  }
              
  // De-allocate the N_Vector and the cvode_mem structures
  if(y == nullptr) {
    free(y);
  } else {
    N_VDestroy_Serial(y);
  }
  if(cvode_mem == nullptr) {
    free(cvode_mem);
  } else {
    CVodeFree(&cvode_mem);
  }  
  
  return output;
}

/*
 *
 Functions when the model is written in R
 *
 */
 
// Interface between cvode integrator and the model function written in R
int cvode_to_R(double t, N_Vector y, N_Vector ydot, void* inputs) {
  data_R* data = (data_R*) inputs;
  data = static_cast<data_R*>(inputs);
  // Interpolate the forcings
  NumericVector forcings = interpolate_list(data->forcings_data, t);
  // Extract the states from the NV_Ith_S container
  NumericVector states;
  for(auto i = 0; i < data->neq ; i++) states[i] = NV_Ith_S(y,i);
  // Run the model
  List output = data->model(wrap(t), states, data->parameters, forcings); 
  // Return the states to the NV_Ith_S
  NumericVector derivatives = output[0];
  for(auto i = 0; i < data->neq; i++)  NV_Ith_S(ydot,i) = derivatives[i];
  return 0; 
}

//' Simulates the model when it is written as a R function
//' @export
// [[Rcpp::export]]
NumericMatrix cvode_R(NumericVector times, NumericVector states, 
                        NumericVector parameters, List forcings_data, 
                        List settings, Function model) {
  // Store all inputs in the data struct
  int neq = states.size();
  data_R data_model = {parameters, forcings_data, states.size(), model};
  
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
    stop("Please choose bdf or adams as method");
  }
  // Shut up Sundials (errors not printed to the screen)
  int flag = CVodeSetErrFile(cvode_mem, NULL);
  if(flag < 0.0) stop("Error in the CVodeSetErrFile function");  
  
  // Initialize the Sundials solver. Here we pass initial N_Vector, the interface function and the initial time
  flag = CVodeInit(cvode_mem, cvode_to_R, times[0], y);
  if(flag < 0.0) stop("Error in the CVodeInit function"); 

  // Tell Sundials the tolerance settings for error control
  flag = CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"]);
  if(flag < 0.0) stop("Error in the CVodeSStolerances function");

  // Tell Sundials the number of state variables, so that I can allocate memory for the Jacobian
  flag = CVDense(cvode_mem, neq);
  if(flag < 0.0) stop("Error in the CVDense function");

  // Give Sundials a pointer to the struct where all the user data is stored. It will be passed (untouched) to the interface as void pointer
  flag = CVodeSetUserData(cvode_mem, &data_model);
  if(flag < 0.0) stop("Error in the CVodeSetUserData function");

  // Set maximum number of steps
  flag = CVodeSetMaxNumSteps(cvode_mem, settings["maxsteps"]);
  if(flag < 0.0) stop("Error in the CVodeSetUserData function");
  
  // Set maximum order of the integration
  flag = CVodeSetMaxOrd(cvode_mem, settings["maxord"]); 
  if(flag < 0.0) stop("Error in the CVodeSetMaxOrd function");
  
  // Set the initial step size
  flag = CVodeSetInitStep(cvode_mem, settings["hini"]);  
  if(flag < 0.0) stop("Error in the CVodeSetInitStep function");
  
  // Set the minimum step size
  flag = CVodeSetMinStep(cvode_mem, settings["hmin"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMinStep function");
  
  // Set the maximum step size
  flag = CVodeSetMaxStep(cvode_mem, settings["hmax"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMaxStep function");  
  
  // Set the maximum number of error test fails
  flag = CVodeSetMaxErrTestFails(cvode_mem, settings["maxerr"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMaxErrTestFails function");  
  
  // Set the maximum number of nonlinear iterations per step
  flag = CVodeSetMaxNonlinIters(cvode_mem, settings["maxnonlin"]);  
  if(flag < 0.0) stop("Error in the CVodeSetMaxNonlinIters function");  
  
  // Set the maximum number of convergence failures
  flag = CVodeSetMaxConvFails(cvode_mem, settings["maxconvfail"]);   
  if(flag < 0.0) stop("Error in the CVodeSetMaxConvFails function");
  
  /*
   * 
   Make a first call to the model to check that everything is ok and retrieve the number of observed variables
   *
   */
  NumericVector forcings(0);
  if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[0]);
  List first_call(0);
  try {
    first_call  = model(times[0], states, parameters, forcings); 
  } catch(std::exception &ex) {
      forward_exception_to_r(ex);
  } 
  /*
   * 
   Fill up the output matrix with the values for the initial time
   *
   */
  
  NumericVector observed;
  int nder = 0;
  if(first_call.size() == 2) {
    observed = first_call[1];
    nder = observed.size();
  }
  NumericMatrix output(times.size(), neq + nder + 1);
  output(0,0) = times[0];
  for(auto i = 0; i < neq; i++) 
      output(0,i+1) = states[i];
  if(nder  > 0) 
      for(auto i = 0; i < nder; i++)  
          output(0,i + 1 + neq) = observed[i];
  // Initialize counter to facilitate interruption of long simulations
  clock_t t_initial, t_varying; // initial and varying times
  t_initial = clock(); 
  double elapsed_time;
  
  
  /*
   * 
   Main time loop. Each timestep call cvode. Handle exceptions and fill up output
   *
   */
  double t = times[0];
  for(unsigned i = 1; i < times.size(); i++) {
    try {
      flag = CVode(cvode_mem, times[i], y, &t, CV_NORMAL);
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
    // Calculate the time that has elapsed since beginning simulation (in seconds)
    t_varying = clock();
    elapsed_time = ((float)t_varying - (float)t_initial)/CLOCKS_PER_SEC;
    // If the time exceeds maximum, throw an exception and finish the simulation
    if(elapsed_time >= as<double>(settings["maxtime"]))   
      throw std::runtime_error("The simulation time exceeded the maximum.");
    } catch(std::exception &ex) {
      forward_exception_to_r(ex);
    } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
    }

     // Write to the output matrix the new values of state variables and time
    output(i,0) = times[i];
    for(auto h = 0; h < neq; h++) output(i,h + 1) = NV_Ith_S(y,h);
  }

  // If we have observed variables we call the model function again
  if(nder > 0 && flag >= 0.0) {
    for(unsigned i = 1; i < times.size(); i++) {
      // Get forcings values at time 0.
      NumericVector forcings;
      if(forcings_data.size() > 0) forcings = interpolate_list(forcings_data, times[i]);
      // Get the state variables into the_data_states
      NumericVector simulated_states(neq);
      for(auto j = 0; j < neq; j++) simulated_states[j] = output(i,j + 1);
      // Call the model function to retrieve total number of outputs and initial values for derived variables
      List model_call  = model(wrap(times[i]), simulated_states, parameters, forcings); 
      observed  = model_call[1];
      // Derived variables already stored by the interface function
      for(auto j = 0; j < nder; j++)  output(i,j + 1 + neq) = observed[j]; 
    } 
  }
               
  // De-allocate the N_Vector and the cvode_mem structures
  if(y == nullptr) {
    free(y);
  } else {
    N_VDestroy_Serial(y);
  }
  if(cvode_mem == nullptr) {
    free(cvode_mem);
  } else {
    CVodeFree(&cvode_mem);
  }  

  return output;  
}


