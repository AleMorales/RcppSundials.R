library(RcppSundials)
library(microbenchmark)

example_model_R <- function(t, y, parms, forcs) {
  derivatives = -y*parms[1]
  return(list(derivatives))
}
simulationR = cvode_R(times = 0:1, 
                     states = rep(1,4), 
                     parameters = 0.1, 
                     forcings_data = NULL,
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = example_model_R) 

compiled_example_model = getNativeSymbolInfo(name = "RcppSundials_example_model",
                                             PACKAGE = "RcppSundials")$address

simulationCpp = cvode_Cpp(times = 0:1, 
                     states = rep(1.0,4), 
                     parameters = 0.1, 
                     forcings_data = NULL,
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                      model = compiled_example_model) 

microbenchmark(testR = cvode_R(times = 0:1, 
                     states = rep(1,4), 
                     parameters = 0.1, 
                     forcings_data = NULL,
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = example_model_R))