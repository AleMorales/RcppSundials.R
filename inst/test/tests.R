library(RcppSundials)
library(microbenchmark)

# Use the Cpp version that is wrapped in R
simulationR = cvode_R(times = c(0.0, 1.0), 
                     states = rep(1,4), 
                     parameters = 0.1, 
                     forcings_data = vector("list",0),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = example_model)

compiled_example_model = getNativeSymbolInfo(name = "RcppSundials_example_model",
                                             PACKAGE = "RcppSundials")$address
# Use the native Cpp version 
simulationCpp = cvode_Cpp(times = 0:1, 
                     states = rep(1.0,4), 
                     parameters = 0.1, 
                     forcings_data = NULL,
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                      model = compiled_example_model) 

microbenchmark(testR = cvode_R(times = 1:3600, 
                     states = rep(1,150), 
                     parameters = 0.1, 
                     forcings_data = NULL,
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = example_model),
               testCpp = cvode_Cpp(times = 1:3600, 
                     states = rep(1,150), 
                     parameters = 0.1, 
                     forcings_data = NULL,
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = compiled_example_model))