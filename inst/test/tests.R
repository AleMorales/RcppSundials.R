library(RcppSundials)
library(microbenchmark)

example_model_R = function(t, states, parameters, forcings) {
   derivatives = -states*parameters[1];
   observed = forcings[1];
   return(list(derivatives, observed));
}

# Use the Cpp version that is wrapped in R
simulationR = cvode_R(times = 1:5, 
                     states = rep(1,4), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = example_model_R)

compiled_example_model = getNativeSymbolInfo(name = "RcppSundials_example_model",
                                             PACKAGE = "RcppSundials")$address
# Use the native Cpp version 
simulationCpp = cvode_Cpp(times = 1:5, 
                     states = rep(1.0,4), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                      model = compiled_example_model) 

stl_example_model = getNativeSymbolInfo(name = "example_model_stl",
                                             PACKAGE = "RcppSundials")$address
simulationstl = cvode_Cpp_stl(times = 1:5, 
                     states = rep(1.0,4), 
                     parameters = 0.1, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                      model = stl_example_model) 


system.time(testR <- cvode_R(times = 1:36000, 
                     states = rep(1,1500), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = example_model))
            
system.time({for(i in 1:50) testCpp = cvode_Cpp(times = 1:3600, 
                     states = rep(1,150), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = compiled_example_model)})
             
             
             
system.time(testStl <- cvode_Cpp_stl(times = 1:36000, 
                     states = rep(1,1500), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1), 
                     model = stl_example_model))