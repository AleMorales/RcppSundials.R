library(RcppSundials)
library(microbenchmark)

example_model_R = function(t, states, parameters, forcings) {
   derivatives = -states*parameters[1];
   observed = forcings[1];
   return(list(derivatives, observed));
}

example_jacobian_R = function(t, states, parameters, forcings) {
   output = matrix(0, length(states), length(states));
   diag(output) = -parameters[1];
   return(output);
}

# Use the Cpp version that is wrapped in R
simulationR = cvode_R(times = 1:5, 
                     states = rep(1,4), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1, jacobian = 0), 
                     model = example_model_R, jacobian = example_jacobian_R)

stl_example_model = getNativeSymbolInfo(name = "example_model_stl",
                                             PACKAGE = "RcppSundials")$address
stl_example_jacobian = getNativeSymbolInfo(name = "example_jacobian_stl",
                                             PACKAGE = "RcppSundials")$address

simulationstl = cvode_Cpp_stl(times = 1:5, 
                     states = rep(1.0,4), 
                     parameters = 0.1, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1, jacobian = 0), 
                      model = stl_example_model, jacobian = "wat?") 

stl_example_dae = getNativeSymbolInfo(name = "example_dae_stl",
                                             PACKAGE = "RcppSundials")$address

simulationstl = ida_Cpp_stl(times = c(0,0.4,4,40,400), 
                     states = c(1,0,0), 
                     derivatives = c(-0.04,0.04,0),
                     parameters = 0.1, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-4, 
                      atol = 1e-6, maxsteps = 1e4, maxord = 5, hini = 0, 
                      hmin = 0, hmax = 0, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 0, jacobian = 0), 
                      model = stl_example_dae, jacobian = "wat?") 


system.time(testR <- cvode_R(times = 1:3600, 
                     states = rep(1,150), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1, jacobian = 1), 
                     model = example_model_R, jacobian = example_jacobian_R))

             
system.time(testStl <- cvode_Cpp_stl(times = 1:3600, 
                     states = rep(1,150), 
                     parameters = 0.1, 
                     forcings_data = list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3, 
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", maxtime = 1, jacobian = 0), 
                     model = stl_example_model, jacobian = stl_example_jacobian))