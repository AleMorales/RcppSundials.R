  library(RcppSundials)
  library(microbenchmark)

stl_example_model = getNativeSymbolInfo(name = "example_model_stl",
                                             PACKAGE = "RcppSundials")$address
stl_example_jacobian = getNativeSymbolInfo(name = "example_jacobian_stl",
                                             PACKAGE = "RcppSundials")$address

simulationstl = wrap_cvodes(times = 1:5, 
                     states_ = rep(1.0,5), 
                     parameters_ = 0.1, 
                     forcings_data_ =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6,
                      atol = rep(1e-6,5), maxsteps = 1e3, maxord = 5, hini = 0,
                      hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", jacobian = 0,
                      minimum = -1e-4, positive = 1,
                      which_states = 1:5, which_observed = 1,
                      stability = T), 
                      model_ = stl_example_model, jacobian_ = stl_example_jacobian) 

cvode_calc_derivs(stl_example_model, t = 0,  states = rep(1.0,5), parameters = 0.1, 
            forcings_data =list(cbind(1:3600,1:3600)))
jac = cvode_calc_jac(stl_example_jacobian, t = 0,  states = rep(1.0,5), parameters = 0.1, 
            forcings_data =list(cbind(1:3600,1:3600)))

stl_example_dae = getNativeSymbolInfo(name = "example_dae_stl",
                                             PACKAGE = "RcppSundials")$address

ida_calc_res(stl_example_dae, t = 0,  states = c(1,0,0), derivatives = c(-0.04,0.04,0),
                parameters = 0.1, forcings_data =list(cbind(1:3600,1:3600)))


simulationstl = ida_Cpp_stl(times = c(0,0.4,4,40,400), 
                     states = c(1,0,0), 
                     derivatives = c(-0.04,0.04,0),
                     parameters = 0.1, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-10, 
                      atol = 1e-10, maxsteps = 1e4, maxord = 5, hini = 0, 
                      hmin = 0, hmax = 0, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", jacobian = 0,
                      minimum = -1e-4, positive = 1, which_states = 1:3, which_observed = 1), 
                      model = stl_example_dae, jacobian = "wat?") 


# Check for segfaults associated to memory leaks
for(i in 1:100) {
  test = cvode_Cpp_stl(times = 1:1e3, 
                     states = rep(1.0,1e2), 
                     parameters = 0.0001, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e4, maxord = 5, hini = 0, 
                      hmin = 0, hmax = 0, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", jacobian = 0,
                      minimum = -1e-4, positive = 0,which_states = 1:5,which_observed = integer()), 
                      model = stl_example_model, jacobian = stl_example_jacobian)
  cvode_calc_derivs(stl_example_model, t = 0,  states = rep(1.0,5), parameters = 0.1, 
            forcings_data =list(cbind(1:3600,1:3600)))
  jac = cvode_calc_jac(stl_example_jacobian, t = 0,  states = rep(1.0,5), parameters = 0.1, 
            forcings_data =list(cbind(1:3600,1:3600)))
  
ida_calc_res(stl_example_dae, t = 0,  states = c(1,0,0), derivatives = c(-0.04,0.04,0),
                parameters = 0.1, forcings_data =list(cbind(1:3600,1:3600)))

simulationstl = ida_Cpp_stl(times = c(0,0.4,4,40,400), 
                     states = c(1,0,0), 
                     derivatives = c(-0.04,0.04,0),
                     parameters = 0.1, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-10, 
                      atol = 1e-10, maxsteps = 1e4, maxord = 5, hini = 0, 
                      hmin = 0, hmax = 0, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", jacobian = 0,
                      minimum = -1e-4, positive = 1, which_states = 1:3, which_observed = 1), 
                      model = stl_example_dae, jacobian = "wat?")   
}


             
microbenchmark::microbenchmark(test = cvode_Cpp_stl(times = 1:1e3, 
                     states = rep(1000.0,1e2), 
                     parameters = 0.001, 
                     forcings_data =list(cbind(1:3600,1:3600)),
                     settings = list(rtol = 1e-6, 
                      atol = 1e-6, maxsteps = 1e4, maxord = 5, hini = 0, 
                      hmin = 0, hmax = 0, maxerr = 5, maxnonlin = 10,
                      maxconvfail = 10, method = "bdf", jacobian = 1,
                      minimum = -1e-4, positive = 0,which_states = 1:1e2,which_observed = integer(1)), 
                      model = stl_example_model, jacobian = stl_example_jacobian), times = 10)