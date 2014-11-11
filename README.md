RcppSundials
============

This package consists of a wrapper around Sundials numerical library using Rcpp.

The necessary Sundials files are included such that the user does not need to have to install Sundials in the system.

The package is built in such a way that the system of ordinary differential equations may be implemented in R or in C++ using either the Rcpp data types or containers from the standard library and the Armadillo library.

Note that this package uses a different interface to that required by the deSolve package. The main different is that interpolation of external forcings is performed by RcppSundials whereas deSolve does not. Also, when writing the model in C++, inputs are not passed as global variables. The package is also designed in such a way that third packages can call directly the C++ functions in charge of performing the simulation from within C++, which means that computationally intensive tasks such as parameter optimization or sensitivity analyses can all be done in C++ without going through the R interpreter.
