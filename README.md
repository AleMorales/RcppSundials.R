RcppSundials
============

This package consists of a wrapper around Sundials numerical library using Rcpp.

The necessary Sundials files are included such that the user does not need to have to install Sundials in the system.

The package is built in such a way that the system of ordinary differential equations must be implemented in in C++ using containers from the standard library and the Armadillo library.

Note that this package uses a different interface to that required by the deSolve package. The main different is that interpolation of external forcings is performed by RcppSundials whereas deSolve does not. Also, inputs cannot not passed as global variables (as in the C and Fortran models written for deSolve). The package is also designed in such a way that third packages can call directly the C++ functions in charge of performing the simulation without going through R, which means that computationally intensive tasks such as parameter optimization or sensitivity analyses may be performed faster than otherwise.
