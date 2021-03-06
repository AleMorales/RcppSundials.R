# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Returns the solution to the IVP problem at specified timepoints. 
#' Object returned is a matrix where nrow = number timepoints &
#' ncol = no. states + no. observed + [(no. states)x(no. parameters)]
#' the third term above is only returned when sensitivity equations are being
#' calculated
#' @export
wrap_cvodes <- function(times, states_, parameters_, forcings_data_, settings, model_, jacobian_) {
    .Call('RcppSundials_wrap_cvodes', PACKAGE = 'RcppSundials', times, states_, parameters_, forcings_data_, settings, model_, jacobian_)
}

#' Allows calling the model that calculates the time derivatives
#' @export
cvode_calc_derivs <- function(model_, t, states, parameters, forcings_data_) {
    .Call('RcppSundials_cvode_calc_derivs', PACKAGE = 'RcppSundials', model_, t, states, parameters, forcings_data_)
}

#' Allows calling the function to calculate the Jacobian matrix of the model
#' @export
cvode_calc_jac <- function(jacobian_, t, states, parameters, forcings_data_) {
    .Call('RcppSundials_cvode_calc_jac', PACKAGE = 'RcppSundials', jacobian_, t, states, parameters, forcings_data_)
}

#' Simulates the model when it is written as a C++ function using stl containers
#' @export
ida_Cpp_stl <- function(times, states_, derivatives_, parameters_, forcings_data_, settings, model_, jacobian_) {
    .Call('RcppSundials_ida_Cpp_stl', PACKAGE = 'RcppSundials', times, states_, derivatives_, parameters_, forcings_data_, settings, model_, jacobian_)
}

#' Allows calling the model that calculates the time derivatives
#' @export
ida_calc_res <- function(model_, t, states, derivatives, parameters, forcings_data_) {
    .Call('RcppSundials_ida_calc_res', PACKAGE = 'RcppSundials', model_, t, states, derivatives, parameters, forcings_data_)
}

#' Allows calling the function to calculate the Jacobian matrix of the model
#' @export
ida_calc_jac <- function(jacobian_, t, states, parameters, forcings_data_) {
    .Call('RcppSundials_ida_calc_jac', PACKAGE = 'RcppSundials', jacobian_, t, states, parameters, forcings_data_)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('RcppSundials_RcppExport_registerCCallable', PACKAGE = 'RcppSundials')
})
