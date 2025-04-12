//! This file contains wrappers for the GSL library.

const gsl_sf_fact  = @cImport(@cInclude("gsl/gsl_sf_fact.h" ));
const gsl_sf_gamma = @cImport(@cInclude("gsl/gsl_sf_gamma.h"));

/// Double factorial function on unsigned integers.
pub fn dfact(n: u32) f64 {
    return gsl_sf_gamma.gsl_sf_doublefact(n);
}

/// Gamma function.
pub fn gamma(x: f64) f64 {
    return gsl_sf_gamma.gsl_sf_gamma(x);
}

/// Regularized lower incomplete gamma function.
pub fn gammainc(a: f64, x: f64) f64 {
    return gsl_sf_gamma.gsl_sf_gamma_inc_P(a, x);
}
