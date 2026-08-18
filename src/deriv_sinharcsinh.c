#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <string.h>
#include "deriv.h"


SEXP c_deriv_sinharcsinh(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp) {

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // Setup Input Parameters (pointers + strides) for all parameters of the distribution.
    InputVector x_vec     = make_input_var(x_sexp, NULL, N);
    InputVector mu_vec    = make_input_var(params_sexp, "mu",    N);
    InputVector sigma_vec = make_input_var(params_sexp, "sigma", N);
    InputVector nu_vec    = make_input_var(params_sexp, "nu",    N);
    InputVector tau_vec   = make_input_var(params_sexp, "tau",   N);

    // Setting up flags for score and hessian; used for fast bitwise operations
    // to see what needs to be calculated. 'req_scores' means 'required scores',
    // analogously 'req_hessian' is used to check what is required to get the
    // requested hessian elements.
    ParameterFlags req_scores  = get_flags_msnt(score_sexp, false);
    ParameterFlags req_hessian = get_flags_msnt(hessian_sexp, true);

    // Initialize empty 'score' double pointers
    double* s_mu    = NULL;
    double* s_sigma = NULL;
    double* s_tau   = NULL;
    double* s_nu    = NULL;

    // If score for mu, sigma, nu, or tau is requested: Initialize double
    // pointer for the corresponding `score_sexp` list element to be modified
    // and returned to R.
    if (req_scores & PARAM_MU)    s_mu    = getListElement(score_sexp, "mu");
    if (req_scores & PARAM_SIGMA) s_sigma = getListElement(score_sexp, "sigma");
    if (req_scores & PARAM_NU)    s_nu    = getListElement(score_sexp, "nu");
    if (req_scores & PARAM_TAU)   s_tau   = getListElement(score_sexp, "tau");



    // Must be private for OMP
    double x, mu, sigma, nu, tau;
    double z, asinhz, exp_tauasinhz, exp_minusnuasinhz, z2p1sqrtinv;
    double z2, nu2, tau2, sigmainv, r, c, h;
    double dldr, dldc, dldz, dcdz, drdz, dzdm, dldm;
    // Main loop
    for (int i = 0; i < N; i++) {
        // Extracting numeric elements
        x     = get_val(&x_vec, i);
        mu    = get_val(&mu_vec, i);
        sigma = get_val(&sigma_vec, i);
        nu    = get_val(&nu_vec, i);
        tau   = get_val(&tau_vec, i);

        z                 = (x - mu) / sigma;
        asinhz            = asinh(z);
        exp_tauasinhz     = exp( tau * asinhz);
        exp_minusnuasinhz = exp(-nu  * asinhz);

        // Performing a series of vector operations used multiple times below
        z2       = z * z;
        tau2     = tau * tau;
        nu2      = nu * nu;
        sigmainv = 1.0 / sigma;

        r = 0.5 * (exp_tauasinhz        - exp_minusnuasinhz);
        c = 0.5 * (exp_tauasinhz * tau  + exp_minusnuasinhz * nu);
        h = 0.5 * (exp_tauasinhz * tau2 - exp_minusnuasinhz * nu2);

        // Partial derivatives used everywhere
        dldr = -r;
        dldc = 1 / c;

        if (req_scores & PARAM_MU) {
            z2p1sqrtinv = 1.0 / sqrt(z2 + 1.0);

            dldz        = -z / (1.0 + z2);
            dcdz        = h * pow(1.0 + z2, -0.5);
            drdz        = c * pow(1.0 + z2, -0.5);
            dzdm        = -sigmainv;
            dldm        = sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv);

            s_mu[i] = (dldr * drdz + dldc * dcdz + dldz) * dzdm;
        }
    }


    if (req_scores & PARAM_MU) {
        Rprintf(" ---- score requires mu\n");
    }
    if (req_hessian & PARAM_MU) {
        Rprintf(" ---- hessian  requires mu\n");
    }
    if (req_hessian & PARAM_SIGMA) {
        Rprintf(" ---- hessian  requires sigma\n");
    }
    if (req_hessian & PARAM_MU_SIGMA) {
        Rprintf(" ---- hessian  requires mu:sigma\n");
    }
    if (req_hessian & PARAM_MU_MU) {
        Rprintf(" ---- hessian  requires mu:mu\n");
    }
    // R //   ## Calculating a series of vectors used over and over again
    // R //   ## when calculating the score(s).
    // R //   z                 <- (x - d$mu) / d$sigma
    // R //   asinhz            <- asinh(z)
    // R //   exp_tauasinhz     <- exp(d$tau * asinhz)
    // R //   exp_minusnuasinhz <- exp(-d$nu * asinhz)

    // R //   ## Performing a series of vector operations used multiple times below
    // R //   z2       <- z^2
    // R //   tau2     <- d$tau^2
    // R //   nu2      <- d$nu^2
    // R //   sigmainv <- 1 / d$sigma

    // R //   r <- 0.5 * (exp_tauasinhz          - exp_minusnuasinhz)
    // R //   c <- 0.5 * (exp_tauasinhz * d$tau  + exp_minusnuasinhz * d$nu)
    // R //   h <- 0.5 * (exp_tauasinhz * tau2   - exp_minusnuasinhz * nu2)

    // R //   ## Partial derivatives used everywhere
    // R //   dldr <- -r
    // R //   dldc <- 1 / c

    // R //   src_mu <-  function() {
    // R //     z2p1sqrtinv <- 1 / sqrt(z2 + 1)

    // R //     dldz <- -z / (1 + z2)
    // R //     dcdz <- h * (1 + z2)^(-0.5)
    // R //     drdz <- c * (1 + z2)^(-0.5)
    // R //     dzdm <- -sigmainv

    // R //     dldm <- sigmainv * z2p1sqrtinv * (-h / c + r * c + z * z2p1sqrtinv)
    // R //     return((dldr * drdz + dldc * dcdz + dldz) * dzdm)
    // R //   }
    // R //   src_sigma <- function() {
    // R //     z2p1sqrtinv <- 1 / sqrt(z2 + 1)

    // R //     dldz <- -z / (1 + z2)
    // R //     dcdz <- h * z2p1sqrtinv
    // R //     drdz <- c * z2p1sqrtinv
    // R //     dzdd <- -z * sigmainv
    // R //     return((dldr * drdz + dldc * dcdz + dldz) * dzdd - sigmainv)
    // R //   }
    // R //   src_nu <- function() {
    // R //       drdv <- 0.5 * asinhz * exp_minusnuasinhz
    // R //       dcdv <- 0.5 * (1 - d$nu * asinhz) * exp_minusnuasinhz
    // R //       return(dldr * drdv + dldc * dcdv)
    // R //   }
    // R //   src_tau <- function() {
    // R //       drdt <- 0.5 * asinhz * exp_tauasinhz
    // R //       dcdt <- 0.5 * (1 + d$tau * asinhz) * exp_tauasinhz
    // R //       return(dldr * drdt + dldc * dcdt)
    // R //   }

    return R_NilValue;

////    // Pre-calculate composite flags for fast branch pruning inside loop
////    bool calc_sigma = (req_scores & PARAM_SIGMA) || (req_hessian & PARAM_SIGMA);
////    bool calc_mu    = (req_scores & PARAM_MU)    || (req_hessian & PARAM_MU);
////    bool calc_nu    = (req_scores & PARAM_NU)    || (req_hessian & PARAM_NU);
}
