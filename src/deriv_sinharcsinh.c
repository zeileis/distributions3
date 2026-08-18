#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <string.h>
#include "deriv.h"


SEXP c_deriv_sinharcsinh(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp) {

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // Setup Input Parameters (pointers + strides) for all parameters of the distribution.
    InputVector x_var     = make_input_var(x_sexp, NULL, N);
    InputVector mu_var    = make_input_var(params_sexp, "mu",    N);
    InputVector sigma_var = make_input_var(params_sexp, "sigma", N);
    InputVector nu_var    = make_input_var(params_sexp, "nu",    N);

    // Setting up flags for score and hessian; used for fast bitwise operations
    // to see what needs to be calculated.
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


    ///double* s_mu    = REAL(s_mu_sexp);
    ///double* s_sigma = REAL(s_sigma_sexp);

    ///SEXP s_sigma_sexp = getListElement(score_sexp, "sigma");
    ///double* s_sigma = REAL(s_sigma_sexp);

    for (int i = 0; i < N; i++) {
        Rprintf(" ---- i = %d\n", i);
        Rprintf("      mu.stride = %d\n", mu_var.stride);
        Rprintf("      mu[%d] = %.5f\n", i, get_val(&mu_var, i));
        Rprintf("      s_mu[%d] = %.5f\n", i, s_mu[i]);
        s_mu[i] = 100. / i;
        if (s_sigma != NULL) {
            s_sigma[i] = 500. / i;
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

    return R_NilValue;

////    // Pre-calculate composite flags for fast branch pruning inside loop
////    bool calc_sigma = (req_scores & PARAM_SIGMA) || (req_hessian & PARAM_SIGMA);
////    bool calc_mu    = (req_scores & PARAM_MU)    || (req_hessian & PARAM_MU);
////    bool calc_nu    = (req_scores & PARAM_NU)    || (req_hessian & PARAM_NU);
}
