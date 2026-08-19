#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <string.h>
#include "deriv.h"

/* Bitwise Flags for Required Parameters
 *
 * Used to store which elements of a score/hessian are required.
 * Bitwise flags which allow for |/& operations. Note that
 * this 'only' allows a 32-bit shift, i.3., `1U << 31` is the
 * absolute maximum. If Other parameters are required write
 * a different type definition (different enum type).
 *
 * In other words, this only allows for up to a maximum of
 * 7 parameters as (N^2 - N) / 2 + N) with N = 7 results in 28.
 */
typedef enum {
    PARAM_NONE  = 0,

    // First derivative/main parameters (used for score and hessian)
    PARAM_MU    = 1U << 0,
    PARAM_SIGMA = 1U << 1,

    // Cross-derivatives (used for hessian)
    PARAM_MU_MU         = 1U << 2,
    PARAM_MU_SIGMA      = 1U << 3,
    PARAM_SIGMA_SIGMA   = 1U << 4,
    // .. up to 1U << 31
} ParameterFlags;

/* Setting Required Parameter Flags
 *
 * @param x a named list where the names define which score or hessian
 *        is required.
 * @param hessian set false for score, and true for hessian.
 *
 * Returns an enum object used to check which derivatives are needed.
 */
static inline ParameterFlags get_flags(SEXP x, bool hessian) {
    ParameterFlags req = PARAM_NONE;

    // Score
    if (!hessian) {
        if (getListElementSEXP(x, "mu")    != R_NilValue) req |= PARAM_MU;
        if (getListElementSEXP(x, "sigma") != R_NilValue) req |= PARAM_SIGMA;
    // Hessian
    } else {
        if (getListElementSEXP(x, "mu:mu")       != R_NilValue) req |= (PARAM_MU | PARAM_MU_MU);
        if (getListElementSEXP(x, "mu:sigma")    != R_NilValue) req |= (PARAM_MU | PARAM_SIGMA | PARAM_MU_SIGMA);
        if (getListElementSEXP(x, "sigma:sigma") != R_NilValue) req |= (PARAM_SIGMA | PARAM_SIGMA_SIGMA);
    }

    return req;
}


SEXP c_deriv_Normal(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP expected_sexp, SEXP ncores) {

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // Expected hessian?
    bool expected = LOGICAL(expected_sexp)[0];

    // Setup Input Parameters (pointers + strides) for all parameters of the distribution.
    InputVector x_vec     = make_input_var(x_sexp, NULL, N);
    InputVector mu_vec    = make_input_var(params_sexp, "mu",    N);
    InputVector sigma_vec = make_input_var(params_sexp, "sigma", N);

    // Setting up flags for score and hessian; used for fast bitwise operations
    // to see what needs to be calculated. 'req_scores' means 'required scores',
    // analogously 'req_hessian' is used to check what is required to get the
    // requested hessian elements.
    ParameterFlags req_scores  = get_flags(score_sexp, false);
    ParameterFlags req_hessian = get_flags(hessian_sexp, true);

    // If score for mu, sigma, nu, or tau is requested: Initialize double
    // pointer for the corresponding `score_sexp` list element to be modified
    // and returned to R.
    double* s_mu        = NULL;
    double* s_sigma     = NULL;
    if (req_scores & PARAM_MU)           s_mu    = getListElement(score_sexp, "mu");
    if (req_scores & PARAM_SIGMA)        s_sigma = getListElement(score_sexp, "sigma");

    // As for score: Initialize double pointer for required elements corresponding
    // to the elements in hessian_sexp to be modified and returned to R.
    double* h_mu_mu     = NULL;
    double* h_mu_sigma  = NULL;
    if (req_hessian & PARAM_MU_MU)       h_mu_mu    = getListElement(hessian_sexp, "mu:mu");
    if (req_hessian & PARAM_MU_SIGMA)    h_mu_sigma = getListElement(hessian_sexp, "mu:sigma");

    double* h_sigma_sigma = NULL;
    if (req_hessian & PARAM_SIGMA_SIGMA) h_sigma_sigma    = getListElement(hessian_sexp, "sigma:sigma");

    // Main loop
    #if _OPENMP
    //Rprintf("number of threads: %d\n", nthreads);
    #pragma omp parallel for num_threads(nthreads) if (N > 50000)
    #endif
    for (int i = 0; i < N; i++) {
        // Extracting numeric elements
        double x     = get_val(&x_vec, i);
        double mu    = get_val(&mu_vec, i);
        double sigma = get_val(&sigma_vec, i);

        if (ISNA(x) || ISNA(mu) || ISNA(sigma)) {
            if (req_scores & PARAM_MU)           s_mu[i] = NA_REAL;
            if (req_scores & PARAM_SIGMA)        s_sigma[i] = NA_REAL;
            if (req_hessian & PARAM_MU_MU)       h_mu_mu[i] = NA_REAL;
            if (req_hessian & PARAM_MU_SIGMA)    h_mu_sigma[i] = NA_REAL;
            if (req_hessian & PARAM_SIGMA_SIGMA) h_sigma_sigma[i] = NA_REAL;
            continue;
        }

        double xmmu   = x - mu;

        double sigmainv  = 1.0 / sigma;
        double sigma2inv = sigmainv * sigmainv;

        double xmmu2, sigma3inv;
        if ((req_scores & PARAM_SIGMA) || (req_hessian & PARAM_MU_SIGMA))                    sigma3inv = sigma2inv * sigmainv;
        if ((req_scores & PARAM_SIGMA) || (!expected && (req_hessian & PARAM_SIGMA_SIGMA)))  xmmu2 = xmmu * xmmu;

        // Calculate and store score if requested
        if (req_scores & PARAM_MU)      s_mu[i]    = xmmu * sigma2inv;
        if (req_scores & PARAM_SIGMA)   s_sigma[i] = xmmu2 * sigma3inv - sigmainv;

        // Calculate and store hessian if requested
        if (expected) {
            if (req_hessian & PARAM_MU_MU)         h_mu_mu[i]       = -1.0 * sigma2inv;
            if (req_hessian & PARAM_MU_SIGMA)      h_mu_sigma[i]    = 0.0;
            if (req_hessian & PARAM_SIGMA_SIGMA)   h_sigma_sigma[i] = -2.0 * sigma2inv;
        } else {
            if (req_hessian & PARAM_MU_MU)         h_mu_mu[i]       = -1.0 * sigma2inv;
            if (req_hessian & PARAM_MU_SIGMA)      h_mu_sigma[i]    = -2.0 * xmmu * sigma3inv;
            if (req_hessian & PARAM_SIGMA_SIGMA) {
                double sigma4inv = sigma2inv * sigma2inv;
                h_sigma_sigma[i] = -3.0 * xmmu2 * sigma4inv + sigma2inv;
            }
        }

    }

    return R_NilValue;
}
