#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "deriv.h"

/* Bitwise Flags for Required Parameters
 *
 * ZTPoisson has a single parameter (lambda), so only the lambda
 * score/hessian and the lambda:lambda cross-derivative exist.
 */
typedef enum {
    PARAM_NONE   = 0,

    // First derivative/main parameters (used for score and hessian)
    PARAM_LAMBDA = 1U << 0,

    // Cross-derivatives (used for hessian)
    PARAM_LAMBDA_LAMBDA = 1U << 1,
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
        if (getListElementSEXP(x, "lambda") != R_NilValue) req |= PARAM_LAMBDA;
    // Hessian
    } else {
        if (getListElementSEXP(x, "lambda") != R_NilValue) req |= (PARAM_LAMBDA | PARAM_LAMBDA_LAMBDA);
    }

    return req;
}

SEXP c_deriv_ZTPoisson(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP expected_sexp, SEXP ncores) {

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // Epsilon (sqrt of machine double eps)
    double eps = 1.490116e-08;

    // Expected hessian?
    bool expected = LOGICAL(expected_sexp)[0];

    // Setup Input Parameters (pointers + strides)
    InputVector x_vec      = make_input_var(x_sexp, NULL, N);
    InputVector lambda_vec = make_input_var(params_sexp, "lambda", N);

    // Required flags for score and hessian
    ParameterFlags req_scores  = get_flags(score_sexp, false);
    ParameterFlags req_hessian = get_flags(hessian_sexp, true);

    // Output pointers (NULL if not requested)
    double* s_lambda = NULL;
    if (req_scores & PARAM_LAMBDA)           s_lambda        = getListElement(score_sexp, "lambda");

    double* h_lambda_lambda = NULL;
    if (req_hessian & PARAM_LAMBDA_LAMBDA)   h_lambda_lambda = getListElement(hessian_sexp, "lambda");

    // Main loop
    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) if (N > 20000)
    #endif
    for (int i = 0; i < N; i++) {
        double x      = get_val(&x_vec, i);
        double lambda = get_val(&lambda_vec, i);

        // Missing lambda always resulting in NaN
        if (ISNAN(lambda) || (lambda < 0 || !R_FINITE(lambda))) {
            if (req_scores & PARAM_LAMBDA)         s_lambda[i]        = R_NaN;
            if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = R_NaN;
            continue;
        }

        // If x is missing: Neither score nor observed Hessian can be calculated
        if (!expected & ISNAN(x)) {
            if (req_scores & PARAM_LAMBDA)           s_lambda[i]        = R_NaN;
            if (req_hessian & PARAM_LAMBDA_LAMBDA)   h_lambda_lambda[i] = R_NaN;
            continue;
        }

        // Calculating score and hessian
        double q  = exp(-lambda);
        double om = 1 - q;          // 1 - q
        double r  = q / om;         // 1 / (1 - q)
        bool   x0 = round(x) == 0;

        // Calculating score:  x / lambda - 1 - q / (1 - q)
        if (req_scores & PARAM_LAMBDA) {
            if (ISNAN(x)) {
                s_lambda[i] = R_NaN;
            } else {
                s_lambda[i] = x0 ? 0.0: x / lambda - 1 - r;
            }
        }

        // Expected hessian
        if (req_hessian & PARAM_LAMBDA_LAMBDA) {
            if (expected) {
                h_lambda_lambda[i] = q / (om * om) - 1 / (lambda * om);
            } else {
                // Observed case: x < 1 or non-integer gives 0 (outside support).
                // Nothing has to be done as s_lambda and h_lambda_lambda were initialized with zeros.
                if (x < 1 || fabs(x - round(x)) > eps) {
                    h_lambda_lambda[i] = 0.0;
                } else {
                    // hessian: -x / lambda^2 + q / (1 - q)^2
                    h_lambda_lambda[i] = -x / (lambda * lambda) + q / (om * om);
                }
            }
        }
    }

    return R_NilValue;
}
