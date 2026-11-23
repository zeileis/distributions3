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
 * HurdlePoisson has two parameters: lambda (Poisson part) and pi (probability
 * of a non-zero observation). Note that the hessian names follow the parameter
 * order of the distribution, i.e., the cross-derivative is named "pi:lambda".
 */
typedef enum {
    PARAM_NONE  = 0,

    // First derivative/main parameters (used for score and hessian)
    PARAM_PI     = 1U << 0,
    PARAM_LAMBDA = 1U << 1,

    // Cross-derivatives (used for hessian)
    PARAM_PI_PI         = 1U << 2,
    PARAM_PI_LAMBDA     = 1U << 3,
    PARAM_LAMBDA_LAMBDA = 1U << 4,
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
        if (getListElementSEXP(x, "pi")     != R_NilValue) req |= PARAM_PI;
        if (getListElementSEXP(x, "lambda") != R_NilValue) req |= PARAM_LAMBDA;
    // Hessian
    } else {
        if (getListElementSEXP(x, "pi")        != R_NilValue) req |= (PARAM_PI | PARAM_PI_PI);
        if (getListElementSEXP(x, "pi:lambda") != R_NilValue) req |= (PARAM_PI | PARAM_LAMBDA | PARAM_PI_LAMBDA);
        if (getListElementSEXP(x, "lambda")    != R_NilValue) req |= (PARAM_LAMBDA | PARAM_LAMBDA_LAMBDA);
    }

    return req;
}

SEXP c_deriv_HurdlePoisson(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP expected_sexp, SEXP ncores) {

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
    InputVector pi_vec     = make_input_var(params_sexp, "pi",     N);
    InputVector lambda_vec = make_input_var(params_sexp, "lambda", N);

    // Required flags for score and hessian
    ParameterFlags req_scores  = get_flags(score_sexp, false);
    ParameterFlags req_hessian = get_flags(hessian_sexp, true);

    // Output pointers (NULL if not requested)
    double* s_pi     = NULL;
    double* s_lambda = NULL;
    if (req_scores & PARAM_PI)     s_pi     = getListElement(score_sexp, "pi");
    if (req_scores & PARAM_LAMBDA) s_lambda = getListElement(score_sexp, "lambda");

    double* h_pi_pi         = NULL;
    double* h_pi_lambda     = NULL;
    double* h_lambda_lambda = NULL;
    if (req_hessian & PARAM_PI_PI)         h_pi_pi         = getListElement(hessian_sexp, "pi");
    if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda     = getListElement(hessian_sexp, "pi:lambda");
    if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda = getListElement(hessian_sexp, "lambda");

    // Main loop
    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) if (N > 20000)
    #endif
    for (int i = 0; i < N; i++) {
        double x      = get_val(&x_vec, i);
        double pi     = get_val(&pi_vec, i);
        double lambda = get_val(&lambda_vec, i);

        bool x0 = round(x) <= 0 || fabs(x - round(x)) > eps;  // Check if x equals 0

        // Missing pi or lambda: Set all to NaN and continue
        if (ISNAN(pi) || ISNAN(lambda)) {
            if (req_scores & PARAM_PI)             s_pi[i]            = R_NaN;
            if (req_scores & PARAM_LAMBDA)         s_lambda[i]        = R_NaN;
            if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = R_NaN;
            if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda[i]     = R_NaN;
            if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = R_NaN;
            continue;
        }

        // Calculating score. If x is missing, store NaN.
        if (ISNAN(x)) {
            if (req_scores & PARAM_PI)      s_pi[i]            = R_NaN;
            if (req_scores & PARAM_LAMBDA)  s_lambda[i]        = R_NaN;
        } else {
            if (req_scores & PARAM_PI)      s_pi[i]     = x0 ? -1. / (1. - pi) : 1. / pi;
            if (req_scores & PARAM_LAMBDA)  s_lambda[i] = x0 ?  0. : x / lambda - 1 - 1 / expm1(lambda);
        }

        // Expected Hessian: sum over the support, independent of x. Closed form:
        //   E[H_pp] = -(1 - pi)^{-1} - pi^{-1}          (P(0) = 1 - pi)
        //   E[H_pl] = 0
        //   E[H_ll] = pi * (q / (1 - q)^2 - 1 / (lambda (1 - q)))
        // i.e. pi times the zero-truncated Poisson expected hessian, because
        // the non-zero part has probability pi. The score is not computed here.
        if (expected) {
            if (pi < 0 || pi > 1 || lambda <= 0 || !R_FINITE(lambda)) {
                if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = R_NaN;
                if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda[i]     = R_NaN;
                if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = R_NaN;
            } else {
                if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda[i]     = 0.0;
                if (req_hessian & PARAM_LAMBDA_LAMBDA) {
                    double om = -expm1(-lambda);   // 1 - q, q = 1 - om (one transcendental call)
                    double q  = 1 - om;
                    h_lambda_lambda[i] = pi * (q / (om * om) - 1 / (lambda * om));
                }

                // Special handling of pi = 0 and pi = 1
                if (pi > eps && pi < (1 - eps)) {
                    if (req_hessian & PARAM_PI_PI)  h_pi_pi[i] = -1 / (1 - pi) - 1 / pi;
                } else {
                    // pi < 0.5 equals True if pi is ~0, and False if pi is ~1
                    if (req_hessian & PARAM_PI_PI)  h_pi_pi[i] = pi < 0.5 ? -1 : R_NaN;
                }
            }
        // Observed case: x < 0 or non-integer gives 0 (outside support)
        } else {
            // h_pi_lambda is always 0.0; nothing to do
            if (x0) {
                if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = -1. / pow(1. - pi, 2);
                // h_lambda_lambda is 0 if x = 0; nothing to do
            } else {
                if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = -1. / (pi * pi);
                if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = -x / (lambda * lambda) + exp(lambda) / pow(expm1(lambda), 2);
            }
        }
    }

    return R_NilValue;
}
