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
    PARAM_PI     = 1U << 0,
    PARAM_LAMBDA = 1U << 1,

    // Cross-derivatives (used for hessian)
    PARAM_PI_PI         = 1U << 2,
    PARAM_PI_LAMBDA     = 1U << 3,
    PARAM_LAMBDA_LAMBDA = 1U << 4,
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

/* Expected Hessian for one observation (ZIPoisson)
 *
 * E[H] = sum_{k = 0}^{kmax} H(k) * P(X = k), where kmax is the (1 - 1e-12)
 * quantile, as in the R implementation (quantile(d, 1 - 1e-12)).
 *
 *   P(X = 0) = p0 = pi + (1 - pi) * exp(-lambda)
 *   P(X = k) = (1 - pi) * dpois(k, lambda),  k > 0
 *
 * The observed Hessian at k = 0 is the x0 branch of the observed code. For
 * k > 0 it is H(k) with pi:pi = -1/(1 - pi)^2 (constant in k),
 * lambda:lambda = -k / lambda^2 (linear in k) and pi:lambda = 0.
 *
 * Caller must ensure 0 <= pi < 1 and 0 <= lambda < Inf. The k = 0 term is
 * written in closed form, the k > 0 terms are summed in a loop. Only the
 * requested accumulators are computed.
 */
static inline void expected_hessian(double pi, double lambda, ParameterFlags req,
                                        double *pp, double *pl, double *ll) {
    double onempi  = 1 - pi;            // one minus pi
    double q       = exp(-lambda);      // P(Poisson = 0)
    double p0      = pi + onempi * q;   // P(ZIP = 0)
    double onempiq = onempi * q;        // one minus pi times q

    // Truncation point: smallest k with P(X > k) <= 1e-12. For the ZIP,
    // P(X > k) = (1 - pi) * P_Poisson(X > k), hence the upper-tail quantile
    // of the Poisson at 1e-12 / (1 - pi). If that is >= 1, kmax = 0.
    double tail = 1e-12 / onempi;
    int kmax = (tail < 1) ? (int) qpois(tail, lambda, 0, 0) : 0;

    // k = 0 term. The p0 factors in the observed Hessian cancel with the
    // pmf p0, which gives the simplified closed forms below (no overflow
    // for small p0).
    double acc_pp = 0.0, acc_pl = 0.0, acc_ll = 0.0;
    if (req & PARAM_PI_PI)         acc_pp = -(1 - q) * (1 - q) / p0;
    if (req & PARAM_PI_LAMBDA)     acc_pl = q + (1 - q) * onempiq / p0;
    if (req & PARAM_LAMBDA_LAMBDA) acc_ll = onempiq - onempiq * onempiq / p0;

    // k > 0 terms: observed Hessian times (1 - pi) * dpois(k, lambda).
    // pi:lambda is zero for x > 0, so it needs no loop.
    if (req & (PARAM_PI_PI | PARAM_LAMBDA_LAMBDA)) {
        double lambda2 = lambda * lambda;
        for (int k = 1; k <= kmax; k++) {
            double pm = onempi * dpois((double) k, lambda, 0);
            if (pm == 0.0) continue;    // avoids 0 * Inf when lambda = 0
            if (req & PARAM_PI_PI)         acc_pp += -pm / (onempi * onempi);
            if (req & PARAM_LAMBDA_LAMBDA) acc_ll += -((double) k) * pm / lambda2;
        }
    }

    *pp = acc_pp;
    *pl = acc_pl;
    *ll = acc_ll;
}

SEXP c_deriv_ZIPoisson(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP expected_sexp, SEXP ncores) {

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // Epsilon (sqrt of machine doulbe eps)
    double eps = 1.490116e-08;

    // Expected hessian?
    bool expected = LOGICAL(expected_sexp)[0];

    // Setup Input Parameters (pointers + strides) for all parameters of the distribution.
    InputVector x_vec      = make_input_var(x_sexp, NULL, N);
    InputVector pi_vec     = make_input_var(params_sexp, "pi",     N);
    InputVector lambda_vec = make_input_var(params_sexp, "lambda", N);

    // Setting up flags for score and hessian; used for fast bitwise operations
    // to see what needs to be calculated. 'req_scores' means 'required scores',
    // analogously 'req_hessian' is used to check what is required to get the
    // requested hessian elements.
    ParameterFlags req_scores  = get_flags(score_sexp, false);
    ParameterFlags req_hessian = get_flags(hessian_sexp, true);

    // Initialize double pointer for the corresponding `score_sexp` list
    // element to be modified and returned to R.
    double* s_pi        = NULL;
    double* s_lambda     = NULL;
    if (req_scores & PARAM_PI)           s_pi     = getListElement(score_sexp, "pi");
    if (req_scores & PARAM_LAMBDA)       s_lambda = getListElement(score_sexp, "lambda");

    // As for score: Initialize double pointer for required elements corresponding
    // to the elements in hessian_sexp to be modified and returned to R.
    double* h_pi_pi     = NULL;
    double* h_pi_lambda = NULL;
    if (req_hessian & PARAM_PI_PI)       h_pi_pi     = getListElement(hessian_sexp, "pi");
    if (req_hessian & PARAM_PI_LAMBDA)   h_pi_lambda = getListElement(hessian_sexp, "pi:lambda");

    double* h_lambda_lambda = NULL;
    if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda    = getListElement(hessian_sexp, "lambda");

    // Main loop
    #if _OPENMP
    //Rprintf("number of threads: %d\n", nthreads);
    #pragma omp parallel for num_threads(nthreads) if (N > 20000)
    #endif
    for (int i = 0; i < N; i++) {
        // Extracting numeric elements
        double x      = get_val(&x_vec, i);
        double pi     = get_val(&pi_vec, i);
        double lambda = get_val(&lambda_vec, i);

        // x (only for the observed case; the expected Hessian does not use x)
        // or parameter missing? Store NA_REAL and continue
        if ((!expected && ISNAN(x)) || ISNAN(pi) || ISNAN(lambda)) {
            if (req_scores & PARAM_PI)             s_pi[i]            = R_NaN;
            if (req_scores & PARAM_LAMBDA)         s_lambda[i]        = R_NaN;
            if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = R_NaN;
            if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda[i]     = R_NaN;
            if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = R_NaN;
            continue;
        }

        // If x is < 0 or not an integer, continue (result = 0)
        if (x < 0 || fabs(x - round(x)) > eps) { continue; }

        double q  = exp(-lambda);
        double p0 = pi + (1 - pi) * q;
        bool   x0 = round(x) == 0;

        double onempi = 1 - pi;         // one minus pi
        double onempiq = onempi * q;    // one minus pi times q

        // Calculating score if requested
        if (x0) {
            if (req_scores & PARAM_PI)      s_pi[i]     = (1 - q) / p0;
            if (req_scores & PARAM_LAMBDA)  s_lambda[i] = -onempiq / p0;
        } else {
            if (req_scores & PARAM_PI)      s_pi[i]     = -1 / onempi;
            if (req_scores & PARAM_LAMBDA)  s_lambda[i] = x / lambda - 1;
        }

        // Expected Hessian: sum over the support, independent of x. The score
        // is not computed in this mode (the R wrapper never requests it).
        if (expected) {
            // Outside the domain the sum is undefined: return NaN, and do not
            // call the Rmath routines (they may warn, which is not thread-safe).
            if (!(pi >= 0 && pi <= 1 && lambda >= 0 && R_FINITE(lambda))) {
                if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = R_NaN;
                if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda[i]     = R_NaN;
                if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = R_NaN;
                continue;
            }

            double pp, pl, ll;
            expected_hessian(pi, lambda, req_hessian, &pp, &pl, &ll);
            if (req_hessian & PARAM_PI_PI)         h_pi_pi[i]         = pp;
            if (req_hessian & PARAM_PI_LAMBDA)     h_pi_lambda[i]     = pl;
            if (req_hessian & PARAM_LAMBDA_LAMBDA) h_lambda_lambda[i] = ll;
        // Observed hessian
        } else {
            // Observed Hessian
            if (x0) {
                if (req_hessian & PARAM_PI_PI)          h_pi_pi[i]         = -pow((1 - q) / p0, 2);
                if (req_hessian & PARAM_PI_LAMBDA)      h_pi_lambda[i]     = q / p0 + (1 - q) * onempiq / pow(p0, 2);
                if (req_hessian & PARAM_LAMBDA_LAMBDA)  h_lambda_lambda[i] = onempiq / p0 - pow(onempiq / p0, 2);
            } else {
                if (req_hessian & PARAM_PI_PI)          h_pi_pi[i]         = -1 / pow(onempi, 2);
                if (req_hessian & PARAM_LAMBDA_LAMBDA)  h_lambda_lambda[i] = -x / pow(lambda, 2);
                // h_pi_lambda[i] is always 0, no calculation needed
            }
        }
    }

    return R_NilValue;
}
