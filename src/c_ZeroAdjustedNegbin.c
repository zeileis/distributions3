
#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <stdbool.h>
#include <stdlib.h>

/* Probability density function of the Zero-Adjusted Negative Binomial distribution
 *
 * - N: Integer, length of the result
 * - x: Numeric vector of length 1 or N with quantiles
 * - mu,sigma,nu: Parameters of the distribution, numeric of length 1 or N
 * - ret_log: boolean. Should the log-pdf be returned?
 * - ncores: integer, number of cores to use (requires OMP)
 *
 * For y = 0: P(Y=0) = nu
 * For y > 0: P(Y=y) = (1-nu) * P_NBI(y) / (1 - P_NBI(0))
 *           where P_NBI is the negative binomial PMF with size = 1/sigma, mu = mu
 *
 * Return: Numeric vector of length N with probabilities
 */
SEXP c_dZeroAdjustedNegbin(SEXP N, SEXP x, SEXP mu, SEXP sigma, SEXP nu, SEXP ret_log, SEXP ncores) {

    int n = asInteger(N); // length of return vector
    int return_log = asInteger(ret_log);

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    double* xptr = REAL(x);
    double* muptr = REAL(mu);
    double* sigmaptr = REAL(sigma);
    double* nuptr = REAL(nu);

    int i;
    double xi, size_i, mu_i, sigma_i, nu_i;
    double fy0, fy, logfy;

    // Stride variables. Set 1 if we got a vector, else 0.
    int s_x     = (LENGTH(x) == 1)     ? 0 : 1;
    int s_mu    = (LENGTH(mu) == 1)    ? 0 : 1;
    int s_sigma = (LENGTH(sigma) == 1) ? 0 : 1;
    int s_nu    = (LENGTH(nu) == 1)    ? 0 : 1;

    // Allocating vector for return
    SEXP rval; PROTECT(rval = allocVector(REALSXP, n));
    double* rvalptr = REAL(rval);

    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) private(xi, size_i, mu_i, sigma_i, nu_i, fy0, fy, logfy)
    #endif
    for (i = 0; i < n; i++) {
        xi = xptr[i * s_x];
        mu_i = muptr[i * s_mu];
        sigma_i = sigmaptr[i * s_sigma];
        nu_i = nuptr[i * s_nu];
        size_i = 1.0 / sigma_i;

        // Handle negative/non-integer/infinite x
        if (!R_FINITE(xi) || xi < 0.0) {
            rvalptr[i] = return_log ? R_NegInf : 0.0;
            continue;
        }

        // x == 0: P(Y=0) = nu
        if (xi == 0.0) {
            rvalptr[i] = return_log ? log(nu_i) : nu_i;
            continue;
        }

        // x > 0: P(Y=x) = (1-nu) * P_NBI(x) / (1 - P_NBI(0))
        // Compute in log scale for numerical stability
        fy0 = dnbinom_mu(0.0, size_i, mu_i, 1);     // log P_NBI(0)
        fy  = dnbinom_mu(xi, size_i, mu_i, 1);      // log P_NBI(x)
        logfy = log(1.0 - nu_i) + fy - log(1.0 - exp(fy0));

        rvalptr[i] = return_log ? logfy : exp(logfy);
    }

    UNPROTECT(1);
    return rval;
}


/* Cumulative distribution function of the Zero-Adjusted Negative Binomial distribution
 *
 * - N: Integer, length of the result
 * - q: Numeric vector of length 1 or N with quantiles
 * - mu,sigma,nu: Parameters of the distribution
 * - lower_tail,log_p: Flags for CDF transforms
 * - ncores: integer, number of cores
 *
 * F(q) = nu                                      if q = 0
 * F(q) = nu + (1-nu)*(F_NBI(q) - F_NBI(0))/(1-F_NBI(0))  if q > 0
 *
 * Return: Numeric vector of probabilities
 */
SEXP c_pZeroAdjustedNegbin(SEXP N, SEXP q, SEXP mu, SEXP sigma, SEXP nu, SEXP lower_tail, SEXP log_p, SEXP ncores) {

    int n = asInteger(N);
    int lowertail = asInteger(lower_tail);
    int logp = asInteger(log_p);

    #if _OPENMP
    int nthreads = asInteger(ncores);
    #endif

    double* qptr = REAL(q);
    double* muptr = REAL(mu);
    double* sigmaptr = REAL(sigma);
    double* nuptr = REAL(nu);

    int i;
    double qi, size_i, mu_i, sigma_i, nu_i;
    double cdf0, cdf1, cdf;

    // Stride variables
    int s_q     = (LENGTH(q) == 1)     ? 0 : 1;
    int s_mu    = (LENGTH(mu) == 1)    ? 0 : 1;
    int s_sigma = (LENGTH(sigma) == 1) ? 0 : 1;
    int s_nu    = (LENGTH(nu) == 1)    ? 0 : 1;

    // Allocating vector for return
    SEXP rval; PROTECT(rval = allocVector(REALSXP, n));
    double* rvalptr = REAL(rval);

    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) private(qi, size_i, mu_i, sigma_i, nu_i, cdf0, cdf1, cdf)
    #endif
    for (i = 0; i < n; i++) {
        qi = qptr[i * s_q];
        mu_i = muptr[i * s_mu];
        sigma_i = sigmaptr[i * s_sigma];
        nu_i = nuptr[i * s_nu];
        size_i = 1.0 / sigma_i;

        // Handle boundary cases
        if (qi < 0.0) {
            rvalptr[i] = lowertail ? (logp ? R_NegInf : 0.0) : (logp ? 0.0 : 1.0);
            continue;
        }
        if (!R_FINITE(qi)) {
            rvalptr[i] = lowertail ? (logp ? 0.0 : 1.0) : (logp ? R_NegInf : 0.0);
            continue;
        }

        // Compute CDF
        cdf0 = pnbinom_mu(0.0, size_i, mu_i, 1, 0);
        cdf1 = pnbinom_mu(floor(qi), size_i, mu_i, 1, 0);

        if (qi == 0.0) {
            cdf = nu_i;
        } else {
            cdf = nu_i + (1.0 - nu_i) * (cdf1 - cdf0) / (1.0 - cdf0);
        }

        // Apply transformations
        if (!lowertail) cdf = 1.0 - cdf;
        rvalptr[i] = logp ? log(cdf) : cdf;
    }

    UNPROTECT(1);
    return rval;
}


/* Quantile function for Zero-Adjusted Negative Binomial
 *
 * - N: Integer, length of result
 * - p: Numeric vector of length 1 or N with probabilities
 * - mu,sigma,nu: Parameters
 * - lower_tail,log_p: Flags
 * - ncores: integer
 *
 * Uses algebraic inversion (no root-finding needed):
 * pnew = (p - nu) / (1 - nu) - 1e-10
 * cdf0 = F_NBI(0)
 * pnew2 = cdf0*(1-pnew) + pnew
 * q = q_NBI(pnew2)
 *
 * Return: Numeric vector of quantiles
 */
SEXP c_qZeroAdjustedNegbin(SEXP N, SEXP p, SEXP mu, SEXP sigma, SEXP nu, SEXP lower_tail, SEXP log_p, SEXP ncores) {

    int n = asInteger(N);
    int lowertail = asInteger(lower_tail);
    int logp = asInteger(log_p);

    #if _OPENMP
    int nthreads = asInteger(ncores);
    #endif

    double* pptr = REAL(p);
    double* muptr = REAL(mu);
    double* sigmaptr = REAL(sigma);
    double* nuptr = REAL(nu);

    int i;
    double pi, size_i, mu_i, sigma_i, nu_i;
    double pnew, cdf0, pnew2;

    // Stride variables
    int s_p     = (LENGTH(p) == 1)     ? 0 : 1;
    int s_mu    = (LENGTH(mu) == 1)    ? 0 : 1;
    int s_sigma = (LENGTH(sigma) == 1) ? 0 : 1;
    int s_nu    = (LENGTH(nu) == 1)    ? 0 : 1;

    // Allocating vector for return
    SEXP rval; PROTECT(rval = allocVector(REALSXP, n));
    double* rvalptr = REAL(rval);

    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) private(pi, size_i, mu_i, sigma_i, nu_i, pnew, cdf0, pnew2)
    #endif
    for (i = 0; i < n; i++) {
        pi = pptr[i * s_p];

        // Apply log_p and lower_tail transformations
        pi = logp ? exp(pi) : pi;
        pi = lowertail ? pi : 1.0 - pi;

        mu_i = muptr[i * s_mu];
        sigma_i = sigmaptr[i * s_sigma];
        nu_i = nuptr[i * s_nu];
        size_i = 1.0 / sigma_i;

        // Handle boundary and invalid cases
        if (ISNAN(pi) || pi < 0.0 || pi > 1.0) {
            rvalptr[i] = R_NaN;
            continue;
        }
        if (pi == 0.0) {
            rvalptr[i] = 0.0;
            continue;
        }
        if (pi == 1.0) {
            rvalptr[i] = R_PosInf;
            continue;
        }

        // Algebraic inversion of ZANBI quantile
        pnew = (pi - nu_i) / (1.0 - nu_i) - 1e-10;
        cdf0 = pnbinom_mu(0.0, size_i, mu_i, 1, 0);
        pnew2 = cdf0 * (1.0 - pnew) + pnew;
        if (pnew2 < 0.0) pnew2 = 0.0;

        rvalptr[i] = qnbinom_mu(pnew2, size_i, mu_i, 1, 0);
    }

    UNPROTECT(1);
    return rval;
}
