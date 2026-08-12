
#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <stdbool.h>
#include <stdlib.h>

/* Probability density function of the SHASH distribution
 *
 * - n: Integer, lengt of the result. This is, at the same
 *   time also the length of the next couple of vectors in case they
 *   are not equal to 1.
 * - q: Numeric vector of length 1 or n with quantiles.
 * - mu,sigma,nu,tau: The parameters of the distribution, all numeric
 *   of length 1 or n.
 * - logp: boolean. Should the log-pdf be returned?
 *
 * Note: Not faster than the R implementation, but this was my
 * test function for what comes next, especially qsash seems
 * slow in bsae R. Optional: for small samples OMP is also useless
 * and even increases executation time for N = 10,000 (too much overhead)
 *
 * Return
 * ------
 * Returns a double (numeric) vector of length(n) with the probability
 * density of the shash distributions defined by mu, sigma, nu, tau.
 *
 * Reto Stauffer, August 2026
 */
SEXP c_pshash(SEXP N, SEXP q, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP lower_tail, SEXP log_p, SEXP ncores) {

    int n = asInteger(N); // length of return vector
    int logp = asInteger(log_p);
    int lowertail = asInteger(lower_tail);

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    double* qptr = REAL(q);
    double* muptr = REAL(mu);
    double* sigmaptr = REAL(sigma);
    double* nuptr = REAL(nu);
    double* tauptr = REAL(tau);

    int i;
    double z, tmp;

    // Stride varaibles. Set 1 if we got a vector, else 0.
    // Used inside the loop to modify i (use i * 1 or i * 0).
    int s_q     = (LENGTH(q) == 1)     ? 0 : 1;
    int s_mu    = (LENGTH(mu) == 1)    ? 0 : 1;
    int s_sigma = (LENGTH(sigma) == 1) ? 0 : 1;
    int s_nu    = (LENGTH(nu) == 1)    ? 0 : 1;
    int s_tau   = (LENGTH(tau))        ? 0 : 1;

    // Allocating vector for return
    SEXP rval; PROTECT(rval = allocVector(REALSXP, n));
    double* rvalptr = REAL(rval);

    #if _OPENMP
    ///Rprintf(" Using OMP with nthreads = %d\n", nthreads);
    #pragma omp parallel for num_threads(nthreads) private(tmp, z)
    #endif
    for (i = 0; i < n; i++) {
        z   = (qptr[i * s_q] - muptr[i * s_mu]) / sigmaptr[i * s_sigma];
        tmp = asinh(z);
        tmp = exp(tauptr[i * s_tau] * tmp) - exp(-nuptr[i * s_nu] * tmp);
        rvalptr[i] = Rf_pnorm5(0.5 * tmp, 0., 1., lowertail, 0); // x, mean, sd, lower_tail, log_p = 0
        if (logp) rvalptr[i] = log(rvalptr[i]);
    }

    UNPROTECT(1);
    return rval;
}

