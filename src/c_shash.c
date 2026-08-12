
#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <stdbool.h>
#include <stdlib.h>

/* Distribution function of the SHASH distribution
 *
 * - n: Integer, length of the result. This is, at the same
 *   time also the length of the next couple of vectors in case they
 *   are not equal to 1.
 * - q: Numeric vector of length 1 or n with quantiles.
 * - mu,sigma,nu,tau: The parameters of the distribution, all numeric
 *   of length 1 or n.
 * - log_p: boolean. Should the log-pdf be returned?
 * - ncores: integer, number of cores to use (requires OMP to be
 *   available)
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

/* Probability density function of the SHASH distribution
 *
 * - n: Integer, length of the result. This is, at the same
 *   time also the length of the next couple of vectors in case they
 *   are not equal to 1.
 * - x: Numeric vector of length 1 or n with quantiles.
 * - mu,sigma,nu,tau: The parameters of the distribution, all numeric
 *   of length 1 or n.
 * - ret_log: boolean. Should the log be returned?
 * - ncores: integer, number of cores to use (requires OMP to be
 *   available)
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
SEXP c_dshash(SEXP N, SEXP x, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP ret_log, SEXP ncores) {

    int n = asInteger(N); // length of return vector
    int return_log = asInteger(ret_log);

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    double* xptr = REAL(x);
    double* muptr = REAL(mu);
    double* sigmaptr = REAL(sigma);
    double* nuptr = REAL(nu);
    double* tauptr = REAL(tau);

    int i;
    double z, r, c, loglik;
    double asinhz, exp_tauasinhz, exp_minusnuasinhz;

    // Stride varaibles. Set 1 if we got a vector, else 0.
    // Used inside the loop to modify i (use i * 1 or i * 0).
    int s_x     = (LENGTH(x) == 1)     ? 0 : 1;
    int s_mu    = (LENGTH(mu) == 1)    ? 0 : 1;
    int s_sigma = (LENGTH(sigma) == 1) ? 0 : 1;
    int s_nu    = (LENGTH(nu) == 1)    ? 0 : 1;
    int s_tau   = (LENGTH(tau))        ? 0 : 1;

    // Allocating vector for return
    SEXP rval; PROTECT(rval = allocVector(REALSXP, n));
    double* rvalptr = REAL(rval);

    #if _OPENMP
    ///Rprintf(" Using OMP with nthreads = %d\n", nthreads);
    #pragma omp parallel for num_threads(nthreads) private(z, asinhz, exp_tauasinhz, exp_minusnuasinhz, r, c, loglik)
    #endif
    for (i = 0; i < n; i++) {
        z                 = (xptr[i * s_x] - muptr[i * s_mu]) / sigmaptr[i * s_sigma];
        rvalptr[i] = z;
        asinhz            = asinh(z);
        exp_tauasinhz     = exp(tauptr[i * s_tau] * asinhz);
        exp_minusnuasinhz = exp(-nuptr[i * s_nu] * asinhz);

        r = 0.5 * (exp_tauasinhz - exp_minusnuasinhz);
        c = 0.5 * (tauptr[i * s_tau] * exp_tauasinhz + nuptr[i * s_nu] * exp_minusnuasinhz);

        loglik = -log(sigmaptr[i * s_sigma]) - 0.5 * log(2 * M_PI) - 0.5 * log(1 + z * z) + log(c) - 0.5 * r * r;
        rvalptr[i] = (return_log) ? loglik : exp(loglik);
    }

    UNPROTECT(1);
    return rval;
}



struct qshash_args {
    double r;
    double nu;
    double tau;
};

// Objective function for root-finding: 0.5 * (exp(tau*w) - exp(-nu*w)) - r = 0
static double qshash_obj(double w, void *info) {
    struct qshash_args *args = (struct qshash_args *)info;
    return 0.5 * (exp(args->tau * w) - exp(-args->nu * w)) - args->r;
}


/* Self-contained Brent's zeroin algorithm
 * This was written by AI */
static double local_zeroin(double ax, double bx, double (*f)(double x, void *info), void *info, double tol) {
    double a = ax, b = bx, c = a;
    double fa = f(a, info), fb = f(b, info), fc = fa;
    double d = 0.0, e = 0.0;
    double eps = 2.2204460492503131e-16; // Machine epsilon
    int max_iter = 1000;
    int iter = 0;

    if (fa == 0.0) return a;
    if (fb == 0.0) return b;

    while (iter < max_iter) {
        iter++;

        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa;
            d = b - a; e = d;
        }

        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double tol1 = 2.0 * eps * fabs(b) + 0.5 * tol;
        double xm = 0.5 * (c - b);

        if (fabs(xm) <= tol1 || fb == 0.0) {
            return b;
        }

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            double s = fb / fa;
            double p, q;

            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (p > 0.0) q = -q;
            p = fabs(p);

            double min1 = 3.0 * xm * q - fabs(tol1 * q);
            double min2 = fabs(e * q);

            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;

        if (fabs(d) > tol1) {
            b += d;
        } else {
            b += (xm > 0.0) ? tol1 : -tol1;
        }

        fb = f(b, info);
    }
    return b;
}


/* Quantile function for SHASH
 *
 * - N: Integer, length of the result. This is, at the same
 *   time also the length of the next couple of vectors in case they
 *   are not equal to 1.
 * - p: Numeric vector of length 1 or n with probabilities.
 * - mu,sigma,nu,tau: The parameters of the distribution, all numeric
 *   of length 1 or n.
 * - lower_tail,log_p: Indicate whether or not 'p' is given as
 *   upper tail or on the log scale, required to first transform 'p' before
 *   calculating the quantiles.
 * - ncores: integer, number of cores to use (requires OMP to be
 *   available)
 *
 * Return
 * ------
 * Numeric vecotr of length N. Can contain missing values (if p outside
 * interval [0, 1]) as well as -Inf/+Inf if p falls on 0 or 1.
 *
 * Reto Stauffer, August 2026
 */
SEXP c_qshash(SEXP N, SEXP p, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP lower_tail, SEXP log_p, SEXP ncores) {

    int n         = asInteger(N); // length of return vector
    int lowertail = asInteger(lower_tail);
    int logp      = asInteger(log_p);

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    double* pptr = REAL(p);
    double* muptr = REAL(mu);
    double* sigmaptr = REAL(sigma);
    double* nuptr = REAL(nu);
    double* tauptr = REAL(tau);

    int i;

    // Stride varaibles. Set 1 if we got a vector, else 0.
    // Used inside the loop to modify i (use i * 1 or i * 0).
    int s_p     = (LENGTH(p) == 1)     ? 0 : 1;
    int s_mu    = (LENGTH(mu) == 1)    ? 0 : 1;
    int s_sigma = (LENGTH(sigma) == 1) ? 0 : 1;
    int s_nu    = (LENGTH(nu) == 1)    ? 0 : 1;
    int s_tau   = (LENGTH(tau))        ? 0 : 1;

    // Allocating vector for return
    SEXP rval; PROTECT(rval = allocVector(REALSXP, n));
    double* rvalptr = REAL(rval);

    double ax    = -20.0;
    double bx    = 20.0;
    double tol   = 1e-12;

    double p_i, r, w;

    #if _OPENMP
    //Rprintf(" Using OMP with nthreads = %d\n", nthreads);
    #pragma omp parallel for num_threads(nthreads) private(p_i, r, w)
    #endif
    for (i = 0; i < n; i++) {
        p_i = pptr[i * s_p];
        p_i = (logp) ? exp(p_i) : p_i;
        p_i = (lowertail) ? p_i : 1. - p_i;
        rvalptr[i] = p_i;

        // Handling cases where p is outside range or on the boundary
        if (ISNAN(p_i) || p_i < 0.0 || p_i > 1.0) { rvalptr[i] = R_NaN; continue; }
        if (p_i == 0.0) { rvalptr[i] = R_NegInf; continue; }
        if (p_i == 1.0) { rvalptr[i] = R_PosInf; continue; }

        // Calculating quantile from standard normal distribution
        r = Rf_qnorm5(p_i, 0.0, 1.0, 1, 0);

        if (nuptr[i * s_nu] == tauptr[i * s_tau]) {
            // Closed-form exact solution when nu == tau (symmetric case)
            w = asinh(r) / tauptr[i * s_tau];
            rvalptr[i] = muptr[i * s_mu] + sigmaptr[i * s_sigma] * sinh(w);
        } else {
            // General asymmetric case: root finding on transformed space w = asinh(z)
            struct qshash_args args = { r, nuptr[i * s_nu], tauptr[i * s_tau] };

            // R_zeroin2 is R's internal zeroin routine
            w = local_zeroin(ax, bx, qshash_obj, (void *)&args, tol);
            rvalptr[i] = muptr[i * s_mu] + sigmaptr[i * s_sigma] * sinh(w);
        }
    }

    UNPROTECT(1);
    return rval;
}

