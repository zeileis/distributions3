
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>
#include "distributions3.h"


/* Calculating moments based on a PDF vector using trapezoidal rule (continuous distributions)
 *
 * @param i    index of distribution (row index for 'at').
 * @param atconst logical, whether 'at' is a vector used for all distributions
 *             (i.e., rows in 'p') or a matrix (same dimension as p).
 * @param dim  dimension of matrices 'at' and 'p'.
 * @param at   numeric matrix of grid position (n x w->length) where the PDF was
 *             provided. The weights are stored on 'w'.
 * @param p    numeric matrix with PDF evaluated at 'at', same dimension as 'at'.
 * @param what 1 = mean, 2 = variance, 3 = skewness, 4 = kurtosis
 *
 * Returns single numeric value.
 */
double c_moments_calculate_trapezoidal(int i, double* p, int np, double* q, int nq, int* dim, int what, doubleVec* width, doubleVec* dens, doubleVec* mids) {

    // Extracting matrix dimensions for convenience
    int n = dim[0], k = dim[1];

    bool q_is_vector = (nq == k) ? true : false;
    bool p_is_vector = (np == k) ? true : false;

    // cdfint: integral of the CDF to check if the density (pdf; p) provided looks OKish
    double cdfint = 0.0, q_lo, q_hi, p_lo, p_hi;
    for (int j = 0; j < (k - 1); j++) {
        q_lo = q_is_vector ? q[j]     : q[i + n * j];
        q_hi = q_is_vector ? q[j + 1] : q[i + n * (j + 1)];

        p_lo = p_is_vector ? p[j]     : p[i + n * j];
        p_hi = p_is_vector ? p[j + 1] : p[i + n * (j + 1)];

        width->values[j] = q_hi - q_lo;
        mids->values[j]  = (q_hi + q_lo) * 0.5;
        dens->values[j]  = (p_hi - p_lo);

        //testing// Rprintf(" --- j=%d: q = [%10.5f, %10.5f] mid = %10.5f ... w = %10.5f, dens = %10.5f\n", j, q_lo, q_hi, mids->values[j], width->values[j], dens->values[j]);
        cdfint += dens->values[j];
    }
    //testing// Rprintf(" --- [cdf integral] cdfint = %.5f\n", cdfint);

    // TODO(R): Calculating the area under the curve which - analytically - should be 1.0.
    //          Here I am just stopping if the area is <= 1e-9 and makes no sense at all, however,
    //          we _could_ check if the area under the curve (distribution) is >= .95 or similar
    //          to ensure that the grid ('at') covers enough of the distribution to get reliable
    //          estimates of the central moments.
    if (cdfint <= 1e-9) { Rf_error("Error: The area under the 'pdf' curve is zero or negative."); }

    //// Calculating mean
    double mean = 0.0;
    for (int j = 0; j < (k - 1); j++) {
        mean += mids->values[j] * dens->values[j];
    }
    mean = mean / cdfint; // Normalize
    if (what == 1) { return mean; }


    // Calculating variance
    double variance = 0.0;
    for (int j = 0; j < (k - 1); j++) {
        variance += pow(mids->values[j] - mean, 2.0) * dens->values[j];
    }
    variance = variance / cdfint; // Normalize
    if (what == 2) { return variance; }

    // Calculating skewness or kurtosis
    double powexp = (double)what;   // Define power parameter (3 for skewness, 4 = kurtosis)
    double res = 0.0;

    for (int j = 0; j < (k - 1); j++) {
       res += pow(mids->values[j] - mean, powexp) * dens->values[j];
    }
    res = res / cdfint / pow(variance, powexp / 2.0); // Normalize and divide
    if (what == 4) { res = res - 3.0; }

    return res;
}


double c_moments_calculate_discrete(int i, double* p, int np, double* q, int nq, int* dim, int what, doubleVec* dens, doubleVec* mids) {

    // Extracting matrix dimensions for convenience
    int n = dim[0], k = dim[1];

    bool q_is_vector = (nq == k) ? true : false;
    bool p_is_vector = (np == k) ? true : false;

    // Normalization using Trapezoidal Rule
    // cdfint: integral of the CDF to check if the density (pdf; p) provided looks OKish
    double cdfint = 0.0, p_lo, p_hi;
    for (int j = 0; j < k; j++) {
        mids->values[j] = q_is_vector ? q[j] : q[i + n * j]; // not actually 'mid' but quanile (count data value)
        if (j == 0) {
            dens->values[j] = p_is_vector ? p[j] : p[i + n * j];
        } else {
            p_lo = p_is_vector ? p[j - 1] : p[i + n * (j - 1)];
            p_hi = p_is_vector ? p[j]     : p[i + n * j];
            dens->values[j]  = (p_hi - p_lo); // * 1.0 (we assume width 1, count data)
        }

        //testing// Rprintf(" --- j=%d: mid (count) = %10.5f ... dens = %10.5f\n", j, mids->values[j], dens->values[j]);
        cdfint         += dens->values[j];
    }
    //testing// Rprintf(" --- [cdf integral] cdfint = %.5f\n", cdfint);

    // TODO(R): Calculating the area under the curve which - analytically - should be 1.0.
    //          Here I am just stopping if the area is <= 1e-9 and makes no sense at all, however,
    //          we _could_ check if the area under the curve (distribution) is >= .95 or similar
    //          to ensure that the grid ('at') covers enough of the distribution to get reliable
    //          estimates of the central moments.
    if (cdfint <= 1e-9) { Rf_error("Error: The area under the 'pdf' curve is zero or negative."); }

    //// Calculating mean
    double mean = 0.0;
    for (int j = 0; j < (k - 1); j++) {
        mean += mids->values[j] * dens->values[j];
    }
    mean = mean / cdfint; // Normalize
    if (what == 1) { return mean; }


    // Calculating variance
    double variance = 0.0;
    for (int j = 0; j < (k - 1); j++) {
        variance += pow(mids->values[j] - mean, 2.0) * dens->values[j];
    }
    variance = variance / cdfint; // Normalize
    if (what == 2) { return variance; }

    // Calculating skewness or kurtosis
    double powexp = (double)what;   // Define power parameter (3 for skewness, 4 = kurtosis)
    double res = 0.0;

    for (int j = 0; j < (k - 1); j++) {
       res += pow(mids->values[j] - mean, powexp) * dens->values[j];
    }
    res = res / cdfint / pow(variance, powexp / 2.0); // Normalize and divide
    if (what == 4) { res = res - 3.0; }

    return res;
}


/*
 * Numerically integrating the density to evaluate the cummulative distribution function (CDF)
 * given we know nothing else than the PDF of a distribution.
 * While 'p' is always a matrix of dimension ('number of distributions x grid to evalate on (at)')
 * 'at' can be a matrix of the very same dimension (if the grid differs between distributions) or
 * a numeric vector with 'grid to evaluate on' re-used for all distributions. This is used for
 * discrete distributions where the grid we work on are typically counts (0, 1, 2, 3, 4, ...) and
 * therefore identical for all distributions evaluated.
 *
 * at      numeric vector (matrix), grid points along the x-axis of the distribution.
 * p       numeric vector (matrix) with the density at grid points 'x'.
 * int     integer of length 2, dimension of matrices at/p
 * x       numeric, where to evaluate the CDF.
 * whatint integer which defines what to calculate (hardcoded).
 *         1 = mean, 2 = variance, ... TODO(R) extend.
 *
 * Dimensions
 * ----------
 * Input 'dim' contains the dimension of the matrices, whereof we will
 * use dim[0] = n and dim[1] = k, i.e., assuming matrix dimension (n x k).
 */
SEXP c_moments_numeric(SEXP p, SEXP q, SEXP dim, SEXP discrete, SEXP whatint) {

    // Setting up pointers for at, pat, dim, x
    double *pptr = REAL(p);
    double *qptr = REAL(q);

    // Getting matrix dimension
    int *dimptr = INTEGER(dim);
    int n = dimptr[0], k = dimptr[1];

    // Number of elements in p and q. Used to check if either of them are
    // matrices or vectors (also handled as 'vectors' if single-row matrix).
    int np = length(p), nq = length(q);

    // Are we dealing with a discrete distribution or not? If not (continuous
    // or continuous approximation of a discrete distribution) we use the
    // trapezoidal rule, else we assume the interval in 'at' is 1.0 (count data)
    // and can sum up.
    bool is_discrete = INTEGER(discrete)[0] > 0 ? true : false;

    // What to calculate?
    int what = INTEGER(whatint)[0];

    // Loop index
    int i;
    int nprotected = 0;

    //// Setting up results vector (length n)
    SEXP rval;
    PROTECT(rval = allocVector(REALSXP, n));
    double *rvalptr = REAL(rval);
    nprotected += 1;

    // Allocating three vectors of length 'k - 1' (number of intervals) to
    // store interval width, approximated density, and interval mids. Used by
    // c_moments_calculate_*(). Note: If is_discrete these vectors are of length 'k',
    // else of length 'k - 1'.
    doubleVec width, dens, mids;
    int kvec = is_discrete ? k :k - 1;

    // Interval width, interval density, and interval mids.
    width.values = (double*)malloc((kvec) * sizeof(double));   width.length = kvec;
    dens.values  = (double*)malloc((kvec) * sizeof(double));   dens.length  = kvec;
    mids.values  = (double*)malloc((kvec) * sizeof(double));   mids.length  = kvec;

    // Looping trough all observations i = 1, ..., n;
    for (i = 0; i < n; i++) {
        // Calculate requested moment
        if (is_discrete) {
            rvalptr[i] = c_moments_calculate_discrete(i, pptr, np, qptr, nq, dimptr, what, &dens, &mids);
        } else {
            rvalptr[i] = c_moments_calculate_trapezoidal(i, pptr, np, qptr, nq, dimptr, what, &width, &dens, &mids);
        }
    }

    // Free memory
    free(width.values); free(dens.values); free(mids.values);

    // Releasing protected object(s) and return
    UNPROTECT(nprotected);
    return rval;
}

