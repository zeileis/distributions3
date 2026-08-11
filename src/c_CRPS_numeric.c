
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>

/* The following graph tries to show a CDF, defined by
 * a series of 'n' points where the quantile (q) is known
 * (input q) as well as the corresponding value from the
 * cummulative distribution function (input p) which define
 * a series of 'bins' from quantiles close to 0 up to close to 1.
 *
 * The CRPS is defined as:
 *
 *           integral((cdf(x) - H(x))^2, dx)
 *
 * ... which approximated in this function using the trapezoidal rule.
 *
 *
 *
 *              1.0 |                   --------
 *                  |                 /
 *                  |                /
 * cdf(x)     p[up] |               +
 *                  |              /
 *                  |             /
 *            p[lo] | ------+-----
 *              0.0 ------------------------------ > x
 *                          |       |       |
 *                        q[lo]   q[up]     x
 *
 * Three different cases are possible:
 * if    p[pi_up] < px[i]       width * (p[lo]^2 + p[up]^2) / 2
 * if    p[pi_lo] > px[i]       width * ((1. - p[lo])^2 + (1. - p[up])^2 / 2
 * else the bin contains the observation, wherefore the contribution
 * to the CRPS is:
 *      width * (p[lo]^2 + p[x]^2) / 2
 *      width * ((1. - p[x])^2 + (1. - p[up])^2 / 2
 * 
 * x: observation; if length(x) == 1 but length px > 1 the same observation will
 *    be used for all observations.
 * px: cdf at observation (vector of length(y))
 * p: cdf along a grid; matrix of dimension (length(y) x K), or vector of length(y)
 * q: quantile along a grid; matrix of dimension (length(y) x K)
 * continuous: integer, 1 = continuous distribution, 0 = discrete
 *
 * December 2022
 * Reto Stauffer
 */
SEXP c_CRPS_numeric(SEXP x, SEXP px, SEXP p, SEXP q, SEXP continuous) {
    int i, ix, j;     /* Loop/element incides */
    int J;            /* Dimension of matrices (if any; p or q) */
    int pi_lo, pi_up; /* p index lower upper grid */
    int qi_lo, qi_up; /* q index lower upper grid */

    /* x    observations
     * nx   number of observations
     * px   the probability at observation 'x' given a distribution.
     * n    length of px; number of points to be evaluated
     * It is possible that nx = 1 and px > 1 if we evaluate multiple
     * distributions at one very specific 'x'.
     *
     */
    int n  = length(px); /* Number of observations */
    int nx = length(x);  /* NOTE: can be of length 1 or length n */

    double *xptr  = REAL(x);
    double *pxptr = REAL(px);
    double *pptr = REAL(p); /* pdf matrix */
    double *qptr = REAL(q); /* q matrix */

    /* Result */
    SEXP rval;
    PROTECT(rval = allocVector(REALSXP, n));
    double *rvalptr = REAL(rval);

    double w, h, h_lo, h_up, tmp, tmp_lo, tmp_up;

    /* Continuous or discrete? */
    int *cptr = INTEGER(continuous); /* integer; continuous or discrete dist? */
    bool cont = true;
    if (cptr[0] == 0.) { cont = false; }

    /* Check if 'p' and/or 'q' are matrices.
     * They can be either a vector of the same length as 'x' (vector length n)
     * or a matrix of dimension (K x n) where the rows correspond to the number
     * of observations; the columns to the different quantiles.
     *
     * In the 'continuous' mode we typically have:
     *    - different quantiles for all observations (q is a matrix)
     *    - same probabilities for all observations (p is a vector)
     * In the 'discrete' mode we typically have:
     *    - same quantiles for all observations (q is a vector)
     *    - different probabilities for all observations (p is a matrix)
     * However, the function could also deal with p and q being a matrix.
     */
    bool ismatrix_p = false;
    bool ismatrix_q = false;
    if (ncols(p) > 1) { ismatrix_p = true; }
    if (ncols(q) > 1) { ismatrix_q = true; }

    /* Dimension of the matrix (p or q) */
    if (ismatrix_p) { J = ncols(p); } else if (ismatrix_q) { J = ncols(q); } else { J = nrows(p); }

    /* Initialize return vector with 0s */
    for (i = 0; i < n; i++) { rvalptr[i] = 0.0; }

    /* Indices description
     * i     {0, n}, i'th observation to be processed
     * ni    0 or i, index for xptr. If we only have one observation, ix = 0, else ix = i
     * j     {1, J}, column index of matrix p/q if they are matrices.
     *
     * Looping over all observations first */
    for (i = 0; i < n; i++) {

        if (nx == 1) { ix = 0; } else { ix = i; }

        /* Looping along the grid (along probabilities/quantiles) defined by j = {1, J}.
         * We don't start at 0 as we always need to include lower/upper bound {j-1, j}. */
        for (j = 1; j < J; j++) {

            /* Calculating position of lower and upper element along grid for p */
            if (ismatrix_p) {
                pi_lo = (j - 1) * n + i;       pi_up = j * n + i;
            } else {
                pi_lo = j - 1;                 pi_up = j;
            }

            /* When calculating the CRPS on 'Empirical' distributions the
             * probability can be NA in case the empirical distribution contains
             * missing values. As we pre-sort quantiles and probabilities, we only
             * have to check p_up (probability at the upper end of the current bin).
             */
            if (ISNA(pptr[pi_up])) { continue; }

            /* Calculating position of lower and upper element along grid for p */
            if (ismatrix_q) {
                qi_lo = (j - 1) * n + i;       qi_up = j * n + i;
            } else {
                qi_lo = j - 1;                 qi_up = j;
            }

            /* upper end of the 'bin' is lower than the observation */
            if (qptr[qi_up] <= xptr[ix]) {
                w     = qptr[qi_up] - qptr[qi_lo];
                if (cont) {
                    h_up  = pptr[pi_up];
                    h_lo  = pptr[pi_lo];
                    h     = (h_up * h_up + h_lo * h_lo) * 0.5;
                } else {
                    h     = pptr[pi_lo] * pptr[pi_lo];
                }
                tmp    = w * h;
            /* lower end of the 'bin' is larger than the observation */
            } else if (qptr[qi_lo] >= xptr[ix]) {
                w     = qptr[qi_up] - qptr[qi_lo];
                if (cont) {
                    h_up  = (1. - pptr[pi_up]);
                    h_lo  = (1. - pptr[pi_lo]);
                    h     = (h_up * h_up + h_lo * h_lo) * 0.5;
                } else {
                    h     = (1. - pptr[pi_lo]);
                    h     = h * h;
                }
                tmp    = w * h;
            /* the 'bin' contains the observation itself */
            } else {

                /* Lower */
                w      = xptr[ix] - qptr[qi_lo];
                if (cont) {
                    h_up   = pxptr[i];
                    h_lo   = pptr[pi_lo];
                    h      = (h_up * h_up + h_lo * h_lo) * 0.5;
                } else {
                    h      = pptr[pi_lo] * pptr[pi_lo];
                }
                tmp_lo = w * h;

                /* Upper */
                w      = qptr[qi_up] - xptr[ix];
                if (cont) {
                    h_up   = (1. - pptr[pi_up]);
                    h_lo   = (1. - pxptr[i]);
                    h      = (h_up * h_up + h_lo * h_lo) * 0.5;
                } else {
                    h      = (1. - pxptr[ix]);
                    /*Rprintf("   h = %.5f    (1 - %.5f) for ix = %d\n", h, pxptr[ix], ix);*/
                    h      = h * h;
                }
                tmp_up = w * h;

                /* Total contribution */
                tmp    = tmp_lo + tmp_up;
            }
            rvalptr[i] = rvalptr[i] + tmp;
        } /* end of loop over K */

        /* Updating 'outer bounds' if needed (if the observation lies outside the quantiles drawn) */
        /* If observation is < than largest quantile. We therefore first need to define the index
         * of the lowest quantile and the highest quantile of the current observation; re-use
         * qi_lo and qi_up for that. */
        if (ismatrix_q) {
            qi_lo = i; qi_up = n * (J - 1) + i;
        } else {
            qi_lo = 0; qi_up = (J - 1);
        }

        if (xptr[ix] < qptr[qi_lo]) {
            tmp   = qptr[qi_lo] - xptr[ix];
        /* If observation is > than largest quantile */
        } else if (xptr[ix] > qptr[qi_up]) {
            tmp     = xptr[ix] - qptr[qi_up];
        /* Else no adjustment needed */
        } else {
            tmp = 0;
        }
        rvalptr[i] = rvalptr[i] + tmp;
    } /* end of loop over n */

    UNPROTECT(1);
    return rval;
}

