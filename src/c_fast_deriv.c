#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <string.h>


// Fast element lookup by name for standard named R lists
SEXP getListElement(SEXP list, const char *str) {
    if (list == R_NilValue || !isNewList(list)) {
        return R_NilValue;
    }

    SEXP names = getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue) {
        return R_NilValue;
    }

    int n = LENGTH(list);
    for (int i = 0; i < n; i++) {
        const char *name = CHAR(STRING_ELT(names, i));
        if (strcmp(name, str) == 0) {
            return VECTOR_ELT(list, i);
        }
    }

    return R_NilValue;
}

// Struct to hold vector pointers and indexing strides (0 for scalar/length 1, 1 for vector/length N)
typedef struct {
    const double *ptr;
    int stride;
} InputVector;

//// Helper to extract a double vector and compute stride
//static inline InputVector make_input_var(SEXP list, const char *name, int expected_N) {
//    InputVector var = {NULL, 0};
//    SEXP elt = PROTECT(getListElement(list, name)); // assumes standard R getListElement or manual loop
//    if (elt != R_NilValue && (isNumeric(elt) || isInteger(elt))) {
//        var.ptr = REAL(elt);
//        var.stride = (LENGTH(elt) == expected_N) ? 1 : 0;
//    }
//    UNPROTECT(1);
//    return var;
//}
static inline InputVector make_input_var(SEXP sexp, const char *name, int expected_N) {
    InputVector var = {NULL, 0};

    if (sexp == R_NilValue) {
        return var;
    }

    SEXP elt = sexp;

    // If a name is provided, extract the element from the list
    if (name != NULL) {
        elt = getListElement(sexp, name);
        if (elt == R_NilValue) return var;
    }

    // Process target SEXP object
    if (elt != R_NilValue && (isNumeric(elt) || isInteger(elt))) {
        var.ptr = REAL(elt);
        var.stride = (LENGTH(elt) == expected_N) ? 1 : 0;
    }

    return var;
}

/* Get Value from InputVec
 *
 * Helper to safely index InputVector (`var`). If the input
 * vector is of length 1, stride equals 0 (int) and the function
 * returns the first element (index 0), else the i-th element is
 * returned.
 */
static inline double get_val(InputVector *var, int i) {
    return var->ptr[i * var->stride];
}

/* Get Maximum Length
 *
 * Whlie x is a SEXP vector, params is a named list (SEXP)
 * with double vectors. Both 'x' as well as all elements in
 * 'params' is allowed to either be of length 1 or N.
 *
 * Returns the maximum length of 'x' and all vectors in
 * the list 'params'.
 */
static int get_max_N(SEXP x, SEXP params) {
    int max_N = LENGTH(x);

    int num_params = LENGTH(params);
    for (int p = 0; p < num_params; p++) {
        SEXP param = VECTOR_ELT(params, p);
        if (param != R_NilValue) {
            int len = LENGTH(param);
            if (len > max_N) {
                max_N = len;
            }
        }
    }

    return max_N;
}


/* Validate Vector Lengths
 *
 * Checks that all vectors, i.e., 'x' as well as all vectors
 * in the 'params' list, are either of length 1 or N (maximum
 * length over all vectors).
 *
 * Returns the maximum length of 'x' and all vectors in
 * the list 'params'.
 */
int validate_lengths(SEXP x, SEXP params) {
    int N = get_max_N(x, params);

    // Check x, just be of length 1 or N
    int len = LENGTH(x);
    if (len != 1 && len != N) {
        Rf_error("[C] Invalid length for 'x': expected 1 or %d, got %d.", N, len);
    }

    // Check elements of params list
    int num_params = LENGTH(params);
    SEXP names = getAttrib(params, R_NamesSymbol);

    // Check vectors in list, must be of length 1 or N
    for (int p = 0; p < num_params; p++) {
        SEXP elt = VECTOR_ELT(params, p);
        int  len = LENGTH(elt);
        if (len != 1 && len != N) {
            Rf_error("[C] Invalid length for parameter '%s': expected 1 or %d, got %d.",
                CHAR(STRING_ELT(names, p)), N, len);
        }
    }

    // Return N, maximum length of 'x' and vectors in 'params'
    return N;
}


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

    // First derivative/main parameters
    PARAM_MU    = 1U << 0,
    PARAM_SIGMA = 1U << 1,
    PARAM_NU    = 1U << 2,
    PARAM_TAU   = 1U << 3,
    // .. up to 1U << 6 (N = 7)

    // Cross-derivatives
    PARAM_MU_MU         = 1U << 7,
    PARAM_MU_SIGMA      = 1U << 8,
    PARAM_MU_TAU        = 1U << 9,
    PARAM_MU_NU         = 1U << 10,
    PARAM_SIGMA_SIGMA   = 1U << 11,
    PARAM_SIGMA_NU      = 1U << 12,
    PARAM_SIGMA_TAU     = 1U << 13,
    PARAM_NU_NU         = 1U << 14,
    PARAM_NU_TAU        = 1U << 15,
    PARAM_TAU_TAU       = 1U << 16
    // .. up to 1U << 27 (N = 7 plus all cross-derivatives)
} ParameterFlags;

/* Setting Required Parameter Flags
 *
 * @param x a named list where the names define which score or hessian
 *        is required.
 * @param hessian set false for score, and true for hessian.
 *
 * Returns an enum object used to check which derivatives are needed.
 */
static inline ParameterFlags get_score_flags(SEXP x, bool hessian) {
    ParameterFlags req = PARAM_NONE;

    // Score
    if (!hessian) {
        if (getListElement(x, "mu")    != R_NilValue) req |= PARAM_MU;
        if (getListElement(x, "sigma") != R_NilValue) req |= PARAM_SIGMA;
        if (getListElement(x, "nu")    != R_NilValue) req |= PARAM_NU;
        if (getListElement(x, "tau")   != R_NilValue) req |= PARAM_TAU;
    // Hessian
    } else {
        if (getListElement(x, "mu:mu")       != R_NilValue) req |= (PARAM_MU | PARAM_MU_MU);
        if (getListElement(x, "mu:sigma")    != R_NilValue) req |= (PARAM_MU | PARAM_SIGMA | PARAM_MU_SIGMA);
        if (getListElement(x, "mu:tau")      != R_NilValue) req |= (PARAM_MU | PARAM_TAU | PARAM_MU_TAU);
        if (getListElement(x, "mu:nu")       != R_NilValue) req |= (PARAM_MU | PARAM_NU | PARAM_MU_NU);
        if (getListElement(x, "sigma:sigma") != R_NilValue) req |= (PARAM_SIGMA | PARAM_SIGMA_SIGMA);
        if (getListElement(x, "sigma:nu")    != R_NilValue) req |= (PARAM_SIGMA | PARAM_NU | PARAM_SIGMA_NU);
        if (getListElement(x, "sigma:tau")   != R_NilValue) req |= (PARAM_SIGMA | PARAM_TAU | PARAM_SIGMA_TAU);
        if (getListElement(x, "nu:nu")       != R_NilValue) req |= (PARAM_NU | PARAM_NU_NU);
        if (getListElement(x, "nu:tau")      != R_NilValue) req |= (PARAM_NU | PARAM_TAU | PARAM_NU_TAU);
        if (getListElement(x, "tau:tau")     != R_NilValue) req |= (PARAM_TAU | PARAM_TAU_TAU);
    }

    return req;
}

SEXP c_fast_derivatives(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp) {

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // 1. Setup Input Parameters (Pointers + Strides)
    //    Add additional parameters as needed
    InputVector x_var     = make_input_var(x_sexp, NULL, N);
    InputVector mu_var    = make_input_var(params_sexp, "mu",    N);
    InputVector sigma_var = make_input_var(params_sexp, "sigma", N);
    InputVector nu_var    = make_input_var(params_sexp, "nu",    N);

    SEXP s_mu    = getListElement(score_sexp, "mu");
    SEXP s_sigma = getListElement(score_sexp, "sigma");
    return R_NilValue;
    SEXP s_nu    = getListElement(score_sexp, "nu");
    SEXP s_tau   = getListElement(score_sexp, "tau");

    int i;
    for (i = 0; i < N; i++) {
        Rprintf(" ---- i = %d\n", i);
        Rprintf("      mu.stride = %d\n", mu_var.stride);
        Rprintf("      mu[%d] = %.5f\n", i, get_val(&mu_var, i));
        Rprintf("      x.stride = %d\n", x_var.stride);
        Rprintf("      x[%d] = %.5f\n", i, get_val(&x_var, i));
        Rprintf("      sigma.stride = %d\n", sigma_var.stride);
        Rprintf("      sigma[%d] = %.5f\n", i, get_val(&sigma_var, i));
        s_mu[i] = i;
        s_tau[i] = i / 10;
    }

    // 2. Setup Score Outputs & Active Flags
    ParameterFlags req_scores  = get_score_flags(score_sexp, false);
    ParameterFlags req_hessian = get_score_flags(hessian_sexp, true);


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
////
////    // 4. Main Execution Loop
////    for (int i = 0; i < N; i++) {
////        double x_i = get_val(&x_var, i);
////        
////        // Extract inputs cleanly via stride multiply
////        double mu_i    = get_val(&mu_var, i);
////        double sigma_i = get_val(&sigma_var, i);
////        double nu_i    = get_val(&nu_var, i);
////
////        // Pre-allocate intermediate scalar derivatives
////        double d_mu = 0.0, d_sigma = 0.0, d_nu = 0.0;
////
////        // Compute base derivatives conditionally
////        if (calc_mu) {
////            d_mu = (x_i - mu_i) / (sigma_i * sigma_i); // Example intermediate
////        }
////        if (calc_sigma) {
////            d_sigma = -1.0 / sigma_i + ((x_i - mu_i) * (x_i - mu_i)) / (sigma_i * sigma_i * sigma_i);
////        }
////        if (calc_nu) {
////            d_nu = nu_i * x_i; // Example intermediate
////        }
////
////        // Fill Scores
////        if (s_mu)    s_mu[i]    = d_mu;
////        if (s_sigma) s_sigma[i] = d_sigma;
////        if (s_nu)    s_nu[i]    = d_nu;
////
////        // Fill Hessians
////        if (h_mu_sigma) h_mu_sigma[i] = -2.0 * (x_i - mu_i) / (sigma_i * sigma_i * sigma_i);
////        if (h_nu_sigma) h_nu_sigma[i] = d_nu * d_sigma;
////    }
////
////    return R_NilValue;
}
