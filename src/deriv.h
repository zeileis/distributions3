// Ensure the header file is only included once during compilation
#pragma once

#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <string.h>

// Fast element lookup by name for standard named R lists
static inline SEXP getListElementSEXP(SEXP list, const char *str) {
    if (list == R_NilValue || !isNewList(list)) return R_NilValue;

    SEXP names = getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue) return R_NilValue;

    int n = LENGTH(list);
    for (int i = 0; i < n; i++) {
        SEXP name_elt = STRING_ELT(names, i);
        if (name_elt != R_NilValue && strcmp(CHAR(name_elt), str) == 0) {
            return VECTOR_ELT(list, i);
        }
    }

    return R_NilValue;
}

static inline double* getListElement(SEXP list, const char *str) {
    SEXP elt = getListElementSEXP(list, str);
    if (elt != R_NilValue && TYPEOF(elt) == REALSXP) return REAL(elt);
    return NULL;
}


// Struct to hold vector pointers and indexing strides (0 for scalar/length 1, 1 for vector/length N)
typedef struct {
    const double *ptr;
    int stride;
} InputVector;

/*
 * Create Input Vector
 *
 * Helper function/type to create easy-to-use input vectors and strides.
 * The strides are used to allow all vectors to either be of length 1 or N.
 * If the length of the vector is 1, the stride equals 0 (int) used to always
 * access the 0-th element in the for loops. Else (length is equal to the
 * 'expected_N') the stride is set 1 (int).
 *
 * Params
 * ------
 * sexp:        a vector (typically for x, q, p) or a named list of double vectors
 *              containing the parameters of the distribution.
 * name:        Null (if sexp is a vector) or name of the list element to be returned.
 * expected_N:  The expected length if the double vector is not of length 1.
 *
 * Return
 * ------
 * Returns an object of class InputVector which consists of a double pointer (.ptr)
 * and an integer (.stride). Function `get_val(&x, i)` is used to extract the i-th
 * element of the input vector which returns element 0 (if the vector is of length 1)
 * or the i-th element if the length is N.
 */
static inline InputVector make_input_var(SEXP sexp, const char *name, int expected_N) {
    InputVector var = {NULL, 0};

    if (sexp == R_NilValue) return var;

    // The sexp vector is our result if name == NULL (sexp is already vector)
    SEXP elt = sexp;

    // If a name is provided, extract the element from the list
    if (name != NULL) {
        elt = getListElementSEXP(sexp, name);
        if (elt == R_NilValue) return var;
    }

    // Process target SEXP object
    if (elt != R_NilValue && (isNumeric(elt) || isInteger(elt))) {
        var.ptr    = REAL(elt);
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
static inline ParameterFlags get_flags_msnt(SEXP x, bool hessian) {
    ParameterFlags req = PARAM_NONE;

    // Score
    if (!hessian) {
        if (getListElementSEXP(x, "mu")    != R_NilValue) req |= PARAM_MU;
        if (getListElementSEXP(x, "sigma") != R_NilValue) req |= PARAM_SIGMA;
        if (getListElementSEXP(x, "nu")    != R_NilValue) req |= PARAM_NU;
        if (getListElementSEXP(x, "tau")   != R_NilValue) req |= PARAM_TAU;
    // Hessian
    } else {
        if (getListElementSEXP(x, "mu:mu")       != R_NilValue) req |= (PARAM_MU | PARAM_MU_MU);
        if (getListElementSEXP(x, "mu:sigma")    != R_NilValue) req |= (PARAM_MU | PARAM_SIGMA | PARAM_MU_SIGMA);
        if (getListElementSEXP(x, "mu:tau")      != R_NilValue) req |= (PARAM_MU | PARAM_TAU | PARAM_MU_TAU);
        if (getListElementSEXP(x, "mu:nu")       != R_NilValue) req |= (PARAM_MU | PARAM_NU | PARAM_MU_NU);
        if (getListElementSEXP(x, "sigma:sigma") != R_NilValue) req |= (PARAM_SIGMA | PARAM_SIGMA_SIGMA);
        if (getListElementSEXP(x, "sigma:nu")    != R_NilValue) req |= (PARAM_SIGMA | PARAM_NU | PARAM_SIGMA_NU);
        if (getListElementSEXP(x, "sigma:tau")   != R_NilValue) req |= (PARAM_SIGMA | PARAM_TAU | PARAM_SIGMA_TAU);
        if (getListElementSEXP(x, "nu:nu")       != R_NilValue) req |= (PARAM_NU | PARAM_NU_NU);
        if (getListElementSEXP(x, "nu:tau")      != R_NilValue) req |= (PARAM_NU | PARAM_TAU | PARAM_NU_TAU);
        if (getListElementSEXP(x, "tau:tau")     != R_NilValue) req |= (PARAM_TAU | PARAM_TAU_TAU);
    }

    return req;
}

