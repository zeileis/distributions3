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


/* Mininum: Helper function mimiking pmin(a, b) in R for single doubles! */
static inline double pmin_double(double a, double b) {
    return (a < b) ? a : b;
}
/* Maximum: Helper function mimiking pmin(a, b) in R for single doubles! */
static inline double pmax_double(double a, double b) {
    return (a > b) ? a : b;
}
