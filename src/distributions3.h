#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* Custom type: stuctured object with a double vector and length.
 * NOTE: .values must be freed by the user! */
typedef struct {
    double* values;
    int length;
} doubleVec;
