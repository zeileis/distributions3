#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <stdbool.h> // for boolean values (true/false)
#include <R_ext/Rdynload.h>

#include "distributions3.h"

SEXP c_CRPS_numeric(SEXP, SEXP, SEXP, SEXP, SEXP);

/* Helper functions used for 2d interpolation and numeric integration */
double interpolate_linear(double xlo, double xhi, double ylo, double yhi, double x);
double integrate_2d(double xlo, double xhi, double ylo, double yhi, double x);

/* Functions calculating numeric moments */
SEXP c_moments_numeric(SEXP p, SEXP q, SEXP dim, SEXP discrete, SEXP whatint);
double c_moments_calculate_trapezoidal(int i, double* p, double* q, int* dim, int what);
double c_moments_calculate_discrete(int i,    double* p, double* q, int* dim, int what);

/* dpq functions for the shash distribution */
SEXP c_pshash(SEXP n, SEXP q, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP lower_tail, SEXP log_p, SEXP ncores);

const R_CallMethodDef CallEntries[] = {
  {"c_CRPS_numeric", (DL_FUNC) &c_CRPS_numeric, 5},
  {"c_moments_numeric", (DL_FUNC) &c_moments_numeric, 5},
  {"c_pshash", (DL_FUNC) &c_pshash, 9},
  {NULL, NULL, 0}
};

void R_init_distributions3(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


