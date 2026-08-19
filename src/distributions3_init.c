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

/* dpq functions for the sinharcsinh distribution */
SEXP c_psinharcsinh(SEXP N, SEXP q, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP lower_tail, SEXP log_p, SEXP ncores);
SEXP c_dsinharcsinh(SEXP N, SEXP x, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP ret_log, SEXP ncores);
SEXP c_qsinharcsinh(SEXP N, SEXP p, SEXP mu, SEXP sigma, SEXP nu, SEXP tau, SEXP lower_tail, SEXP log_p, SEXP cores);
double local_zeroin(double ax, double bx, double (*f)(double x, void *info), void *info, double tol);

/* dpq for zero adjusted negative binomial */
SEXP c_dZeroAdjustedNegbin(SEXP N, SEXP x, SEXP mu, SEXP sigma, SEXP nu, SEXP ret_log, SEXP ncores);
SEXP c_pZeroAdjustedNegbin(SEXP N, SEXP q, SEXP mu, SEXP sigma, SEXP nu, SEXP lower_tail, SEXP log_p, SEXP ncores);
SEXP c_qZeroAdjustedNegbin(SEXP N, SEXP p, SEXP mu, SEXP sigma, SEXP nu, SEXP lower_tail, SEXP log_p, SEXP ncores);

/* deriv2: combines score and hessian function */
SEXP c_deriv_SinhArcsinh(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP ncores);
SEXP c_deriv_Normal(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP expected, SEXP ncores);

static R_CallMethodDef CallEntries[] = {
  {"c_CRPS_numeric", (DL_FUNC) &c_CRPS_numeric, 5},
  {"c_moments_numeric", (DL_FUNC) &c_moments_numeric, 5},

  {"c_psinharcsinh", (DL_FUNC) &c_psinharcsinh, 9},
  {"c_dsinharcsinh", (DL_FUNC) &c_dsinharcsinh, 8},
  {"c_qsinharcsinh", (DL_FUNC) &c_qsinharcsinh, 9},

  {"c_dZeroAdjustedNegbin", (DL_FUNC) &c_dZeroAdjustedNegbin, 7},
  {"c_pZeroAdjustedNegbin", (DL_FUNC) &c_pZeroAdjustedNegbin, 8},
  {"c_qZeroAdjustedNegbin", (DL_FUNC) &c_qZeroAdjustedNegbin, 8},

  {"c_deriv_SinhArcsinh", (DL_FUNC) &c_deriv_SinhArcsinh, 5},
  {"c_deriv_Normal", (DL_FUNC) &c_deriv_Normal, 6},
  {NULL, NULL, 0}
};

void R_init_distributions3(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


