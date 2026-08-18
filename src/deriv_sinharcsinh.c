#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <string.h>
#include "deriv.h"

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

    // First derivative/main parameters (used for score and hessian)
    PARAM_MU    = 1U << 0,
    PARAM_SIGMA = 1U << 1,
    PARAM_NU    = 1U << 2,
    PARAM_TAU   = 1U << 3,

    // Cross-derivatives (used for hessian)
    PARAM_MU_MU         = 1U << 4,
    PARAM_MU_SIGMA      = 1U << 5,
    PARAM_MU_TAU        = 1U << 6,
    PARAM_MU_NU         = 1U << 7,
    PARAM_SIGMA_SIGMA   = 1U << 8,
    PARAM_SIGMA_NU      = 1U << 9,
    PARAM_SIGMA_TAU     = 1U << 10,
    PARAM_NU_NU         = 1U << 11,
    PARAM_NU_TAU        = 1U << 12,
    PARAM_TAU_TAU       = 1U << 13
    // .. up to 1U << 31
} ParameterFlags;

/* Setting Required Parameter Flags
 *
 * @param x a named list where the names define which score or hessian
 *        is required.
 * @param hessian set false for score, and true for hessian.
 *
 * Returns an enum object used to check which derivatives are needed.
 */
static inline ParameterFlags get_flags(SEXP x, bool hessian) {
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


SEXP c_deriv_sinharcsinh(SEXP x_sexp, SEXP params_sexp, SEXP score_sexp, SEXP hessian_sexp, SEXP ncores) {

    #if _OPENMP
    int nthreads = asInteger(ncores); // Number of threads
    #endif

    // Validate vector lengths, returns maximum length
    int N = validate_lengths(x_sexp, params_sexp);

    // Setup Input Parameters (pointers + strides) for all parameters of the distribution.
    InputVector x_vec     = make_input_var(x_sexp, NULL, N);
    InputVector mu_vec    = make_input_var(params_sexp, "mu",    N);
    InputVector sigma_vec = make_input_var(params_sexp, "sigma", N);
    InputVector nu_vec    = make_input_var(params_sexp, "nu",    N);
    InputVector tau_vec   = make_input_var(params_sexp, "tau",   N);

    // Setting up flags for score and hessian; used for fast bitwise operations
    // to see what needs to be calculated. 'req_scores' means 'required scores',
    // analogously 'req_hessian' is used to check what is required to get the
    // requested hessian elements.
    ParameterFlags req_scores  = get_flags(score_sexp, false);
    ParameterFlags req_hessian = get_flags(hessian_sexp, true);


    // If score for mu, sigma, nu, or tau is requested: Initialize double
    // pointer for the corresponding `score_sexp` list element to be modified
    // and returned to R.
    double* s_mu        = NULL;
    double* s_sigma     = NULL;
    double* s_tau       = NULL;
    double* s_nu        = NULL;
    if (req_scores & PARAM_MU)           s_mu    = getListElement(score_sexp, "mu");
    if (req_scores & PARAM_SIGMA)        s_sigma = getListElement(score_sexp, "sigma");
    if (req_scores & PARAM_NU)           s_nu    = getListElement(score_sexp, "nu");
    if (req_scores & PARAM_TAU)          s_tau   = getListElement(score_sexp, "tau");

    // As for score: Initialize double pointer for required elements corresponding
    // to the elements in hessian_sexp to be modified and returned to R.
    double* h_mu_mu     = NULL;
    double* h_mu_sigma  = NULL;
    double* h_mu_nu     = NULL;
    double* h_mu_tau    = NULL;
    if (req_hessian & PARAM_MU_MU)       h_mu_mu    = getListElement(hessian_sexp, "mu:mu");
    if (req_hessian & PARAM_MU_SIGMA)    h_mu_sigma = getListElement(hessian_sexp, "mu:sigma");
    if (req_hessian & PARAM_MU_NU)       h_mu_nu    = getListElement(hessian_sexp, "mu:nu");
    if (req_hessian & PARAM_MU_TAU)      h_mu_tau   = getListElement(hessian_sexp, "mu:tau");

    double* h_sigma_sigma = NULL;
    double* h_sigma_nu    = NULL;
    double* h_sigma_tau   = NULL;
    if (req_hessian & PARAM_SIGMA_SIGMA) h_sigma_sigma    = getListElement(hessian_sexp, "sigma:sigma");
    if (req_hessian & PARAM_SIGMA_NU)    h_sigma_nu       = getListElement(hessian_sexp, "sigma:nu");
    if (req_hessian & PARAM_SIGMA_TAU)   h_sigma_tau      = getListElement(hessian_sexp, "sigma:tau");

    double* h_nu_nu       = NULL;
    double* h_nu_tau      = NULL;
    if (req_hessian & PARAM_NU_NU)       h_nu_nu          = getListElement(hessian_sexp, "nu:nu");
    if (req_hessian & PARAM_NU_TAU)      h_nu_tau         = getListElement(hessian_sexp, "nu:tau");

    double* h_tau_tau     = NULL;
    if (req_hessian & PARAM_TAU_TAU)     h_tau_tau        = getListElement(hessian_sexp, "tau:tau");


    // Pre-calculated boolean flags which indicate which elements to be calculated
    // calc_mu is set 'true' if either the score requires dl/dmu or the calculations
    // of the hessian elements require dl/dmu. Same for sigma, nu, tau.
    bool calc_mu    = (req_scores & PARAM_MU)    || (req_hessian & PARAM_MU);
    bool calc_sigma = (req_scores & PARAM_SIGMA) || (req_hessian & PARAM_SIGMA);
    bool calc_nu    = (req_scores & PARAM_NU)    || (req_hessian & PARAM_NU);
    bool calc_tau   = (req_scores & PARAM_TAU)   || (req_hessian & PARAM_TAU);

    // Main loop
    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int i = 0; i < N; i++) {
        // Extracting numeric elements
        double x     = get_val(&x_vec, i);
        double mu    = get_val(&mu_vec, i);
        double sigma = get_val(&sigma_vec, i);
        double nu    = get_val(&nu_vec, i);
        double tau   = get_val(&tau_vec, i);

        // Calculating elements required multiple times for
        // the calculation of the derivatives (score/hessian)
        double z                 = (x - mu) / sigma;
        double z2                = z * z;
        double tau2              = tau * tau;
        double nu2               = nu * nu;

        //double asinhz            = asinh(z);
        double z2p1sqrt          = hypot(1.0, z); //sqrt(z2 + 1.0);
        double z2p1sqrtinv       = 1.0 / z2p1sqrt;
        double asinhz            = log(z + z2p1sqrt); // faster than asinh(z)
        double exp_tauasinhz     = exp( tau * asinhz);
        double exp_minusnuasinhz = exp(-nu  * asinhz);

        // Performing a series of vector operations used multiple times below
        double sigmainv = 1.0 / sigma;

        double r = 0.5 * (exp_tauasinhz        - exp_minusnuasinhz);
        double c = 0.5 * (exp_tauasinhz * tau  + exp_minusnuasinhz * nu);
        double h = 0.5 * (exp_tauasinhz * tau2 - exp_minusnuasinhz * nu2);

        // Partial derivatives used everywhere
        double dldr = -r;
        double dldc = 1 / c;

        // -------------- shared btw. score/hessian -------------

        // Calculate scores (dldm, dldd, dlds, dldm) in case they
        // are requested as scores, or are required to calculate the hessian
        // elements requested by the user.
        double dldz, dcdz, drdz, dzdm, dldm;
        if (calc_mu) {
            dldz = -z / (1.0 + z2);
            dcdz = h * z2p1sqrtinv;
            drdz = c * z2p1sqrtinv;
            dzdm = -sigmainv;
            dldm = (dldr * drdz + dldc * dcdz + dldz) * dzdm;
        }

        double dzdd, dldd;
        if (calc_sigma) {
            dzdd = -z * sigmainv;
            dldd = (dldr * drdz + dldc * dcdz + dldz) * dzdd - sigmainv;
        }

        double drdv, dcdv, dldv;
        if (calc_nu) {
            drdv = 0.5 * asinhz * exp_minusnuasinhz;
            dcdv = 0.5 * (1.0 - nu * asinhz) * exp_minusnuasinhz;
            dldv = dldr * drdv + dldc * dcdv;
        }

        double drdt, dcdt, dldt;
        if (calc_tau) {
            drdt = 0.5 * asinhz * exp_tauasinhz;
            dcdt = 0.5 * (1.0 + tau * asinhz) * exp_tauasinhz;
            dldt = dldr * drdt + dldc * dcdt;
        }

        // ------------------- store score ----------------------

        // Score 'mu' (corresponds to gamlss.dist::SHASH()$dldm
        if (req_scores & PARAM_MU)      s_mu[i]    = dldm;
        // Score 'sigma' (corresponds to gamlss.dist::SHASH()$dldd
        if (req_scores & PARAM_SIGMA)   s_sigma[i] = dldd;
        // Score 'nu' (corresponds to gamlss.dist::SHASH()$dldv
        if (req_scores & PARAM_NU)      s_nu[i]    = dldv;
        // Score 'tau' (corresponds to gamlss.dist::SHASH()$dldt
        if (req_scores & PARAM_TAU)     s_tau[i]   = dldt;

        // ------------- calculate and store hessian ------------

        // Hessian 'mu:mu' (corresponds to gamlss.dist::SHASH()$d2ldm2)
        if (req_hessian & PARAM_MU_MU)         h_mu_mu[i]         = pmin_double(-dldm * dldm, -1e-15);
        // Hessian 'mu:sigma' (corresponds to gamlss.dist::SHASH()$d2ldmdd)
        if (req_hessian & PARAM_MU_SIGMA)      h_mu_sigma[i]      = -dldm * dldd;
        // Hessian 'mu:tau' (corresponds to gamlss.dist::SHASH()$d2ldmdt)
        if (req_hessian & PARAM_MU_TAU)        h_mu_tau[i]        = -dldm * dldt;
        // Hessian 'mu:nu' (corresponds to gamlss.dist::SHASH()$d2ldmdv)
        if (req_hessian & PARAM_MU_NU)         h_mu_nu[i]         = -dldm * dldv;
        // Hessian 'sigma:sigma' (corresponds to gamlss.dist::SHASH()$d2ldd2)
        if (req_hessian & PARAM_SIGMA_SIGMA)   h_sigma_sigma[i]   = pmin_double(-dldd * dldd, -1e-15);
        // Hessian 'sigma:nu' (corresponds to gamlss.dist::SHASH()$d2ldddv)
        if (req_hessian & PARAM_SIGMA_NU)      h_sigma_nu[i]      = -dldd * dldv;
        // Hessian 'sigma:tau' (corresponds to gamlss.dist::SHASH()$d2ldddt)
        if (req_hessian & PARAM_SIGMA_TAU)     h_sigma_tau[i]     = -dldd * dldt;
        // Hessian 'nu:nu' (corresponds to gamlss.dist::SHASH()$d2ldv2)
        if (req_hessian & PARAM_NU_NU)         h_nu_nu[i]         = pmin_double(-dldv * dldv, -1e-15);
        // Hessian 'nu:tau' (corresponds to gamlss.dist::SHASH()$d2ldvdt)
        if (req_hessian & PARAM_NU_TAU)        h_nu_tau[i]        = -dldv * dldt;
        // Hessian 'tau:tau' (corresponds to gamlss.dist::SHASH()$d2ldt2)
        if (req_hessian & PARAM_TAU_TAU)       h_tau_tau[i]       = pmin_double(-dldt * dldt, -1e-15);

    }

    return R_NilValue;
}
