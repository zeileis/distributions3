# Method for numerically evaluating the CRPS of probability distributions

Method to the [`crps`](https://rdrr.io/pkg/scoringRules/man/scores.html)
generic function from the scoringRules package for numerically
evaluating the (continuous) ranked probability score (CRPS) of any
probability distributions3 object.

## Usage

``` r
crps.distribution(
  y,
  x,
  drop = TRUE,
  elementwise = NULL,
  gridsize = 500L,
  batchsize = 10000L,
  applyfun = NULL,
  cores = NULL,
  method = NULL,
  ...
)

crps.Beta(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Bernoulli(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Binomial(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Erlang(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Exponential(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Gamma(y, x, drop = TRUE, elementwise = NULL, ...)

crps.GEV(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Geometric(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Gumbel(y, x, drop = TRUE, elementwise = NULL, ...)

crps.HyperGeometric(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Logistic(y, x, drop = TRUE, elementwise = NULL, ...)

crps.LogNormal(y, x, drop = TRUE, elementwise = NULL, ...)

crps.NegativeBinomial(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Normal(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Poisson(y, x, drop = TRUE, elementwise = NULL, ...)

crps.StudentsT(y, x, drop = TRUE, elementwise = NULL, ...)

crps.Uniform(y, x, drop = TRUE, elementwise = NULL, ...)

crps.XBetaX(y, x, drop = TRUE, elementwise = NULL, method = "cdf", ...)

crps.GAMLSS(y, x, drop = TRUE, elementwise = NULL, ...)

crps.BAMLSS(y, x, drop = TRUE, elementwise = NULL, ...)
```

## Arguments

- y:

  A distribution object, e.g., as created by
  [`Normal`](https://zeileis.github.io/distributions3/dev/reference/Normal.md)
  or
  [`Binomial`](https://zeileis.github.io/distributions3/dev/reference/Binomial.md).

- x:

  A vector of elements whose CRPS should be determined given the
  distribution `y`.

- drop:

  logical. Should the result be simplified to a vector if possible?

- elementwise:

  logical. Should each distribution in `y` be evaluated at all elements
  of `x` (`elementwise = FALSE`, yielding a matrix)? Or, if `y` and `x`
  have the same length, should the evaluation be done element by element
  (`elementwise = TRUE`, yielding a vector)? The default of `NULL` means
  that `elementwise = TRUE` is used if the lengths match and otherwise
  `elementwise = FALSE` is used.

- gridsize:

  positive integer. Size of the grid used to approximate the CDF for the
  numerical calculation of the CRPS.

- batchsize:

  positive integer. Maximum batch size. Used to split the input into
  batches. Lower values reduce required memory but may increase
  computation time.

- applyfun:

  An optional [`lapply`](https://rdrr.io/r/base/lapply.html)-style
  function with arguments `function(X, FUN, ...)`. It is used to compute
  the CRPS for each element of `y`. The default is to use the basic
  `lapply` function unless the `cores` argument is specified (see
  below).

- cores:

  `NULL` (default) or positive integer. If integer, the `applyfun` is
  set to [`mclapply`](https://rdrr.io/r/parallel/mclapply.html) with the
  desired number of `cores`, except on Windows where
  [`parLapply`](https://rdrr.io/r/parallel/clusterApply.html) with
  `makeCluster(cores)` is used.

- method:

  `NULL` (default) or character. Specifies how the grid for the
  calculation is set up. If `NULL` it is set to `method = "cdf"` if the
  distribution (`y`) is discrete and `gridsize` is large enough to span
  the required range for CRPS calculation. Else (including continuous
  distributions) `method = "quantile"` is used.

- ...:

  Currently not used.

## Value

If `length(y)` equals one or `length(y) = length(x)` a (possibly named)
numeric vector is returned if `drop = TRUE` (default). Else a matrix of
dimension `length(y)` times `length(x)` is returned containing the CRPS.

## Details

The (continuous) ranked probability score (CRPS) for (univariate)
probability distributions can be computed based on the the
object-oriented infrastructure provided by the distributions3 package.
The general `crps.distribution` method does so by using numeric
integration based on the `cdf` and/or `quantile` methods (for more
details see below). If dedicated closed-form CRPS computations are
provided (either by the scoringRules package or distributions3) these
are used as they are both computationally faster and numerically more
precise. For example, the `crps` method for `Normal` objects leverages
[`crps_norm`](https://rdrr.io/pkg/scoringRules/man/scores_norm.html)
rather than relying on numeric integration.

## Numerical approximation strategy

The general method for any `distribution` object uses the following
strategy. By default (if the `method` argument is `NULL`), it
distinguishes distributions whose entire support is continuous, or whose
entire support is discrete or mixed discrete-continuous distribution
using
[`is_continuous`](https://zeileis.github.io/distributions3/dev/reference/is_discrete.md)
and
[`is_discrete`](https://zeileis.github.io/distributions3/dev/reference/is_discrete.md),
respectively.

For continuous and mixed distributions, an equidistant grid of
`gridsize + 5L` probabilities is drawn for which the corresponding
`quantile`s for each distribution `y` are calculated (including the
observation `x`). The calculation of the CRPS then uses a trapezoidal
approximation for the numeric integration. For discrete distributions,
`gridsize` equidistant quantiles (in steps of `1L`) are drawn and the
corresponding probabilities from the `cdf` are calculated (including the
observation `x`) used to numerically integrate for calculating the CRPS.
If `gridsize` is not sufficient to cover the required range, the method
falls back to the procedure used for continuous and mixed distributions
to approximate the CRPS.

If the `method` argument is set to either `"cdf"` or `"quantile"`, then
the specific strategy for setting up the grid of observations and
corresponding probabilities can be enforced. This can be useful if for a
certain distribution class, only a `cdf` or only a `quantile` method is
available or only one of them is numerically stable or computationally
efficient etc.

The numeric approximation requires to set up a matrix of dimension
`length(y) * (gridsize + 5L)` (or `length(y) * (gridsize + 1L)`) which
may be very memory intensive if `length(y)` and/or `gridsize` are large.
Thus, the data is split in batches of (approximately) equal size, not
larger than `batchsize`. Thus, the memory requirement is reduced to
`batchsize * (gridsize + 5)` in each step. Hence, a smaller value of
`batchsize` will reduce memory footprint but will slightly increase
computation time.

The error (deviation between numerical approximation and analytic
solution) has been shown to be in the order of `1e-2` for a series of
distributions tested. Accuracy can be increased by increasing `gridsize`
and will be lower for a smaller `gridsize`.

For parallelization of the numeric computations, a suitable `applyfun`
can be provided that carries out the integration for each element of
`y`. To facilitate setting up a suitable `applyfun` using the basic
parallel package, the argument `cores` is provided for convenience. When
used, `y` is split into `B` equidistant batches; at least `B = cores`
batches or a multiple of `cores` with a maximum size of `batchsize`. On
systems running Windows `parlapply` is used, else `mclapply`.

## Examples

``` r
set.seed(6020)

## three normal distributions X and observations x
library("distributions3")
X <- Normal(mu = c(0, 1, 2), sigma = c(2, 1, 1))
x <- c(0, 0, 1)

## evaluate crps
## using infrastructure from scoringRules (based on closed-form analytic equations)
library("scoringRules")
crps(X, x)
#> [1] 0.4673900 0.6024414 0.6024414

## using general distribution method explicitly (based on numeric integration)
crps.distribution(X, x)
#> [1] 0.4674058 0.6024482 0.6024482

## analogously for Poisson distribution
Y <- Poisson(c(0.5, 1, 2))
crps(Y, x)
#> [1] 0.1631650 0.4762224 0.4991650
crps.distribution(Y, x)
#> [1] 0.1631650 0.4762224 0.4991650
```
