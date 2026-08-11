# Simulate responses from fitted model objects

Default method for simulating new responses from any model object with a
[`prodist`](https://zeileis.github.io/distributions3/dev/reference/prodist.md)
method (for extracting a probability distribution object).

## Usage

``` r
# Default S3 method
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  An object for which a
  [`prodist`](https://zeileis.github.io/distributions3/dev/reference/prodist.md)
  method is available.

- nsim:

  The number of response vectors to simulate. Should be a positive
  integer. Defaults to 1.

- seed:

  An optional random seed that is to be set using
  [`set.seed`](https://rdrr.io/r/base/Random.html) prior to drawing the
  random sample. The previous random seed from the global environment
  (if any) is restored afterwards.

- ...:

  Arguments passed to
  [`simulate.distribution`](https://zeileis.github.io/distributions3/dev/reference/random.md).

## Value

A data frame with an attribute `"seed"` containing the `.Random.seed`
from before the simulation.

## Details

This default method simply combines two building blocks provided in this
package: (1)
[`prodist`](https://zeileis.github.io/distributions3/dev/reference/prodist.md)
for extracting the probability distribution from a fitted model object,
(2)
[`simulate.distribution`](https://zeileis.github.io/distributions3/dev/reference/random.md)
for simulating new observations from this distribution (internally
calling
[`random`](https://zeileis.github.io/distributions3/dev/reference/random.md)).

Thus, this enables simulation from any fitted model object that provides
a `prodist` method. It waives the need to implement a dedicated
[`simulate`](https://rdrr.io/r/stats/simulate.html) method for this
model class.

## Examples

``` r
## Poisson GLM for FIFA 2018 goals
data("FIFA2018", package = "distributions3")
m <- glm(goals ~ difference, data = FIFA2018, family = poisson)

## simulate new goals via glm method
set.seed(0)
g_glm <- simulate(m, n = 3)

## alternatively use the new default method
set.seed(0)
g_default <- simulate.default(m, n = 3)

## same results
all.equal(g_glm, g_default, check.attributes = FALSE)
#> [1] TRUE
```
