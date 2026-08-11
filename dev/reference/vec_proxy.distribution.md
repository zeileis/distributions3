# Methods for including distributions as vctrs in tibbles

Methods for
[`vec_proxy`](https://vctrs.r-lib.org/reference/vec_proxy.html) and
[`vec_restore`](https://vctrs.r-lib.org/reference/vec_proxy.html) from
vctrs in order to include `distribution` objects in
[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) objects.

## Usage

``` r
vec_proxy.distribution(x, ...)

vec_restore.distribution(x, to, ...)
```

## Arguments

- x, to:

  Objects inheriting from `distribution`.

- ...:

  Currently not used.

## Value

The `vec_proxy` method returns a `distribution` object which
additionally inherits of `data.frame` while the `vec_restore` method
restores the original `distribution` classes.

## Details

The methods for
[`vec_proxy`](https://vctrs.r-lib.org/reference/vec_proxy.html) and
[`vec_restore`](https://vctrs.r-lib.org/reference/vec_proxy.html) from
vctrs are needed so that `distribution` objects can be included as a
vector column in (and extracted from)
[`tibble`](https://tibble.tidyverse.org/reference/tibble.html) data
frames. `vec_proxy` simply adds the class `data.frame` which is the
actual underlying data structure used by `distribution` objects. This
way the number of rows etc. can be correctly determined. Conversely,
`vec_restore` strips off the additional `data.frame` class and restores
the original `distribution` classes. Users typically do not need to call
`vec_proxy` and `vec_restore` directly.

## Examples

``` r
## Poisson GLM for FIFA 2018 goals data
data("FIFA2018", package = "distributions3")
m <- glm(goals ~ difference, data = FIFA2018, family = poisson)

## Predict fitted Poisson distributions for teams with ability differences
## of -1, 0, 1 (out-of-sample) using the new data as a data.frame
nd <- data.frame(difference = -1:1)
nd$dist <- prodist(m, newdata = nd)
nd
#>   difference                        dist
#> 1         -1 Poisson(lambda = 0.8181454)
#> 2          0 Poisson(lambda = 1.2370397)
#> 3          1 Poisson(lambda = 1.8704100)

## Do the same using the new data as a tibble
library("tibble")
nt <- tibble(difference = -1:1)
nt$dist <- prodist(m, newdata = nt)
nt
#> # A tibble: 3 × 2
#>   difference dist                    
#>        <int> <Poisson>               
#> 1         -1 Poisson(lambda = 0.8181)
#> 2          0 Poisson(lambda = 1.2370)
#> 3          1 Poisson(lambda = 1.8704)
```
