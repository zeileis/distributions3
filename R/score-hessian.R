## score/hessian generics
score <- function(d, ...) UseMethod("score")
hessian <- function(d, ...) UseMethod("hessian")

## fallback methods based on numeric differentiation
score.distribution <- function(d, x, which = NULL, drop = TRUE, eps = .Machine$double.eps^(1/3), ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- names(unclass(d))
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## compute scores
  scr <- function(par) {
    d1 <- d2 <- d
    d1[[par]] <- d1[[par]] + eps
    d2[[par]] <- d2[[par]] - eps
    (log_pdf(d1, x) - log_pdf(d2, x)) / (2 * eps)
  }

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

hessian.distribution <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, eps = .Machine$double.eps^(1/4), ...) {
  ## numeric differentiation yields observed hessian only
  if (!identical(expected, FALSE)) stop("only the observed hessian is available")
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- names(unclass(d))
  pp <- outer(p, p, paste, sep = ":")
  diag(pp) <- p
  p <- setNames(
    c(diag(pp), pp[upper.tri(pp)], pp[upper.tri(pp)]),
    c(diag(pp), pp[upper.tri(pp)], pp[lower.tri(pp)])
  )[pp]
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## compute scores
  hess <- function(par) {
    par <- strsplit(par, ":", fixed = TRUE)[[1L]]
    par <- rep_len(par, 2L)
    d1 <- d2 <- d
    d1[[par[2L]]] <- d1[[par[2L]]] + eps
    d2[[par[2L]]] <- d2[[par[2L]]] - eps
    (score(d1, x, which = par[1L]) - score(d2, x, which = par[1L])) / (2 * eps)
  }

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}


## Normal methods for score/hessian
score.Normal <- function(d, x, which = NULL, drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## available and selected parameters
  p <- c("mu", "sigma")
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)

  ## compute scores
  scr <- function(par) switch(par,
    "mu"    = (x - d$mu)/(d$sigma^2),
    "sigma" = (x - d$mu)^2/(d$sigma^3) - 1/d$sigma)

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- setNames(scr(which), names(d))
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    dimnames(s) <- list(names(d), which)
  }
  return(s)
}

hessian.Normal <- function(d, x, which = NULL, drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu" = "mu", "sigma:mu" = "mu:sigma", "mu:sigma" = "mu:sigma", "sigma" = "sigma")
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## function for computing Hessian elements (expected or observed)
  hess <- if (expected) {
    function(par) switch(par,
      "mu"    = rep_len(-1/d$sigma^2, n),
      "sigma" = rep_len(-2 /d$sigma^2, n),
      rep.int(0, n))
  } else {
    function(par) switch(par,
      "mu"    = rep_len(-1/d$sigma^2, n),
      "sigma" = -3 * (x - d$mu)^2/(d$sigma^4) + 1/d$sigma^2,
      -2 * (x - d$mu)/d$sigma^3)
  }
  
  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- setNames(hess(w), names(d))
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    dimnames(h) <- list(names(d), w)
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}

## Poisson methods for score/hessian
score.Poisson <- function(d, x, which = "lambda", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter
  which <- match.arg(which, "lambda", several.ok = TRUE)

  ## compute score
  s <- x/d$lambda - 1
  if (!drop) s <- cbind("lambda" = s)
  return(s)
}

hessian.Poisson <- function(d, x, which = "lambda", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter
  which <- match.arg(which, "lambda", several.ok = TRUE)

  ## compute hessian
  h <- if (expected) rep_len(-1/d$lambda, n) else -x/d$lambda^2
  if (!drop) h <- cbind("lambda" = h)
  return(h)
}

## Bernoulli methods for score/hessian
score.Bernoulli <- function(d, x, which = "p", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter
  which <- match.arg(which, "p", several.ok = TRUE)

  ## compute score
  s <- (x - d$p)/(d$p * (1 - d$p))
  if (!drop) s <- cbind("p" = s)
  return(s)
}

hessian.Bernoulli <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter
  which <- match.arg(which, "p", several.ok = TRUE)

  ## compute hessian
  h <- if (expected) rep_len(-1/(d$p * (1 - d$p)), n) else -x/d$p^2 - (1 - x)/(1 - d$p)^2
  if (!drop) h <- cbind("p" = h)
  return(h)
}

## Binomial methods for score/hessian
score.Binomial <- function(d, x, which = "p", drop = TRUE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")

  ## only one parameter supported
  which <- match.arg(which, c("p", "size"), several.ok = TRUE)
  if (!identical(which, "p")) warning("only the scores with respect to 'p' are supported")

  ## compute score
  s <- (x - d$size * d$p)/(d$p * (1 - d$p))
  if (!drop) s <- cbind("p" = s)
  return(s)
}

hessian.Binomial <- function(d, x, which = "p", drop = TRUE, expected = FALSE, ...) {
  ## sanity check
  n <- c(length(d), length(x))
  if (n[1L] != n[2L] && all(n > 1L)) stop("'d' and 'x' must have length 1 or the same length")
  n <- max(n)

  ## only one parameter supported
  which <- match.arg(which, c("p", "size"), several.ok = TRUE)
  if (!identical(which, "p")) warning("only the scores with respect to 'p' are supported")

  ## compute hessian
  h <- if (expected) rep_len(-d$size/(d$p * (1 - d$p)), n) else -x/d$p^2 - (d$size - x)/(1 - d$p)^2
  if (!drop) h <- cbind("p" = h)
  return(h)
}
