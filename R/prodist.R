#' Extracting fitted or predicted probability distributions from models
#' 
#' Generic function with methods for various model classes for extracting
#' fitted (in-sample) or predicted (out-of-sample) probability `distributions3`
#' objects.
#' 
#' To facilitate making probabilistic forecasts based on regression and time
#' series model objects, the function \code{prodist} extracts fitted or
#' predicted probability \code{distribution} objects. Currently, methods are
#' provided for objects fitted by \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, and \code{\link[stats]{arima}} in base R as
#' well as \code{glm.nb} from the \pkg{MASS} package and
#' \code{hurdle}/\code{zeroinfl}/\code{zerotrunc} from the \pkg{pscl} or
#' \pkg{countreg} packages.
#' 
#' All methods essentially
#' proceed in two steps: First, the standard \code{\link[stats]{predict}}
#' method for these model objects is used to compute fitted (in-sample, default)
#' or predicted (out-of-sample) distribution parameters. Typically, this includes
#' the mean plus further parameters describing scale, dispersion, shape, etc.).
#' Second, the \code{distributions} objects are set up using the generator
#' functions from \pkg{distributions3}.
#'
#' Note that these probability distributions only reflect the random variation in
#' the dependent variable based on the model employed (and its associated
#' distributional assumpation for the dependent variable). This does not capture
#' the uncertainty in the parameter estimates.
#' 
#' For both linear regression models and generalized linear models, estimated
#' by \code{lm} and \code{glm} respectively, there is some ambiguity as to which
#' estimate for the dispersion parameter of the model is to be used. While the
#' \code{\link[stats]{logLik}} methods use the maximum-likelihood (ML) estimate
#' implicitly, the \code{summary} methods report an estimate that is standardized
#' with the residual degrees of freedom, n - k (rather than the number of
#' observations, n). The \code{prodist} methods for these objects follow
#' the \code{logLik} method by default but the \code{summary} behavior can be
#' mimicked by setting the \code{sigma} or \code{dispersion} arguments
#' accordingly.
#' 
#' @aliases prodist.lm prodist.glm prodist.negbin prodist.Arima prodist.hurdle prodist.zeroinfl prodist.zerotrunc
#' 
#' @param object A model object.
#' @param ... Arguments passed on to methods, typically for calling the
#' underlying \code{\link[stats]{predict}} methods, e.g., \code{newdata} for
#' \code{\link[stats]{lm}} or \code{\link[stats]{glm}} objects or \code{n.ahead}
#' for \code{\link[stats]{arima}} objects.
#' @param sigma character or numeric or \code{NULL}. Specification of the standard
#' deviation \code{sigma} to be used for the \code{\link{Normal}} distribution in the
#' \code{lm} method. The default \code{"ML"} (or equivalently \code{"MLE"} or \code{NULL})
#' uses the maximum likelihood estimate based on the residual sum of squares divided
#' by the number of observations, n. Alternatively, \code{sigma = "OLS"} uses the
#' least-squares estimate (divided by the residual degrees of freedom, n - k). Finally,
#' a concrete numeric value can also be specified in \code{sigma}.
#' @param dispersion character or numeric or \code{NULL}. Specification of the
#' dispersion parameter in the \code{glm} method. The default \code{NULL}
#' (or equivalently \code{"deviance"}) is to use the \code{\link[stats]{deviance}}
#' divided by the number of observations, n. Alternatively, \code{dispersion = "Chisquared"}
#' uses the Chi-squared statistic divided by the residual degrees of freedom, n - k.
#' Finally, a concrete numeric value can also be specified in \code{dispersion}.
#' 
#' @return An object inheriting from \code{distribution}.
#' 
#' @seealso \code{\link[stats]{predict}}, \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, \code{\link[stats]{arima}}
#' 
#' @keywords distribution
#' 
#' @examples
#' ## Model: Linear regression
#' ## Fit: lm
#' ## Data: 1920s cars data
#' data("cars", package = "datasets")
#' 
#' ## Stopping distance (ft) explained by speed (mph)
#' reg <- lm(dist ~ speed, data = cars)
#' 
#' ## Extract fitted normal distributions (in-sample, with constant variance)
#' pd <- prodist(reg)
#' head(pd)
#' 
#' ## Extract log-likelihood from model object
#' logLik(reg)
#' 
#' ## Replicate log-likelihood via distributions object
#' sum(log_pdf(pd, cars$dist))
#' log_likelihood(pd, cars$dist)
#'
#' ## Compute corresponding medians and 90% interval
#' qd <- quantile(pd, c(0.05, 0.5, 0.95))
#' head(qd)
#' 
#' ## Visualize observations with predicted quantiles
#' plot(dist ~ speed, data = cars)
#' matplot(cars$speed, qd, add = TRUE, type = "l", col = 2, lty = 1)
#' 
#' ## Sigma estimated by maximum-likelihood estimate (default, used in logLik)
#' ## vs. least-squares estimate (used in summary)
#' nd <- data.frame(speed = 50)
#' prodist(reg, newdata = nd, sigma = "ML")
#' prodist(reg, newdata = nd, sigma = "OLS")
#' summary(reg)$sigma
#' 
#' 
#' ## Model: Poisson generalized linear model
#' ## Fit: glm
#' ## Data: FIFA 2018 World Cup data
#' data("FIFA2018", package = "distributions3")
#' 
#' ## Number of goals per team explained by ability differences
#' poisreg <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' summary(poisreg)
#' ## Interpretation: When the ratio of abilities increases by 1 percent,
#' ## the expected number of goals increases by around 0.4 percent
#' 
#' ## Predict fitted Poisson distributions for teams with equal ability (out-of-sample)
#' nd <- data.frame(difference = 0)
#' prodist(poisreg, newdata = nd)
#' 
#' ## Extract fitted Poisson distributions (in-sample)
#' pd <- prodist(poisreg)
#' head(pd)
#' 
#' ## Extract log-likelihood from model object
#' logLik(poisreg)
#' 
#' ## Replicate log-likelihood via distributions object
#' sum(log_pdf(pd, FIFA2018$goals))
#' log_likelihood(pd, FIFA2018$goals)
#' 
#' 
#' ## Model: Autoregressive integrated moving average model
#' ## Fit: arima
#' ## Data: Quarterly approval ratings of U.S. presidents (1945-1974)
#' data("presidents", package = "datasets")
#' 
#' ## ARMA(2,1) model
#' arma21 <- arima(presidents, order = c(2, 0, 1))
#' 
#' ## Extract predicted normal distributions for next two years
#' p <- prodist(arma21, n.ahead = 8)
#' p
#' 
#' ## Compute median (= mean) forecast along with 80% and 95% interval
#' quantile(p, c(0.5, 0.1, 0.9, 0.025, 0.975))
#' 
#' @export
prodist <- function(object, ...) {
  UseMethod("prodist")
}

#' @rdname prodist
#' @export
prodist.lm <- function(object, ..., sigma = "ML") {
  ## predicted means
  mu <- predict(object, ...)

  ## estimated sigma
  if(is.null(sigma)) sigma <- "ML"
  if(is.character(sigma)) {
    sigma <- match.arg(toupper(sigma), c("MLE", "OLS"))
    wts <- if(is.null(object$weights)) 1 else object$weights
    sigma <- switch(sigma,
      "MLE" = sqrt(sum(wts * residuals(object)^2)/nobs(object)),
      "OLS" = sqrt(sum(wts * residuals(object)^2)/df.residual(object)))
  }
  if(!is.numeric(sigma)) stop("unknown specification of sigma")

  ## set up distribution
  Normal(mu = mu, sigma = sigma)
}

#' @rdname prodist
#' @export
prodist.glm <- function(object, ..., dispersion = NULL) {
  ## fitted means
  mu <- predict(object, type = "response", ...)

  ## dispersion parameter phi
  if(is.null(dispersion)) {
    phi <- NULL
    dispersion <- "deviance"
  } else if(is.character(dispersion)) {
    phi <- NULL
    dispersion <- match.arg(gsub("-", "", tolower(dispersion), fixed = TRUE), c("deviance", "chisquared"))
  } else {
    phi <- dispersion
  }
  if(is.null(phi)) { 
    phi <- if(object$family$family %in% c("poisson", "binomial")) {
      1
    } else {
      switch(dispersion,
        "deviance" = deviance(object)/nobs(object),
        "chisquared" = sum((object$weights * object$residuals^2)[object$weights > 0])/df.residual(object))
    }
  }
  
  ## size for binomial distributions is available only in-sample:
  if(object$family$family == "binomial") {
    size <- if("newdata" %in% names(list(...))) 1L else object$prior.weights
  }
  
  ## special case: MASS::negative.binomial()
  if(startsWith(object$family$family, "Negative Binomial")) {
    object$family$family <- "negative.binomial"
    phi <- environment(object$family$variance)$.Theta
  }
  
  ## set up distributions object (if possible)
  switch(object$family$family,
    "gaussian" = Normal(mu = mu, sigma = sqrt(phi)),
    "poisson" = Poisson(lambda = mu),
    "binomial" = Binomial(size = size, p = mu),
    "Gamma" = distributions3::Gamma(shape = 1/phi, rate = 1/(phi * mu)),
    "negative.binomial" = NegativeBinomial(mu = mu, size = phi),
    "inverse.gaussian" = stop("inverse Gaussian distributions3 object not implemented yet"), ## FIXME: could use SuppDists for this
    "quasi" = stop("quasi family is not associated with a full probability distribution"),
    "quasibinomial" = stop("quasi-Poisson family is not associated with a full probability distribution"),
    "quasipoisson" = stop("quasi-Binomial family is not associated with a full probability distribution"),
    stop(sprintf("%s family not supported yet", object$family$family))
  )
}

#' @export
prodist.negbin <- function(object, ...) {
  mu <- predict(object, type = "response", ...)
  NegativeBinomial(mu = mu, size = object$theta)
}

#' @export
prodist.Arima <- function(object, ...) {
  p <- predict(object, ...)
  n <- stats::.preformat.ts(p$pred)
  n <- if(is.null(dim(n))) {
    format(stats::time(p$pred))
  } else {
    t(outer(rownames(n), colnames(n), paste))[t(n) != ""]
  }
  Normal(mu = setNames(as.numeric(p$pred), n), sigma = as.numeric(p$se))
}


## is countreg loaded (after pscl)?
.use_countreg <- function() {
  flavor <- c("pscl", "countreg")
  flavor <- setNames(match(flavor, .packages()), flavor)
  flavor <- names(sort(flavor))[1L]
  identical(flavor, "countreg")
}

#' @export
prodist.hurdle <- function(object, ...) {
  if(.use_countreg()) {
    dist <- predict(object, type = "distribution", ...)
  } else {
    dist <- object$dist$count
    mu <- predict(object, type = "count", ...)
    pi <- 1 - predict(object, type = "prob", at = 0:1, ...)[, 1L]
    dist <- switch(dist,
      "poisson"   = HurdlePoisson(lambda = mu, pi = pi),
      "geometric" = HurdleNegativeBinomial(mu = mu, theta = 1, pi = pi),
      "negbin"    = HurdleNegativeBinomial(mu = mu, theta = as.numeric(object$theta["count"]), pi = pi)
    )
  }
  return(dist)
}

#' @export
prodist.zeroinfl <- function(object, ...) {
  if(.use_countreg()) {
    dist <- predict(object, type = "distribution", ...)
  } else {
    dist <- object$dist
    mu <- predict(object, type = "count", ...)
    pi <- predict(object, type = "zero", ...)
    dist <- switch(dist,
      "poisson"   = ZIPoisson(lambda = mu, pi = pi),
      "geometric" = ZINegativeBinomial(mu = mu, theta = 1, pi = pi),
      "negbin"    = ZINegativeBinomial(mu = mu, theta = as.numeric(object$theta), pi = pi)
    )
  }
  return(dist)
}

#' @export
prodist.zerotrunc <- function(object, ...) {
  dist <- object$dist
  mu <- predict(object, type = "count", ...)
  switch(dist,
    "poisson"   = ZTPoisson(lambda = mu),
    "geometric" = ZTNegativeBinomial(mu = mu, theta = 1),
    "negbin"    = ZTNegativeBinomial(mu = mu, theta = as.numeric(object$theta))
  )
}

#' @rdname prodist
#' @export
prodist.distribution <- function(object, ...) object

## Further examples requiring other packages ---------------
## 
## library("MASS")
## m2 <- glm.nb(goals ~ difference, data = FIFA2018)
## pd2 <- prodist(m2)
## 
## data("cns", package = "faraway")
## cns$CNS <- cns$An + cns$Sp + cns$Other
## cns_logit <- glm(cbind(CNS, NoCNS) ~ Water + Work, data = cns, family = binomial)
## logLik(cns_logit)
## log_likelihood(prodist(cns_logit), cns$CNS)
## 
## library("pscl")
## data("bioChemists", package = "pscl")
## art_hnb <- hurdle(art ~ . | ., data = bioChemists, dist = "negbin")
## logLik(art_hnb)
## log_likelihood(prodist(art_hnb), bioChemists$art)
## all.equal(
##   pdf(prodist(art_hnb), 0:10),
##   predict(art_hnb, type = "prob", at = 0:10),
##   check.attributes = FALSE)
## 
## data("SwissLabor", package = "AER")
## swiss_logit <- glm(participation ~ . + I(age^2), data = SwissLabor, family = binomial)
## 
## 
## Multinomial and ordinal models --------------------------
## FIXME: Multinomial() distribution not yet vectorized
## 
## prodist.polr <- function(object, ...) {
##   p <- predict(object, type = "prob", ...)
##   Multinomial(size = 1, p = p)
## }
## 
## prodist.multinom <- function(object, ...) {
##   p <- predict(object, type = "prob", ...)
##   Multinomial(size = 1, p = p) ## FIXME predict size?
## }
##
## library("nnet")
## data("BankWages", package = "AER")
## bwmale <- subset(BankWages, gender == "male")
## bw_mnl <- multinom(job ~ education + minority, data = bwmale, trace = FALSE)
## 
## library("MASS")
## bw_olm <- polr(job ~ education + minority, data = bwmale, Hess = TRUE)
## summary(bw_olm)
