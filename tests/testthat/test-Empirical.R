# -------------------------------------------------------
# Checking Empirical distribution
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

test_that("Empirical exists", {
    expect_true(is.function(Empirical))
})

test_that("Empirical unexpected input", {
    # All converted to NA
    expect_error(expect_warning(Empirical(LETTERS[1:3])),
        regexp = "at least two finite observations must be available")
    # Only one is finite
    expect_error(Empirical(c(-Inf, 3, Inf)),
        regexp = "at least two finite observations must be available")
})

test_that("Empirical construction function from vector", {
    x <- sample(1:10)
    expect_silent(d <- Empirical(x))
    expect_identical(length(d), 1L)
    expect_null(names(d))
    expect_output(print(d), regexp = "Empirical distribution \\(.*\\)",)

    # Converting to matrix
    expect_silent(m <- as.matrix(d))
    expect_true(is.matrix(m))
    expect_identical(dim(m), c(1L, 10L))
    expect_identical(dimnames(m), list(NULL, paste0("o_", 1:10)))
})

test_that("Empirical construction from unnamed numeric matrix", {
    x <- matrix(sample(1:10), nrow = 2)
    expect_silent(d <- Empirical(x))
    expect_identical(length(d), 2L)
    expect_null(names(d))
    expect_output(print(d), regexp = "Empirical distribution \\(.*\\)",)

    # Converting to matrix
    expect_silent(m <- as.matrix(d))
    expect_true(is.matrix(m))
    expect_identical(dim(m), c(2L, 5L))
    expect_identical(dimnames(m), list(NULL, paste0("o_", 1:5)))
})

test_that("Empirical construction from unnamed character matrix", {
    x <- matrix(as.character(sample(1:10)), nrow = 2)
    expect_silent(d <- Empirical(x))
    expect_identical(length(d), 2L)
    expect_null(names(d))
})

test_that("Empirical construction from named matrix", {
    x <- matrix(as.character(sample(1:10)), nrow = 2,
                dimnames = list(letters[1:2], LETTERS[1:5]))
    expect_silent(d <- Empirical(x))
    expect_identical(length(d), 2L)
    expect_identical(names(d), letters[1:2])

    # Converting to matrix
    expect_silent(m <- as.matrix(d))
    expect_true(is.matrix(m))
    expect_identical(dim(m), c(2L, 5L))
    expect_identical(dimnames(m), list(letters[1:2], LETTERS[1:5]))
    expect_true(!any(is.na(m)))
})

test_that("Empirical construction from data.frame (list)", {
    x <- data.frame(a = 1:3, b = 11:13)
    expect_silent(d <- Empirical(x))
    expect_identical(length(d), 2L)
    expect_identical(names(d), letters[1:2])

    # Converting to matrix
    expect_silent(m <- as.matrix(d))
    expect_identical(m, matrix(as.double(c(1:3, 11:13)), byrow = TRUE,
                               nrow = 2, dimnames = list(letters[1:2], paste0("o_", 1:3))))
})

test_that("Empirical construction from list, unequal number of observations", {
    # Two distributions w/ different number of observations
    x <- list(a = sample(1:10), b = sample(1:10, 4))
    expect_silent(d <- Empirical(x))
    expect_identical(length(d), 2L)
    expect_identical(names(d), letters[1:2])

    # Converting to matrix
    expect_silent(m <- as.matrix(d))
    expect_true(is.matrix(m))
    expect_identical(dim(m), c(2L, 10L))
    expect_identical(dimnames(m), list(letters[1:2], paste0("o_", 1:10)))
    expect_identical(sum(is.na(m)), 6L)
    expect_true(all(is.na(m[2, 5:10])))
})


test_that("Empirical support method works", {
    # Checking defaults
    expect_identical(formals(distributions3:::support.Empirical),
        as.pairlist(alist(d =, drop = TRUE, ... =)))

    d <- Empirical(list(a = 8:20, b = c(-1.2, -1, -0.8)))
    expect_identical(support(d),
        matrix(c(8, 20, -1.2, -0.8), byrow = TRUE, nrow = 2, dimnames = list(letters[1:2], c("min", "max"))))

    expect_identical(support(d[1], drop = FALSE),
        matrix(c(8, 20), byrow = TRUE, nrow = 1, dimnames = list(letters[1], c("min", "max"))))

    expect_identical(support(d[1]), c(8, 20) |> setNames(c("min", "max")))
})

test_that("is_continuous/is_discrete methods both works", {
    # Checking defaults
    expect_identical(formals(distributions3:::is_continuous.Empirical),
        as.pairlist(alist(d =, ... =)))

    d <- Empirical(list(a = 8:20, b = c(-1.2, -1, -0.8)))
    expect_silent(x <- is_continuous(d))
    expect_identical(x, rep(FALSE, 2L) |> setNames(letters[1:2]))
    expect_silent(x <- is_discrete(d))
    expect_identical(x, rep(FALSE, 2L) |> setNames(letters[1:2]))
})


# -------------------------------------------------------------------
# S3 methods for distribution functions and moments
# -------------------------------------------------------------------

set.seed(1234)
a <- c(-2, -1, -1, 1, 1, 2)
b <- rnorm(15, mean = 10, sd = 10)
d <- Empirical(list(a = a, b = b))

test_that("cdf.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::cdf.Empirical),
        as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, ... =)))

    expect_silent(x <- cdf(d, 0))
    expect_identical(x, c(0.5, 2 / 15) |> setNames(letters[1:2]))

    expect_silent(x <- cdf(d, c(-2, 0, 2)))
    expect_identical(x, matrix(c(1/6, 1/2, 6/6, 2/15, 2/15, 4/15),
            nrow = 2, byrow = TRUE, dimnames = list(letters[1:2], paste0("p_", c(-2, 0, 2)))))

    # Upper tail
    expect_silent(x2 <- cdf(d, c(-2, 0, 2), lower.tail = FALSE))
    expect_identical(x2, 1 - x)
})


test_that("pdf.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::pdf.Empirical),
        as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, ... =)))

    expect_silent(x <- pdf(d, 1))
    expect_identical(x, c(2 / 6, 0) |> setNames(letters[1:2]))

    expect_silent(x <- pdf(d, c(-2, 0, 1)))
    expect_equal(x, matrix(c(1/6, 0, 2/6, 0, 0, 0),
            nrow = 2, byrow = TRUE, dimnames = list(letters[1:2], paste0("d_", c(-2, 0, 1)))))

    # Log-pdf
    expect_silent(x2 <- pdf(d, c(-2, 0, 1), log = TRUE))
    expect_identical(x2, log(x))

    # log_pdf function
    expect_identical(formals(distributions3:::log_pdf.Empirical),
        as.pairlist(alist(d =, x =, drop = TRUE, elementwise = NULL, ... =)))
    expect_silent(x3 <- log_pdf(d, c(-2, 0, 1)))
    expect_equal(x2, x3, ignore_attr = TRUE)
    expect_identical(dimnames(x3), list(letters[1:2], paste0("l_", c(-2, 0, 1))))
})


test_that("quantile.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::quantile.Empirical),
        as.pairlist(alist(x =, probs = , drop = TRUE, elementwise = NULL, ... =)))

    expect_silent(x <- quantile(d, 0))
    expect_identical(x, c(min(a), min(b)) |> setNames(letters[1:2]))

    expect_silent(x <- quantile(d, 1))
    expect_identical(x, c(max(a), max(b)) |> setNames(letters[1:2]))

    expect_silent(x <- quantile(d, 0.5))
    expect_identical(x, c(unname(quantile(a, 0.5)), unname(quantile(b, 0.5))) |>
                     setNames(letters[1:2]))

    expect_silent(x <- quantile(d, c(0.2, 0.5)))
    expect_identical(x, c(unname(quantile(a, 0.2)), unname(quantile(b, 0.5))) |>
                     setNames(letters[1:2]))

    expect_silent(x <- quantile(d, c(0.2, 0.5)))
    expect_identical(x, c(unname(quantile(a, 0.2)), unname(quantile(b, 0.5))) |>
                     setNames(letters[1:2]))

    expect_silent(x <- quantile(d, c(0.2, 0.5), elementwise = FALSE))
    expect_identical(x, rbind(quantile(a, c(0.2, 0.5)), quantile(b, c(0.2, 0.5))) |>
        structure(dimnames = list(letters[1:2], paste0("q_", c(0.2, 0.5)))))

    # Lower tail
    expect_silent(x2 <- quantile(d, c(0.2, 0.5), elementwise = FALSE, lower.tail = FALSE))
    expect_identical(x2, 1 - x)

})


test_that("random.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::random.Empirical),
            as.pairlist(alist(x =, n = 1L, drop = TRUE, ... = )))

    # n = 0
    expect_silent(x <- random(d, n = 0L))
    expect_identical(x, vector("double", 0))

    # n = 1 (default)
    expect_silent(x <- random(d)) # Default n = 1
    expect_type(x, "double")
    expect_true(is.vector(x) && length(x) == 2L)
    expect_identical(names(x), letters[1:2])
    expect_true(x[[1L]] %in% a && x[[2L]] %in% b)

    # n = 10
    expect_silent(x <- random(d, n = 10L))
    expect_type(x, "double")
    expect_true(is.matrix(x))
    expect_identical(dim(x), c(2L, 10L))
    expect_identical(dimnames(x), list(letters[1:2], paste0("r_", 1:10)))
    expect_true(all(x[1L, ] %in% a & x[2L, ] %in% b))

})


test_that("mean.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::mean.Empirical),
            as.pairlist(alist(x =, ... =)))

    expect_silent(x <- mean(d))
    expect_type(x, "double")
    expect_identical(x, c(mean(a, na.rm = TRUE), mean(b)) |>
                     setNames(letters[1:2]))

    expect_silent(x <- variance(d))
    expect_type(x, "double")
    expect_equal(x, c(sd(a, na.rm = TRUE)^2, sd(b)^2) |>
                 setNames(letters[1:2]))
})


test_that("skewness.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::skewness.Empirical),
            as.pairlist(alist(x =, type = 1L, ... =)))

    sk1 <- function(x) {
        x <- na.omit(x)
        sqrt(length(x)) * sum((x - mean(x))^3) / sqrt(sum((x - mean(x))^2)^3)
    }
    sk2 <- function(x) {
        x <- na.omit(x)
        sk1(x) * sqrt(length(x) * (length(x) - 1L)) / (length(x) - 2L)
    }
    sk3 <- function(x) {
        x <- na.omit(x)
        sk1(x) * sqrt((1 - 1 / length(x))^3)
    }

    expect_equal(skewness(d, type = 1L), c(sk1(a), sk1(b)) |> setNames(letters[1:2]))
    expect_equal(skewness(d, type = 2L), c(sk2(a), sk2(b)) |> setNames(letters[1:2]))
    expect_equal(skewness(d, type = 3L), c(sk3(a), sk3(b)) |> setNames(letters[1:2]))

    expect_error(skewness(d, type = 0), regexp = "invalid 'type' argument")
    expect_error(skewness(d, type = 4), regexp = "invalid 'type' argument")
    expect_error(skewness(Empirical(1:2), type = 3L),
        regexp = "at least three finite observations per distribution required to calculate type 2 skewness")
    expect_error(skewness(Empirical(c(-Inf, -Inf, 1, 2, Inf)), type = 3L),
        regexp = "at least three finite observations per distribution required to calculate type 2 skewness")
})


test_that("kurtosis.Empirical works as expected", {
    # Checking defaults
    expect_identical(formals(distributions3:::kurtosis.Empirical),
            as.pairlist(alist(x =, type = 3L, ... =)))

    ku1 <- function(x) {
        x <- na.omit(x)
        x2 <- (x - mean(x))^2
        length(x) * sum(x2^2) / sum(x2)^2 - 3
    }
    ku2 <- function(x) {
        x <- na.omit(x)
        n <- length(x)
        ((n + 1) * ku1(x) + 6) * (n - 1) / ((n - 2) * (n - 3))
    }
    ku3 <- function(x) {
        x <- na.omit(x)
        n <- length(x)
        (1 - 1 / n)^2 * (ku1(x) + 3) - 3
    }

    expect_equal(kurtosis(d, type = 1L), c(ku1(a), ku1(b)) |> setNames(letters[1:2]))
    expect_equal(kurtosis(d, type = 2L), c(ku2(a), ku2(b)) |> setNames(letters[1:2]))
    expect_equal(kurtosis(d, type = 3L), c(ku3(a), ku3(b)) |> setNames(letters[1:2]))
    expect_error(kurtosis(d, type = 10), regexp = "invalid 'type' argument")

    expect_error(kurtosis(Empirical(1:3), type = 0L), regexp = "invalid 'type' argument")
    expect_error(kurtosis(Empirical(1:3), type = 4L), regexp = "invalid 'type' argument")
    expect_error(kurtosis(Empirical(1:3), type = 2L),
        regexp = "at least four finite observations per distribution required to calculate type 2 kurtosis")
    expect_error(kurtosis(Empirical(c(-Inf, -Inf, 1, 2, 3)), type = 2L),
        regexp = "at least four finite observations per distribution required to calculate type 2 kurtosis")
})


# -------------------------------------------------------------------
# Additional tests for dedicated dpqr methods
# -------------------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }
suppressPackageStartupMessages(library("scoringRules"))

# dempirical
test_that("dempirical works as expected", {
    expect_identical(formals(dempirical),
        as.pairlist(alist(x =, y =, log = FALSE, na.rm = FALSE, method = NULL, ... = )))

    y <- c(1, 1, NA_real_, 2, 2, 2, 3, 4, 5, NA_real_)
    n <- sum(!is.na(y))

    # Incorrect argument
    expect_error(dempirical(0, y, method = "foo"), regexp = "'arg' should be one of .hist., .density.")

    # Using default 'method = NULL' (discrete)
    expect_identical(dempirical(0, y), NA_real_)
    expect_identical(dempirical(1:2, y), rep(NA_real_, 2L))

    expect_identical(dempirical(0, y, na.rm = TRUE), 0)
    expect_identical(dempirical(1, y, na.rm = TRUE), 2 / n)
    expect_identical(x <- dempirical(0:6, y, na.rm = TRUE),
            sapply(0:6, function(x) { sum(y == x, na.rm = TRUE) / sum(!is.na(y)) }))

    expect_identical(log(x), dempirical(0:6, y, na.rm = TRUE, log = TRUE))

    # Using histogram estimation
    expect_equal(x <- dempirical(1.5, y, method = "hist"), 0.625)
    expect_equal(log(x), dempirical(1.5, y, method = "hist", log = TRUE))

    # Using density estimation
    expect_equal(x <- dempirical(1.5, y, method = "density"), 0.2885984, tolerance = 1e-6)
    expect_equal(log(x), dempirical(1.5, y, method = "density", log = TRUE))
})

# pempirical
test_that("pempirical works as expected", {
    expect_identical(formals(pempirical),
        as.pairlist(alist(q =, y = , lower.tail = TRUE, log.p = FALSE, na.rm = TRUE)))

    set.seed(1234)
    y <- c(rnorm(50), NA_real_)
    expect_identical(x <- pempirical(0, y), 0.76)
    expect_identical(pempirical(0, y, log.p = TRUE), log(x))
    expect_identical(pempirical(0, y, lower.tail = FALSE), 1 - x)
    expect_identical(pempirical(0, y, lower.tail = FALSE, log.p = TRUE), log(1 - x))

    # length of q > 1L
    expect_identical(pempirical((-1):1, y), c(0.22, 0.76, 0.92))

    # Not removing NAs - expecting NA as return
    expect_true(is.na(pempirical(0, y, na.rm = FALSE)))
})

# qempirical
test_that("qempirical works as expected", {
    expect_identical(formals(qempirical),
        as.pairlist(alist(p =, y = , lower.tail = TRUE, log.p = FALSE, na.rm = TRUE, ... =)))

    set.seed(1234)
    y <- c(rnorm(50), NA_real_)

    expect_error(qempirical(0.77, y, na.rm = FALSE), regexp = "missing values and NaN's not allowed if 'na.rm' is FALSE")

    expect_identical(x <- qempirical(0.77, y), unname(quantile(y, 0.77, na.rm = TRUE)))
    expect_identical(qempirical(log(0.77), y, log.p = TRUE), x)
    expect_identical(qempirical(0.77, y, lower.tail = FALSE), 1 - x)
    expect_identical(qempirical(log(0.77), y, lower.tail = FALSE, log.p = TRUE), 1 - x)

    # length of p > 1L
    expect_identical(x <- qempirical(c(0.2, 0.5, 0.8), y), unname(quantile(y, c(0.2, 0.5, 0.8), na.rm = TRUE)))

    # Not removing NAs - expecting error)
    expect_error(qempirical(0.5, y, na.rm = FALSE),
        regexp = "issing values and NaN's not allowed if 'na.rm' is FALSE")
})

# rempirical
test_that("rempirical works as expected", {
    expect_identical(formals(rempirical),
        as.pairlist(alist(n =, y = , na.rm = TRUE)))

    set.seed(1234)
    y <- c(rnorm(50), NA_real_)

    # n = 0L:
    expect_identical(rempirical(0, y), vector("double", 0))

    # n = 1L
    expect_silent(x <- rempirical(1L, na.omit(y)))
    expect_type(x, "double")
    expect_true(is.vector(x) && length(x) == 1L)
    expect_true(x %in% y)

    # n = 10L
    expect_silent(x <- rempirical(10L, na.omit(y)))
    expect_type(x, "double")
    expect_true(is.vector(x) && length(x) == 10L)
    expect_true(all(x %in% y))

    # Allowing to draw NAs
    expect_silent(x <- rempirical(n = 100L, c(1, 2,  rep(NA_real_, 10L)), na.rm = FALSE))
    expect_true(any(is.na(x)))

    # Removing NAs
    expect_silent(x <- rempirical(n = 100L, c(1, 2,  rep(NA_real_, 10L)), na.rm = TRUE))
    expect_true(all(x %in% 1:2))
})





