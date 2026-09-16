# -------------------------------------------------------
# Testing utils: auxiliary/helper functions and methods
# -------------------------------------------------------

if (interactive()) { library("distributions3"); library("testthat") }

test_that("is_distribution() works", {
  expect_true(is_distribution(Normal()))
  expect_false(is_distribution(123))
})

test_that("{methods}.dstribution work", {
  n <- Normal(c(0, 10), c(1, 1))

  expect_null(dim.distribution(n))
  expect_equal(length.distribution(n), 2)
  expect_equal(n[1], Normal(0, 1))
  expect_length(format.distribution(n), 2L)
  expect_null(names.distribution(n))

  expect_silent(n <- `names<-.distribution`(n, c("a", "b")))
  expect_named(n, c("a", "b"))
  expect_silent(names(n) <- NULL)
  expect_equal(as.matrix.distribution(n), as.matrix(data.frame(mu = c(0, 10), sigma = c(1, 1))))
  df <- data.frame(n = 1:2)
  df$n <- n
  expect_equal(as.data.frame.distribution(n), df)
  expect_equal(as.list.distribution(n), list(mu = c(0, 10), sigma = c(1, 1)))
  expect_equal(n, c.distribution(n[1], n[2]))
  expect_equal(
    capture_output(print(summary(n))),
    capture_output({
      cat("Normal distribution: \n")
      print(summary.data.frame(as.data.frame(as.matrix(n))))
    })
  )
})

test_that("apply_dpqr() applied to 'random' works", {
  N <- Normal(c(0, 10, 100), 1)
  N_named <- N
  names(N_named) <- LETTERS[1:length(N)]

  ## length(d) = 1, n = 1, drop = TRUE
  expect_true(is.numeric(random(N[1], 1)))
  expect_null(dim(random(N[1], 1)))
  expect_length(random(N[1], 1), 1)
  expect_equal(
    {
      set.seed(123)
      random(N[1], 1)
    },
    {
      set.seed(123)
      drop(random(N[1], 1))
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N_named[1], 1)
    },
    {
      set.seed(123)
      drop(random(N_named[1], 1))
    }
  )
  expect_null(names(random(N[1], 1)))
  expect_equal(names(random(N_named[1], 1)), "A")

  ## length(d) = 1, n = 1, drop = FALSE
  expect_true(is.numeric(random(N[1], 1, drop = FALSE)))
  expect_equal(dim(random(N[1], 1, drop = FALSE)), c(1L, 1L))
  expect_equal(colnames(random(N[1], 1, drop = FALSE)), "r_1")
  expect_null(rownames(random(N[1], 1, drop = FALSE)))
  expect_equal(colnames(random(N_named[1], 1, drop = FALSE)), "r_1")
  expect_equal(rownames(random(N_named[1], 1, drop = FALSE)), "A")

  ## length(d) > 1, n = 1, drop = TRUE
  expect_true(is.numeric(random(N, 1)))
  expect_null(dim(random(N, 1)))
  expect_length(random(N, 1), length(N))
  expect_equal(
    {
      set.seed(123)
      random(N, 1)
    },
    {
      set.seed(123)
      drop(random(N, 1))
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N_named, 1)
    },
    {
      set.seed(123)
      drop(random(N_named, 1))
    }
  )
  expect_null(names(random(N, 1)))
  expect_equal(names(random(N_named, 1)), c("A", "B", "C"))

  ## length(d) > 1, n = 1, drop = FALSE
  expect_true(is.numeric(random(N, 1, drop = FALSE)))
  expect_equal(dim(random(N, 1, drop = FALSE)), c(3L, 1L))
  expect_equal(colnames(random(N, 1, drop = FALSE)), "r_1")
  expect_null(rownames(random(N, 1, drop = FALSE)))
  expect_equal(colnames(random(N_named, 1, drop = FALSE)), "r_1")
  expect_equal(rownames(random(N_named, 1, drop = FALSE)), c("A", "B", "C"))

  ## length(d) = 1, n > 1, drop = TRUE
  expect_true(is.numeric(random(N[1], 2)))
  expect_null(dim(random(N[1], 2)))
  expect_length(random(N[1], 2), 2)
  expect_equal(
    {
      set.seed(123)
      random(N[1], 2)
    },
    {
      set.seed(123)
      drop(random(N[1], 2))
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N_named[1], 2)
    },
    {
      set.seed(123)
      drop(random(N_named[1], 2))
    }
  )
  expect_null(names(random(N[1], 2)))
  expect_null(names(random(N_named[1], 2)))

  ## length(d) = 1, n > 1, drop = FALSE
  expect_true(is.numeric(random(N[1], 2, drop = FALSE)))
  expect_equal(dim(random(N[1], 2, drop = FALSE)), c(1L, 2L))
  expect_equal(colnames(random(N[1], 2, drop = FALSE)), c("r_1", "r_2"))
  expect_null(rownames(random(N[1], 2, drop = FALSE)))
  expect_equal(colnames(random(N_named[1], 2, drop = FALSE)), c("r_1", "r_2"))
  expect_equal(rownames(random(N_named[1], 2, drop = FALSE)), "A")

  ## length(d) = n > 1, drop = TRUE
  expect_equal(
    {
      set.seed(123)
      random(N[1:2], 2)
    },
    {
      set.seed(123)
      drop(random(N[1:2], 3)[, -3L, drop = FALSE])
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N_named[1:2], 2)
    },
    {
      set.seed(123)
      drop(random(N_named[1:2], 3)[, -3L, drop = FALSE])
    }
  )

  ## length(d) = n > 1, drop = FALSE
  expect_equal(
    {
      set.seed(123)
      random(N[1:2], 2, drop = FALSE)
    },
    {
      set.seed(123)
      drop(random(N[1:2], 3, drop = FALSE)[, -3L, drop = FALSE])
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N_named[1:2], 2, drop = FALSE)
    },
    {
      set.seed(123)
      drop(random(N_named[1:2], 3, drop = FALSE)[, -3L, drop = FALSE])
    }
  )

  ## length(d) > 1, n > 1, drop = TRUE
  expect_true(is.numeric(random(N, 2)))
  expect_equal(dim(random(N, 2)), c(3L, 2L))
  expect_equal(
    {
      set.seed(123)
      random(N, 2)
    },
    {
      set.seed(123)
      drop(random(N, 2))
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N_named, 2)
    },
    {
      set.seed(123)
      drop(random(N_named, 2))
    }
  )
  expect_equal(colnames(random(N, 2)), c("r_1", "r_2"))
  expect_null(rownames(random(N, 2)))
  expect_equal(colnames(random(N_named, 2)), c("r_1", "r_2"))
  expect_equal(rownames(random(N_named, 2)), c("A", "B", "C"))

  ## length(d) > 1, n > 1, drop = FALSE
  expect_true(is.numeric(random(N, 2, drop = FALSE)))
  expect_equal(dim(random(N, 2, drop = FALSE)), c(3L, 2L))
  expect_equal(
    {
      set.seed(123)
      random(N, cbind(2), drop = FALSE)
    },
    {
      set.seed(123)
      random(N, 2, drop = FALSE)
    }
  )
  expect_equal(
    {
      set.seed(123)
      random(N, rbind(2), drop = FALSE)
    },
    {
      set.seed(123)
      random(N, 2, drop = FALSE)
    }
  )
  expect_equal(colnames(random(N, 2, drop = FALSE)), c("r_1", "r_2"))
  expect_null(rownames(random(N, 2, drop = FALSE)))
  expect_equal(colnames(random(N_named, 2, drop = FALSE)), c("r_1", "r_2"))
  expect_equal(rownames(random(N_named, 2, drop = FALSE)), c("A", "B", "C"))
})

test_that("apply_dpqr() applied to 'pdf', 'log_pdf' and 'cdf' works", {
  N <- Normal(c(0, 10, 100), 1)
  N_named <- N
  names(N_named) <- LETTERS[1:length(N)]

  ## length(d) = 1, at = 1, drop = TRUE
  expect_true(is.numeric(pdf(N[1], 0.5)))
  expect_null(dim(pdf(N[1], 0.5)))
  expect_length(pdf(N[1], 0.5), 1)
  expect_equal(pdf(N[1], cbind(0.5)), pdf(N[1], 0.5))
  expect_equal(pdf(N[1], rbind(0.5)), pdf(N[1], 0.5))
  expect_equal(pdf(N[1], 0.5), drop(pdf(N[1], 0.5)))
  expect_equal(pdf(N_named[1], 0.5), drop(pdf(N_named[1], 0.5)))
  expect_null(names(pdf(N[1], 0.5)))
  expect_equal(names(pdf(N_named[1], 0.5)), "A")

  ## length(d) = 1, at = 1, drop = FALSE
  expect_true(is.numeric(pdf(N[1], 0.5, drop = FALSE)))
  expect_equal(dim(pdf(N[1], 0.5, drop = FALSE)), c(1L, 1L))
  expect_equal(pdf(N[1], cbind(0.5), drop = FALSE), pdf(N[1], 0.5, drop = FALSE))
  expect_equal(pdf(N[1], rbind(0.5), drop = FALSE), pdf(N[1], 0.5, drop = FALSE))
  expect_equal(colnames(pdf(N[1], 0.5, drop = FALSE)), "d_0.5")
  expect_null(rownames(pdf(N[1], 0.5, drop = FALSE)))
  expect_equal(colnames(pdf(N_named[1], 0.5, drop = FALSE)), "d_0.5")
  expect_equal(rownames(pdf(N_named[1], 0.5, drop = FALSE)), "A")

  ## length(d) > 1, at = 1, drop = TRUE
  expect_true(is.numeric(pdf(N, 0.5)))
  expect_null(dim(pdf(N, 0.5)))
  expect_length(pdf(N, 0.5), length(N))
  expect_equal(pdf(N, cbind(0.5)), pdf(N, 0.5))
  expect_equal(pdf(N, rbind(0.5)), pdf(N, 0.5))
  expect_equal(pdf(N, 0.5), drop(pdf(N, 0.5)))
  expect_equal(pdf(N_named, 0.5), drop(pdf(N_named, 0.5)))
  expect_null(names(pdf(N, 0.5)))
  expect_equal(names(pdf(N_named, 0.5)), c("A", "B", "C"))

  ## length(d) > 1, at = 1, drop = FALSE
  expect_true(is.numeric(pdf(N, 0.5, drop = FALSE)))
  expect_equal(dim(pdf(N, 0.5, drop = FALSE)), c(3L, 1L))
  expect_equal(pdf(N, cbind(0.5), drop = FALSE), pdf(N, 0.5, drop = FALSE))
  expect_equal(pdf(N, rbind(0.5), drop = FALSE), pdf(N, 0.5, drop = FALSE))
  expect_equal(colnames(pdf(N, 0.5, drop = FALSE)), "d_0.5")
  expect_null(rownames(pdf(N, 0.5, drop = FALSE)))
  expect_equal(colnames(pdf(N_named, 0.5, drop = FALSE)), "d_0.5")
  expect_equal(rownames(pdf(N_named, 0.5, drop = FALSE)), c("A", "B", "C"))

  ## length(d) = 1, at > 1, drop = TRUE
  expect_true(is.numeric(pdf(N[1], c(0.2, 0.5))))
  expect_null(dim(pdf(N[1], c(0.2, 0.5))))
  expect_length(pdf(N[1], c(0.2, 0.5)), 2)
  expect_equal(pdf(N[1], cbind(c(0.2, 0.5))), pdf(N[1], c(0.2, 0.5)))
  expect_equal(pdf(N[1], rbind(c(0.2, 0.5))), pdf(N[1], c(0.2, 0.5)))
  expect_equal(pdf(N, c(0.2, 0.5)), drop(pdf(N, c(0.2, 0.5))))
  expect_equal(pdf(N_named, c(0.2, 0.5)), drop(pdf(N_named, c(0.2, 0.5))))
  expect_null(names(pdf(N[1], c(0.2, 0.5))))
  expect_null(names(pdf(N_named[1], c(0.2, 0.5))))

  ## length(d) = 1, at > 1, drop = FALSE
  expect_true(is.numeric(pdf(N[1], c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(pdf(N[1], c(0.2, 0.5), drop = FALSE)), c(1L, 2L))
  expect_equal(pdf(N[1], cbind(c(0.2, 0.5)), drop = FALSE), pdf(N[1], c(0.2, 0.5), drop = FALSE))
  expect_equal(pdf(N[1], rbind(c(0.2, 0.5)), drop = FALSE), pdf(N[1], c(0.2, 0.5), drop = FALSE))
  expect_equal(colnames(pdf(N[1], c(0.2, 0.5), drop = FALSE)), c("d_0.2", "d_0.5"))
  expect_null(rownames(pdf(N[1], c(0.2, 0.5), drop = FALSE)))
  expect_equal(colnames(pdf(N_named[1], c(0.2, 0.5), drop = FALSE)), c("d_0.2", "d_0.5"))
  expect_equal(rownames(pdf(N_named[1], c(0.2, 0.5), drop = FALSE)), "A")

  ## length(d) = at > 1, drop = TRUE
  expect_true(is.numeric(pdf(N[1:2], c(0.2, 0.5))))
  expect_null(dim(pdf(N[1:2], c(0.2, 0.5))))
  expect_length(pdf(N[1:2], c(0.2, 0.5)), 2)
  expect_equal(
    pdf(N[1:2], cbind(c(0.2, 0.5))),
    matrix(
      rbind(pdf(N[1], c(0.2, 0.5)), pdf(N[2], c(0.2, 0.5))),
      ncol = 2, dimnames = list(NULL, c("d_0.2", "d_0.5"))
    )
  )
  expect_equal(
    pdf(N[1:2], rbind(c(0.2, 0.5))),
    matrix(
      rbind(pdf(N[1], c(0.2, 0.5)), pdf(N[2], c(0.2, 0.5))),
      ncol = 2, dimnames = list(NULL, c("d_0.2", "d_0.5"))
    )
  )
  expect_equal(pdf(N, c(0.2, 0.5)), drop(pdf(N, c(0.2, 0.5))))
  expect_equal(pdf(N_named, c(0.2, 0.5)), drop(pdf(N_named, c(0.2, 0.5))))
  expect_null(names(pdf(N[1:2], c(0.2, 0.5))))
  expect_equal(names(pdf(N_named[1:2], c(0.2, 0.5))), c("A", "B"))

  ## length(d) = at > 1, drop = FALSE
  expect_true(is.numeric(pdf(N[1:2], c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(pdf(N[1:2], c(0.2, 0.5), drop = FALSE)), c(2L, 1L))
  expect_equal(
    pdf(N[1:2], cbind(c(0.2, 0.5)), drop = FALSE),
    matrix(
      rbind(pdf(N[1], c(0.2, 0.5)), pdf(N[2], c(0.2, 0.5))),
      ncol = 2, dimnames = list(NULL, c("d_0.2", "d_0.5"))
    )
  )
  expect_equal(
    pdf(N[1:2], rbind(c(0.2, 0.5)), drop = FALSE),
    matrix(
      rbind(pdf(N[1], c(0.2, 0.5)), pdf(N[2], c(0.2, 0.5))),
      ncol = 2, dimnames = list(NULL, c("d_0.2", "d_0.5"))
    )
  )
  expect_equal(colnames(pdf(N[1:2], c(0.2, 0.5), drop = FALSE)), "density")
  expect_null(rownames(pdf(N[1:2], c(0.2, 0.5), drop = FALSE)))
  expect_equal(colnames(pdf(N_named[1:2], c(0.2, 0.5), drop = FALSE)), "density")
  expect_equal(rownames(pdf(N_named[1:2], c(0.2, 0.5), drop = FALSE)), c("A", "B"))

  ## length(d) > 1, at > 1, drop = TRUE
  expect_true(is.numeric(pdf(N, c(0.2, 0.5))))
  expect_equal(dim(pdf(N, c(0.2, 0.5))), c(3L, 2L))
  expect_equal(pdf(N, cbind(c(0.2, 0.5))), pdf(N, c(0.2, 0.5)))
  expect_equal(pdf(N, rbind(c(0.2, 0.5))), pdf(N, c(0.2, 0.5)))
  expect_equal(pdf(N, c(0.2, 0.5)), drop(pdf(N, c(0.2, 0.5))))
  expect_equal(pdf(N_named, c(0.2, 0.5)), drop(pdf(N_named, c(0.2, 0.5))))
  expect_equal(colnames(pdf(N, c(0.2, 0.5))), c("d_0.2", "d_0.5"))
  expect_null(rownames(pdf(N, c(0.2, 0.5))))
  expect_equal(colnames(pdf(N_named, c(0.2, 0.5))), c("d_0.2", "d_0.5"))
  expect_equal(rownames(pdf(N_named, c(0.2, 0.5))), c("A", "B", "C"))

  ## length(d) > 1, at > 1, drop = FALSE
  expect_true(is.numeric(pdf(N, c(0.2, 0.5), drop = FALSE)))
  expect_equal(dim(pdf(N, c(0.2, 0.5), drop = FALSE)), c(3L, 2L))
  expect_equal(pdf(N, cbind(c(0.2, 0.5)), drop = FALSE), pdf(N, c(0.2, 0.5), drop = FALSE))
  expect_equal(pdf(N, rbind(c(0.2, 0.5)), drop = FALSE), pdf(N, c(0.2, 0.5), drop = FALSE))
  expect_equal(colnames(pdf(N, c(0.2, 0.5), drop = FALSE)), c("d_0.2", "d_0.5"))
  expect_null(rownames(pdf(N, c(0.2, 0.5), drop = FALSE)))
  expect_equal(colnames(pdf(N_named, c(0.2, 0.5), drop = FALSE)), c("d_0.2", "d_0.5"))
  expect_equal(rownames(pdf(N_named, c(0.2, 0.5), drop = FALSE)), c("A", "B", "C"))

  ## 'log_pdf', 'cdf', and 'quantile'
  expect_equal(colnames(log_pdf(N[1:2], c(0.2, 0.5), drop = FALSE)), "logLik")
  expect_equal(colnames(cdf(N[1:2], c(0.2, 0.5), drop = FALSE)), "probability")
  expect_equal(colnames(quantile(N[1:2], c(0.2, 0.5), drop = FALSE)), "quantile")
  expect_equal(colnames(log_pdf(N_named, c(0.2, 0.5), drop = FALSE)), c("l_0.2", "l_0.5"))
  expect_equal(colnames(cdf(N_named, c(0.2, 0.5), drop = FALSE)), c("p_0.2", "p_0.5"))
  expect_equal(colnames(quantile(N_named, c(0.2, 0.5), drop = FALSE)), c("q_0.2", "q_0.5"))
})

##get_deriv_names(p, which = NULL, expand = FALSE)
test_that("get_deriv_names() works as expected", {

  ## Default arguments
  expect_identical(formals(get_deriv_names),
    as.pairlist(alist(p =, which = NULL, expand = FALSE, check = TRUE)),
    info = "arguments/defaults not as expected")

  ## Testing incorrect use
  expect_error(get_deriv_names(p = TRUE),
    regex = "argument 'p' must be a character vector of length > 1L")
  expect_error(get_deriv_names(p = character()),
    regex = "argument 'p' must be a character vector of length > 1L")
  expect_error(get_deriv_names(p = ""),
    regex = "argument 'p' must be a character vector of length > 1L")

  expect_error(get_deriv_names("mu", which = TRUE),
    regex = "argument 'which' must be NULL or a character vector of length > 1L")
  expect_error(get_deriv_names("mu", which = character()),
    regex = "argument 'which' must be NULL or a character vector of length > 1L")
  expect_error(get_deriv_names("mu", which = ""),
    regex = "argument 'which' must be NULL or a character vector of length > 1L")

  expect_error(get_deriv_names("mu", expand = 'foo'),
    regex = "invalid argument type")

  ## Setting p, which = NULL, no expansion
  expect_silent(x <- get_deriv_names(letters[1:3]))
  expect_identical(x, setNames(letters[1:3], letters[1:3]))

  expect_silent(x <- get_deriv_names(letters[1:3], check = FALSE))
  expect_identical(x, setNames(letters[1:3], letters[1:3]))

  ## Setting p + which, checking return type/order
  expect_silent(x <- get_deriv_names(letters[1:3], which = "b"))
  expect_identical(x, c(b = "b"))
  expect_silent(x <- get_deriv_names(letters[1:3], which = c("c", "a")))
  expect_identical(x, c(c = "c", a = "a"))

  ## Checking 'expand = TRUE'
  expect_silent(x <- get_deriv_names(letters[1:3], expand = TRUE))
  tmp <- c("a" = "a", "a:b" = "b:a", "a:c" = "c:a", "a:b" = "a:b",
           "b" = "b", "b:c" = "c:b", "a:c" = "a:c", "b:c" = "b:c", "c" = "c")
  expect_identical(x, tmp)

  ## Checking 'expand = TRUE' with additional 'which' argument (& correct order)
  expect_silent(x <- get_deriv_names(letters[1:3], which = c("b", "b:c", "c:b", "a", "b"), expand = TRUE))
  expect_identical(x, c("b" = "b", "b:c" = "b:c", "b:c" = "c:b", "a" = "a", "b" = "b"))

  ## If which == p it should be returned as a named version of itself (testing shortcut)
  expect_silent(x <- get_deriv_names(letters[1:3], which = letters[1:3]))
  expect_identical(x, setNames(letters[1:3], letters[1:3]))

  ## Ensure it respects the order of 'which'
  expect_silent(x <- get_deriv_names(letters[1:3], which = letters[3:1]))
  expect_identical(x, setNames(letters[3:1], letters[3:1]))

})


##max_length(..., check = TRUE)
test_that("max_length() works as expected", {

  ## Default arguments
  expect_identical(formals(max_length),
    as.pairlist(alist(... =, check = TRUE)),
    info = "arguments/defaults not as expected")

  ## No non-negative input arguments: expecting warning (from max()) and -Inf as return
  expect_warning(x <- max_length())
  expect_identical(x, -Inf)
  expect_warning(x <- max_length(check = FALSE))
  expect_identical(x, -Inf)

  ## All inputs of length 1: all fine
  expect_identical(expect_silent(max_length(1, 2, 3, 4)), 1L)
  expect_identical(expect_silent(max_length(1, 2, 3, 4, check = FALSE)), 1L)
  expect_identical(expect_silent(max_length(a = 1, b = 2, c = 3, d = 4)), 1L)
  expect_identical(expect_silent(max_length(a = 1, b = 2, c = 3, d = 4, check = FALSE)), 1L)

  ## All length 5, all fine
  x <- 1:5
  expect_identical(expect_silent(max_length(x, x, x, x)), 5L)
  expect_identical(expect_silent(max_length(x, x, x, x, check = FALSE)), 5L)
  expect_identical(expect_silent(max_length(a = x, b = x, c = x, d = x)), 5L)
  expect_identical(expect_silent(max_length(a = x, x, c = x, x)), 5L) # mixed named/unnamed
  expect_identical(expect_silent(max_length(a = x, b = x, c = x, d = x, check = FALSE)), 5L)
  expect_identical(expect_silent(max_length(x, b = x, c = x, x, check = FALSE)), 5L) # mixed named/unnamed

  ## Mixed length, but check = FALSE (should run silently, returning max length)
  expect_identical(expect_silent(max_length(1:10, 1:3, 1:6, check = FALSE)), 10L)
  expect_identical(expect_silent(max_length(a = 1:10, b = 1:3, c = 1:6, check = FALSE)), 10L)
  expect_identical(expect_silent(max_length(a = 1:10, 1:3, 1:6, check = FALSE)), 10L) # mixed named/unnamed

  ## Check = TRUE but lengths not matching (i.e., not all of length L or N),
  ## expecting error and appropriate error message
  expect_error(max_length(1:10, 1:3),
    regex = "parameter lengths do not match.*got lengths\\: 1:10 = 10, 1\\:3 = 3") # unnamed
  expect_error(max_length(foo = 1:10, bar = 1:3),
    regex = "parameter lengths do not match.*got lengths\\: foo = 10, bar = 3") # named
  expect_error(max_length(foo = 1:10, 1:3),
    regex = "parameter lengths do not match.*got lengths\\: foo = 10, 1\\:3 = 3") # mixed named/unnamed

})


devtools::load_all("../")
d <- Normal(1:3) |> setNames(letters[1:3])
testfun <- function(x, d) switch(x, "mu" = rep(1, length(x), "sigma" = rep(2, length(x)), rep(3, length(x))))
apply_deriv(d, testfun, "mu")

devtools::load_all("../")
apply_deriv(d, testfun, c("mu"= "mu"))
apply_deriv(d, testfun, c("sigma"= "mu"))
apply_deriv(d, testfun, c("mu" = "fooo"), drop = FALSE)

devtools::load_all("../")
apply_deriv(unname(d), testfun, c("mu"= "mu"))
apply_deriv(unname(d), testfun, c("sigma" = "mu"))
apply_deriv(unname(d), testfun, c("mu" = "fooo"), drop = FALSE)

apply_deriv(d, testfun, c("foo"= "mu"), drop = FALSE) # ERR

testfun('mu')

devtools::load_all("../")
apply_deriv(d, testfun, c("mu"= "mu", "sigma" = "sigma"))
apply_deriv(d, testfun, c("mu" = "foo", "sigma" = "bar"))
apply_deriv(d, testfun, c("mu" = "foo", "sigma" = "bar"), drop = FALSE) # No effect




