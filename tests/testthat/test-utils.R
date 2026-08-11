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
  expect_length(print.distribution(n), 2L)
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
