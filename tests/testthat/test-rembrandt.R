test_that("rembrandt dataset has the expected structure", {
  expect_type(rembrandt, "list")
  expect_setequal(names(rembrandt), c("A", "b", "group", "disease_type"))

  expect_true(inherits(rembrandt$A, "sparseMatrix"))
  expect_equal(dim(rembrandt$A), c(172, 811))
  expect_true(all(is.finite(rembrandt$A@x)))

  expect_length(rembrandt$b, 172)
  expect_true(all(is.finite(rembrandt$b)))

  expect_length(rembrandt$group, 811)
  expect_setequal(unique(rembrandt$group), as.character(1:22))

  expect_length(rembrandt$disease_type, 172)
  expect_true(is.factor(rembrandt$disease_type))
})

test_that("rembrandt fits cleanly through the public API", {
  fit <- sglssnal(rembrandt$A, rembrandt$b, rembrandt$group, lambda = 0.5)
  expect_equal(dim(fit$x), c(811, 1))
  expect_true(all(is.finite(fit$x@x)))

  # also exercise the auto lambda path with sparse A end-to-end
  fit_auto <- sglssnal(rembrandt$A, rembrandt$b, rembrandt$group, nlambda = 5)
  expect_equal(dim(fit_auto$x), c(811, 5))
})
