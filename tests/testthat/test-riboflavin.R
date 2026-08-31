test_that("riboflavin dataset has the expected structure", {
  expect_type(riboflavin, "list")
  expect_setequal(names(riboflavin), c("A", "b", "group"))

  expect_true(is.matrix(riboflavin$A))
  expect_equal(dim(riboflavin$A), c(71, 1199))
  expect_true(all(is.finite(riboflavin$A)))
  expect_false(is.null(rownames(riboflavin$A)))
  expect_false(is.null(colnames(riboflavin$A)))

  expect_length(riboflavin$b, 71)
  expect_true(all(is.finite(riboflavin$b)))

  expect_length(riboflavin$group, 1199)
  expect_length(unique(riboflavin$group), 36)
  expect_true("GO:0006766" %in% riboflavin$group) # vitamin metabolic process
})

test_that("riboflavin fits cleanly through the public API", {
  fit <- sglssnal(riboflavin$A, riboflavin$b, riboflavin$group, lambda = 0.05)
  expect_equal(dim(fit$x), c(1199, 1))
  expect_true(all(is.finite(fit$x)))

  # also exercise the auto lambda path with dense A end-to-end
  fit_auto <- sglssnal(riboflavin$A, riboflavin$b, riboflavin$group, nlambda = 5)
  expect_equal(dim(fit_auto$x), c(1199, 5))
})
