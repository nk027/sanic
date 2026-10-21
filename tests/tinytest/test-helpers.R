library(tinytest)
library(Matrix)

S <- matrix(c(2, 1, 1, 2), 2)
G <- matrix(c(1, 2, 3, 4), 2)
S_sparse <- sparsify(S)

# Matrix predicates ---------------------------------------------------------
expect_true(sanic:::is_sparse(S_sparse))
expect_false(sanic:::is_sparse(S))
expect_true(sanic:::is_square(S))
expect_false(sanic:::is_square(matrix(1, 2, 3)))
expect_true(sanic:::is_symmetric(S))
expect_true(sanic:::is_symmetric(S_sparse))
expect_false(sanic:::is_symmetric(G))
expect_false(sanic:::is_symmetric(matrix(1, 2, 3)))
expect_error(sanic:::is_symmetric(1:4), "matrix")

set.seed(1)
expect_true(sanic:::maybe_symmetric(S))
expect_false(sanic:::maybe_symmetric(G, n_compared = 1000L))
expect_false(sanic:::maybe_symmetric(matrix(1, 2, 3)))
expect_true(is.logical(sanic:::maybe_symmetric(G, checks = FALSE)))
expect_error(sanic:::maybe_symmetric(1:4), "matrix")

# Numeric and integer checks ------------------------------------------------
expect_equal(sanic:::num_check(1), 1)
expect_equal(sanic:::num_check(1L), 1)
expect_equal(sanic:::num_check(2, min = 1, max = 3, fun = sqrt), sqrt(2))
expect_error(sanic:::num_check("1"), "check")
expect_error(sanic:::num_check(c(1, 2)), "check")
expect_error(sanic:::num_check(-1, min = 0), "check")
expect_error(sanic:::num_check(4, max = 3), "check")
expect_equal(sanic:::num_check(1, msg = "custom message"), 1)

expect_identical(sanic:::int_check(2), 2L)
expect_identical(sanic:::int_check(2.9), 2L)
expect_error(sanic:::int_check(-1), "parameters")
expect_error(sanic:::int_check(1, min = 2, msg = "too small"), "too small")
