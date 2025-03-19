# A small "mock" Y: D=3 features, N=4 samples
mockY <- matrix(c(
  10, 12, 15, 20,
  5,   5,  7,  3,
  1,   2,  2,  2
), nrow = 3, ncol = 4, byrow = TRUE)

# Basic X: a design matrix for N=4 samples, p=2 covariates
mockX <- matrix(rnorm(4 * 2), nrow = 4, ncol = 2)
mock_nsample <- 5
set.seed(123)
mock_logWpara <- array(rnorm(3 * 4 * mock_nsample, mean=2, sd=1), 
                       dim = c(3, 4, mock_nsample))
mock_logWpara <- log2(exp(mock_logWpara))
mock_externalscale <- c(1.0, 2.0, 0.5)

test_that("GAMMA works without externalscale (design matrix approach)", {
  user_data <- list(X = mockX, gamma = 0.5, externalscale = NULL)
  result <- GAMMA(mock_logWpara, user_data)
  expect_equal(dim(result), c(3, 5))
  expect_true(all(!is.na(result)))
})

test_that("GAMMA requires X if externalscale is NULL", {
  # If externalscale is NULL and X is missing, we expect an error
  user_data <- list(X = NULL, gamma = 0.5, externalscale = NULL)
  expect_error(GAMMA(mock_logWpara, user_data),
               regexp = "must provide 'X'")
})

test_that("GAMMA works with externalscale", {
  user_data <- list(externalscale = mock_externalscale, gamma = 0.2)
  result <- GAMMA(mock_logWpara, user_data)
  expect_equal(dim(result), c(3, 5))
})

test_that("GAMMA errors if gamma length != externalscale length", {
  user_data <- list(externalscale = mock_externalscale, 
                    gamma = c(0.2, 0.3))  # length 2
  # externalscale has length 3 => mismatch => error
  expect_error(GAMMA(mock_logWpara, user_data),
               regexp = "Length of 'gamma' does not match length of 'externalscale'")
})

test_that("Case B: gamma=0 gives deterministic result", {
  user_data <- list(externalscale = mock_externalscale, gamma = 0)
  result <- GAMMA(mock_logWpara, user_data)
  expect_equal(result, matrix(mock_externalscale, nrow = 3, ncol = 5))
})

test_that("Case B: external scale with gamma vector matching length", {
  user_data <- list(externalscale = mock_externalscale, gamma = c(0.1, 0.2, 0.3))
  result <- GAMMA(mock_logWpara, user_data)
  expect_equal(dim(result), c(3, 5))
})

test_that("GAMMA warns for small log2 values", {
  small_logWpara <- mock_logWpara - 20  # Shift values down
  expect_warning(
    GAMMA(small_logWpara, list(X = mockX)),
    "outside typical log2 range"
  )
})

test_that("Case A baseline offset is correct", {
  user_data <- list(X = mockX, gamma = 0)
  result <- GAMMA(mock_logWpara, user_data)
  expected_baseline <- -colMeans(mock_logWpara, dims = 1)
  expect_equal(result, expected_baseline)  # Gamma=0 => no randomness
})

test_that("Auto-conversion from loge to log2 works", {
  loge_values <- matrix(1, nrow = 3, ncol = 5)  # loge(1) = 0
  log2_values <- log2(exp(loge_values))          # Should be 0
  user_data <- list(externalscale = rep(1, 3), gamma = 0)
  result <- suppressWarnings(GAMMA(loge_values, user_data))
  expect_equal(result, log2_values)
})

test_that("Case A gamma=0 gives deterministic result", {
  user_data <- list(X = mockX, gamma = 0)
  result <- GAMMA(mock_logWpara, user_data)
  expected <- -colMeans(mock_logWpara, dims = 1)  # No Lambdaperp term
  expect_equal(result, expected)
})

test_that("GAMMA warns if values are out of typical log2 range", {
  # Create artificially huge logWpara so that the result might exceed log2 range
  big_logWpara <- mock_logWpara + 200  # shift way up
  
  # We'll capture warnings with expect_warning
  expect_warning(
    result <- GAMMA(big_logWpara, list(X = mockX)),
    regexp = "outside typical log2 range"
  )
  
  # Also expect a second warning about auto-conversion if max_val > 100
  expect_warning(
    result2 <- GAMMA(big_logWpara, list(X = mockX)),
    regexp = "auto-converting from log\\(e\\) to log2"
  )
})

# -----------------------------
# Now the wrapper functions:
# -----------------------------

test_that("default_GAMMA works (no external scale)", {
  result <- default_GAMMA(mock_logWpara, mockX, gamma = 0.5)
  expect_equal(dim(result), c(3, 5))
})

test_that("external_GAMMA works (with external scale)", {
  result <- external_GAMMA(mock_logWpara, mock_externalscale, gamma = 0.3)
  expect_equal(dim(result), c(3, 5))
})

test_that("sensitivity_GAMMA runs multiple gamma values", {
  gamma_vals <- c(0.1, 0.5, 1.0)
  results_list <- sensitivity_GAMMA(mock_logWpara, mockX, gamma_vals)
  expect_equal(length(results_list), length(gamma_vals))
  lapply(results_list, function(r) {
    expect_equal(dim(r), c(3, 5))
  })
})

