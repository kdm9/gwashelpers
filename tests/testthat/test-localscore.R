# Tests for the local score / Lindley process implementation.
# All expected values are derived clean-room from the paper formulae:
#   Fariello et al. (2017) doi:10.1111/mec.14141, Equations 1-4.

# ---------------------------------------------------------------------------
# 4a. Lindley process -- paper Eq. 1 (hand-calculable)
# ---------------------------------------------------------------------------
test_that("4a: Lindley process matches paper Eq. 1 (hand-computed)", {
    scores = c(0.5, -0.3, 0.8, -1.0, 0.2, 0.3)
    # h1 = max(0, 0 + 0.5)       = 0.5
    # h2 = max(0, 0.5 + (-0.3))  = 0.2
    # h3 = max(0, 0.2 + 0.8)     = 1.0
    # h4 = max(0, 1.0 + (-1.0))  = 0.0
    # h5 = max(0, 0.0 + 0.2)     = 0.2
    # h6 = max(0, 0.2 + 0.3)     = 0.5
    expect_equal(lindley(scores), c(0.5, 0.2, 1.0, 0.0, 0.2, 0.5))
})

# ---------------------------------------------------------------------------
# 4b. Lindley process -- non-negativity invariant
# ---------------------------------------------------------------------------
test_that("4b: Lindley process is always non-negative", {
    # All-positive input: non-decreasing, positive
    expect_true(all(lindley(c(1, 2, 3)) >= 0))
    # All-negative input: should be all zeros (process resets immediately)
    expect_equal(lindley(c(-1, -2, -3, -0.5)), c(0, 0, 0, 0))
    # Mixed input: still non-negative
    set.seed(99)
    random_scores = rnorm(200)
    expect_true(all(lindley(random_scores) >= 0))
})

# ---------------------------------------------------------------------------
# 4c. Score transformation (paper section 2.1.2)
# ---------------------------------------------------------------------------
# p = c(0.001, 0.01, 0.1, 0.5), xi = 2
# -log10(p) = c(3, 2, 1, ~0.301)
# X_m = c(1, 0, -1, ~-1.699)  -- only p < 0.01 gives X_m > 0
test_that("4c: Score column equals -log10(p) - xi", {
    p = c(0.001, 0.01, 0.1, 0.5)
    xi = 2
    result = lindley_p(p, xi = xi)
    expect_equal(result$score, -log10(p) - xi, tolerance = 1e-10)
    # Only p < 10^(-xi) = 0.01 gives positive score
    expect_true(result$score[1] > 0)   # p = 0.001 < 0.01
    expect_equal(result$score[2], 0, tolerance = 1e-10)  # p = 0.01 exactly -> score = 0
    expect_true(result$score[3] < 0)   # p = 0.1 > 0.01
})

test_that("4c: Score transformation with xi=1 shifts cut-off to p < 0.1", {
    p = c(0.001, 0.01, 0.1, 0.5)
    result = lindley_p(p, xi = 1)
    # With xi=1, cut-off is p < 0.1 = 10^(-1)
    expect_true(result$score[1] > 0)   # p = 0.001 < 0.1
    expect_true(result$score[2] > 0)   # p = 0.01 < 0.1
    expect_equal(result$score[3], 0, tolerance = 1e-10)  # p = 0.1 exactly
    expect_true(result$score[4] < 0)   # p = 0.5 > 0.1
})

# ---------------------------------------------------------------------------
# 4d. Autocorrelation -- perfectly correlated (monotone) sequence
# ---------------------------------------------------------------------------
# x = 1:5: lag-1 Pearson correlation = 1.0
test_that("4d: autocor returns 1.0 for a perfectly correlated sequence", {
    expect_equal(gwashelpers:::autocor(1:5), 1.0, tolerance = 1e-10)
    expect_equal(gwashelpers:::autocor(as.numeric(1:20)), 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4e. Autocorrelation -- perfectly anti-correlated sequence
# ---------------------------------------------------------------------------
# x = c(1, -1, 1, -1, 1): lag-1 correlation = -1.0, abs = 1.0
test_that("4e: autocor returns 1.0 for a perfectly anti-correlated sequence", {
    expect_equal(gwashelpers:::autocor(c(1, -1, 1, -1, 1)), 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4f. Autocorrelation -- AR(1) simulation recovers phi
# ---------------------------------------------------------------------------
# AR(1) with phi=0.3 on raw values (not pnorm-transformed) should give
# autocor (Pearson lag-1) approximately equal to phi=0.3 (within tolerance 0.1).
test_that("4f: autocor recovers phi from an AR(1) process", {
    set.seed(42)
    n = 500
    x = numeric(n)
    x[1] = rnorm(1)
    for (i in 2:n) x[i] = 0.3 * x[i - 1] + rnorm(1)
    rho = gwashelpers:::autocor(x)
    expect_true(rho > 0.25 && rho < 0.35,
        info = paste("Expected rho approx 0.3, got", rho))
})

# ---------------------------------------------------------------------------
# 4g. Threshold formula -- hand-computed from paper Eq. 4
# ---------------------------------------------------------------------------
# Coefficient table (a.poly.coef, b.poly.coef), R 1-indexed columns = xi:
#
#   a.poly.coef:
#     rho^3 row: c(-5.5, 2.47, 2.04, 0.22)  -- xi=1,2,3,4
#     rho^2 row: c(6.76, -4.16, -5.76, -4.08)
#     rho   row: c(-5.66, -1.82, 1.04, 1.16)
#     const row: c(-2.51, -4.58, -6.95, -9.16)
#
#   b.poly.coef:
#     rho^2 row: c(-1.22, 0.37, 2.55, 3.45)
#     rho   row: c(3.17, 2.14, -0.02, -0.98)
#     const row: c(-1.99, -2.35, -2.31, -2.33)
#
# For xi=2 (R column 2):
#   a coefficients: rho^3=2.47, rho^2=-4.16, rho=-1.82, intercept=-4.58
#   b coefficients: rho^2=0.37, rho=2.14, intercept=-2.35
#
# At rho=0, xi=2, n=10000, alpha=0.05:
#   a0 = log(10000) + (-4.58) = 9.21034 - 4.58 = 4.63034
#   b0 = -2.35
#   log(-log(1-0.05)) = log(-log(0.95)) = log(0.051293) = -2.97015
#   thres = (-2.97015 - 4.63034) / (-2.35) = -7.60049 / -2.35 approx 3.234
test_that("4g-A: lindley_thresh near rho=0 matches analytical value (xi=2)", {
    # Uniform random p-values are nearly uncorrelated (autocor approx 0)
    set.seed(999)
    p = runif(10000)
    rho = gwashelpers:::autocor(p)
    expect_true(abs(rho) < 0.05,
        info = paste("autocor =", rho, "-- expected approx 0"))

    # Hand-computed threshold with rho=0, n=10000, xi=2, alpha=0.05
    # Using xi=2 column: a.poly.coef[,2]: 2.47,-4.16,-1.82,-4.58
    #                    b.poly.coef[,2]: 0.37, 2.14, -2.35
    a0 = log(10000) + (-4.58)  # polynomial terms vanish at rho=0
    b0 = -2.35
    expected0 = (log(-log(1 - 0.05)) - a0) / b0  # approx 3.234

    actual = lindley_thresh(p, xi = 2, alpha = 0.05)
    expect_equal(actual, expected0, tolerance = 0.1)
})

test_that("4g-B: lindley_thresh exactly reproduces the polynomial formula", {
    # The function output must match a manual re-implementation of the formula.
    # Uses the xi=2 column of the coefficient matrices (R 1-indexed column 2).
    set.seed(42)
    p = runif(10000)
    rho = gwashelpers:::autocor(p)
    nmarker = length(p)

    # xi=2 coefficients: a.poly.coef[,2] and b.poly.coef[,2]
    a = log(nmarker) + 2.47 * rho^3 + (-4.16) * rho^2 + (-1.82) * rho + (-4.58)
    b = 0.37 * rho^2 + 2.14 * rho + (-2.35)
    expected = (log(-log(1 - 0.05)) - a) / b

    actual = lindley_thresh(p, xi = 2, alpha = 0.05)
    expect_equal(actual, expected, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4h. Threshold scales with nmarker (log-linear)
# ---------------------------------------------------------------------------
# a = log(M) + poly(rho): larger M -> larger a -> larger threshold.
test_that("4h: threshold increases with number of markers", {
    # Use uncorrelated p-values so autocor is near 0 for both
    set.seed(123)
    p_small = runif(1000)
    p_large = runif(100000)

    thres_small = lindley_thresh(p_small, xi = 2, alpha = 0.05)
    thres_large = lindley_thresh(p_large, xi = 2, alpha = 0.05)

    expect_true(thres_large > thres_small,
        info = paste("thres_small =", thres_small, "; thres_large =", thres_large))
})

# ---------------------------------------------------------------------------
# 4i. xi validity check
# ---------------------------------------------------------------------------
test_that("4i: lindley_thresh returns NULL for invalid xi", {
    set.seed(1)
    p = runif(100)
    expect_null(suppressMessages(lindley_thresh(p, xi = 5)))
    expect_null(suppressMessages(lindley_thresh(p, xi = 0)))
})

# ---------------------------------------------------------------------------
# 4j. windows() -- paper section 2.1.1 excursion interval labelling
# ---------------------------------------------------------------------------
test_that("4j: windows() labels consecutive TRUE runs correctly", {
    state = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
    expected = c(NA_real_, 1, 1, NA_real_, 2, NA_real_, 3, 3)
    result = windows(state)
    expect_equal(result, expected)
})

test_that("4j: windows() returns all NA for all-FALSE input", {
    expect_equal(windows(c(FALSE, FALSE, FALSE)), c(NA_real_, NA_real_, NA_real_))
})

test_that("4j: windows() labels a single TRUE run as window 1", {
    expect_equal(windows(c(TRUE, TRUE, TRUE)), c(1, 1, 1))
})

# ---------------------------------------------------------------------------
# 4k. lindley_p() output structure
# ---------------------------------------------------------------------------
test_that("4k: lindley_p() returns a data.frame with correct structure", {
    set.seed(42)
    p = runif(50, 0.001, 1)
    result = lindley_p(p, xi = 2)

    expect_s3_class(result, "data.frame")
    expect_named(result, c("score", "lindley_score", "lindley_thresh"))
    expect_equal(nrow(result), length(p))
    # Threshold is a scalar replicated to every row
    expect_true(length(unique(result$lindley_thresh)) == 1)
    # Lindley scores are non-negative
    expect_true(all(result$lindley_score >= 0))
})

test_that("4k: lindley_p() rejects out-of-range p-values", {
    expect_error(lindley_p(c(0.5, 1.5), xi = 2))
    expect_error(lindley_p(c(-0.1, 0.5), xi = 2))
})

# ---------------------------------------------------------------------------
# 4l. Integration test -- detectable peak in middle of sequence
# ---------------------------------------------------------------------------
# p-values: 10 background markers (p=0.3), 5 very significant (p=0.0001),
# 10 background markers. The Lindley score should peak in the signal cluster.
test_that("4l: lindley_p() detects a clear signal cluster (integration test)", {
    p = c(rep(0.3, 10), rep(0.0001, 5), rep(0.3, 10))
    result = lindley_p(p, xi = 1)

    # Peak of Lindley process should be in the middle cluster (positions 11-15)
    peak_pos = which.max(result$lindley_score)
    expect_true(peak_pos >= 11 && peak_pos <= 15,
        info = paste("Peak at position", peak_pos))

    # Lindley score at peak should exceed scores at the flanks
    expect_true(result$lindley_score[peak_pos] > result$lindley_score[1])
    expect_true(result$lindley_score[peak_pos] > result$lindley_score[25])
})
