# Tests for the local score / Lindley process implementation.
# All expected values are derived clean-room from the paper formulae:
#   Fariello et al. (2017) doi:10.1111/mec.14141, Equations 1–4.

# ---------------------------------------------------------------------------
# 4a. Lindley process — paper Eq. 1 (hand-calculable)
# ---------------------------------------------------------------------------
# scores = c(0.5, -0.3, 0.8, -1.0, 0.2, 0.3)
# h1 = max(0, 0 + 0.5)       = 0.5
# h2 = max(0, 0.5 + (-0.3))  = 0.2
# h3 = max(0, 0.2 + 0.8)     = 1.0
# h4 = max(0, 1.0 + (-1.0))  = 0.0
# h5 = max(0, 0.0 + 0.2)     = 0.2
# h6 = max(0, 0.2 + 0.3)     = 0.5
test_that("4a: Lindley process matches paper Eq. 1 (hand-computed)", {
    scores <- c(0.5, -0.3, 0.8, -1.0, 0.2, 0.3)
    expect_equal(lindley(scores), c(0.5, 0.2, 1.0, 0.0, 0.2, 0.5))
    expect_equal(max(lindley(scores)), 1.0)
})

# ---------------------------------------------------------------------------
# 4b. Lindley process — non-negativity invariant
# ---------------------------------------------------------------------------
test_that("4b: Lindley process is always non-negative", {
    # All-positive input: non-decreasing, positive
    expect_true(all(lindley(c(1, 2, 3)) >= 0))
    # All-negative input: should be all zeros (process resets immediately)
    expect_equal(lindley(c(-1, -2, -3, -0.5)), c(0, 0, 0, 0))
    # Mixed input: still non-negative
    set.seed(99)
    random_scores <- rnorm(200)
    expect_true(all(lindley(random_scores) >= 0))
})

# ---------------------------------------------------------------------------
# 4c. Score transformation — paper §2.1.2
# ---------------------------------------------------------------------------
# p = c(0.001, 0.01, 0.1, 0.5), xi = 2
# -log10(p) = c(3, 2, 1, ~0.301)
# X_m = c(1, 0, -1, ~-1.699)  — only p < 0.01 gives X_m > 0
test_that("4c: Score column equals -log10(p) - xi", {
    p <- c(0.001, 0.01, 0.1, 0.5)
    xi <- 2
    result <- lindley_p(p, xi = xi)
    expect_equal(result$score, -log10(p) - xi, tolerance = 1e-10)
    # Only p < 10^(-xi) = 0.01 gives positive score
    expect_true(result$score[1] > 0)   # p = 0.001 < 0.01
    expect_equal(result$score[2], 0, tolerance = 1e-10)  # p = 0.01 exactly -> score = 0
    expect_true(result$score[3] < 0)   # p = 0.1 > 0.01
})

test_that("4c: Score transformation with xi=1 shifts cut-off to p < 0.1", {
    p <- c(0.001, 0.01, 0.1, 0.5)
    result <- lindley_p(p, xi = 1)
    # With xi=1, cut-off is p < 0.1 = 10^(-1)
    expect_true(result$score[1] > 0)   # p = 0.001 < 0.1
    expect_true(result$score[2] > 0)   # p = 0.01 < 0.1
    expect_equal(result$score[3], 0, tolerance = 1e-10)  # p = 0.1 exactly
    expect_true(result$score[4] < 0)   # p = 0.5 > 0.1
})

# ---------------------------------------------------------------------------
# 4d. Autocorrelation — perfectly correlated (monotone) sequence
# ---------------------------------------------------------------------------
# x = 1:5: lag-1 Pearson correlation = 1.0
test_that("4d: autocor returns 1.0 for a perfectly correlated sequence", {
    expect_equal(gwashelpers:::autocor(1:5), 1.0, tolerance = 1e-10)
    expect_equal(gwashelpers:::autocor(as.numeric(1:20)), 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4e. Autocorrelation — perfectly anti-correlated sequence
# ---------------------------------------------------------------------------
# x = c(1, -1, 1, -1, 1): lag-1 correlation = -1.0, abs = 1.0
test_that("4e: autocor returns 1.0 for a perfectly anti-correlated sequence", {
    expect_equal(gwashelpers:::autocor(c(1, -1, 1, -1, 1)), 1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4f. Autocorrelation — AR(1) simulation recovers phi
# ---------------------------------------------------------------------------
# AR(1) with phi=0.3 should give autocor ≈ 0.3 (within tolerance 0.1)
test_that("4f: autocor recovers phi from an AR(1) process", {
    set.seed(42)
    n <- 500
    x <- numeric(n)
    x[1] <- rnorm(1)
    for (i in 2:n) x[i] <- 0.3 * x[i - 1] + rnorm(1)
    rho <- gwashelpers:::autocor(x)
    expect_true(rho > 0.2 && rho < 0.4,
        info = paste("Expected rho ≈ 0.3, got", rho))
})

# ---------------------------------------------------------------------------
# 4g. Threshold formula — hand-computed from paper Eq. 4
# ---------------------------------------------------------------------------
# Parameters: xi=2, cor≈0.25, nmarker=10000, alpha=0.05
#
# Coefficient matrices (xi=2 is column 2):
#   a_poly: c(2.04, -5.76, 1.04, -6.95)  [rho^3, rho^2, rho, intercept]
#   b_poly: c(2.55, -0.02, -2.31)         [rho^2, rho, intercept]
#
# a = log(10000) + 2.04*(0.25^3) + (-5.76)*(0.25^2) + 1.04*0.25 + (-6.95)
#   = 9.21034 + 0.03188 - 0.36 + 0.26 - 6.95
#   = 2.19222
#
# b = 2.55*(0.25^2) + (-0.02)*0.25 + (-2.31)
#   = 0.159375 - 0.005 - 2.31
#   = -2.15563
#
# log(-log(1 - 0.05)) = log(-log(0.95)) = log(0.051293) = -2.97015
#
# thres = (-2.97015 - 2.19222) / (-2.15563)
#       = -5.16237 / -2.15563
#       ≈ 2.395
test_that("4g: lindley_thresh matches hand-computed Gumbel threshold", {
    # Build an AR(1) p-value vector targeting autocor ≈ 0.25 with 10000 markers.
    # Use pnorm() to map AR(1) realisations into (0,1) as p-value surrogates.
    set.seed(7)
    n <- 10000
    phi <- 0.25  # AR(1) coefficient chosen so autocor ≈ 0.25
    x <- numeric(n)
    x[1] <- rnorm(1)
    for (i in 2:n) x[i] <- phi * x[i - 1] + sqrt(1 - phi^2) * rnorm(1)
    pval <- pnorm(x)

    # Verify the autocorrelation is close enough to 0.25 for the formula
    rho <- gwashelpers:::autocor(pval)
    expect_true(abs(rho - 0.25) < 0.05,
        info = paste("autocor =", rho, "— expected ≈ 0.25"))

    # Compute threshold and compare with hand-computed value
    thres <- lindley_thresh(pval, xi = 2, alpha = 0.05)
    expect_equal(thres, 2.395, tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# 4h. Threshold scales with nmarker (log-linear)
# ---------------------------------------------------------------------------
# a = log(M) + poly(rho): larger M → larger a → larger threshold.
test_that("4h: threshold increases with number of markers", {
    # Use uncorrelated p-values so autocor is near 0 for both
    set.seed(123)
    p_small <- runif(1000)
    p_large <- runif(100000)

    thres_small <- lindley_thresh(p_small, xi = 2, alpha = 0.05)
    thres_large <- lindley_thresh(p_large, xi = 2, alpha = 0.05)

    expect_true(thres_large > thres_small,
        info = paste("thres_small =", thres_small, "; thres_large =", thres_large))
})

# ---------------------------------------------------------------------------
# 4i. xi validity check
# ---------------------------------------------------------------------------
test_that("4i: lindley_thresh returns NULL for invalid xi", {
    set.seed(1)
    p <- runif(100)
    expect_null(lindley_thresh(p, xi = 5))
    expect_null(lindley_thresh(p, xi = 0))
    expect_null(suppressMessages(lindley_thresh(p, xi = 5)))
})

# ---------------------------------------------------------------------------
# 4j. windows() — paper §2.1.1 excursion interval labelling
# ---------------------------------------------------------------------------
test_that("4j: windows() labels consecutive TRUE runs correctly", {
    state <- c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
    expected <- c(NA, 1, 1, NA, 2, NA, 3, 3)
    result <- windows(state)
    expect_equal(result, expected)
})

test_that("4j: windows() returns all NA for all-FALSE input", {
    expect_equal(windows(c(FALSE, FALSE, FALSE)), c(NA, NA, NA))
})

test_that("4j: windows() labels a single TRUE run as window 1", {
    expect_equal(windows(c(TRUE, TRUE, TRUE)), c(1, 1, 1))
})

# ---------------------------------------------------------------------------
# 4k. lindley_p() output structure
# ---------------------------------------------------------------------------
test_that("4k: lindley_p() returns a data.frame with correct structure", {
    set.seed(42)
    p <- runif(50, 0.001, 1)
    result <- lindley_p(p, xi = 2)

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
# 4l. Integration test — detectable peak in middle of sequence
# ---------------------------------------------------------------------------
# p-values: 10 background markers (p≈0.3), 5 very significant (p=0.0001),
# 10 background markers. The Lindley score should peak in the signal cluster.
test_that("4l: lindley_p() detects a clear signal cluster (integration test)", {
    p <- c(rep(0.3, 10), rep(0.0001, 5), rep(0.3, 10))
    result <- lindley_p(p, xi = 1)

    # Peak of Lindley process should be in the middle cluster (positions 11–15)
    peak_pos <- which.max(result$lindley_score)
    expect_true(peak_pos >= 11 && peak_pos <= 15,
        info = paste("Peak at position", peak_pos))

    # Lindley score at peak should exceed scores at the flanks
    expect_true(result$lindley_score[peak_pos] > result$lindley_score[1])
    expect_true(result$lindley_score[peak_pos] > result$lindley_score[25])
})
