# Tests for the local score / Lindley process implementation.
# All expected values are derived clean-room from the paper formulae:
#   Fariello et al. (2017) doi:10.1111/mec.14141, Equations 1-4.

###################
# Lindley process #
###################

test_that("Lindley process matches paper Eq. 1 (hand-computed)", {
    scores = c(0.5, -0.3, 0.8, -1.0, 0.2, 0.3)
    # h1 = max(0, 0 + 0.5)       = 0.5
    # h2 = max(0, 0.5 + (-0.3))  = 0.2
    # h3 = max(0, 0.2 + 0.8)     = 1.0
    # h4 = max(0, 1.0 + (-1.0))  = 0.0
    # h5 = max(0, 0.0 + 0.2)     = 0.2
    # h6 = max(0, 0.2 + 0.3)     = 0.5
    expect_equal(lindley(scores), c(0.5, 0.2, 1.0, 0.0, 0.2, 0.5))
})

test_that("Lindley process is always non-negative", {
    # All-positive input: non-decreasing, positive
    expect_true(all(lindley(c(1, 2, 3)) >= 0))
    # All-negative input: should be all zeros (process resets immediately)
    expect_equal(lindley(c(-1, -2, -3, -0.5)), c(0, 0, 0, 0))
    # Mixed input: still non-negative
    set.seed(99)
    random_scores = rnorm(200)
    expect_true(all(lindley(random_scores) >= 0))
})

test_that("Score equals -log10(p) - xi", {
    # p = c(0.001, 0.01, 0.1, 0.5), xi = 2
    # -log10(p) = c(3, 2, 1, ~0.301)
    # X_m = c(1, 0, -1, ~-1.699)  -- only p < 0.01 gives X_m > 0
    p = c(0.001, 0.01, 0.1, 0.5)
    xi = 2
    result = local_score(p, xi = xi)
    expect_equal(result$score, -log10(p) - xi, tolerance = 1e-10)
    # Only p < 10^(-xi) = 0.01 gives positive score
    expect_true(result$score[1] > 0)   # p = 0.001 < 0.01
    expect_equal(result$score[2], 0, tolerance = 1e-10)  # p = 0.01 exactly -> score = 0
    expect_true(result$score[3] < 0)   # p = 0.1 > 0.01
})

test_that("Score transformation with xi=1 shifts cut-off to p < 0.1", {
    p = c(0.001, 0.01, 0.1, 0.5)
    result = local_score(p, xi = 1)
    # With xi=1, cut-off is p < 0.1 = 10^(-1)
    expect_true(result$score[1] > 0)   # p = 0.001 < 0.1
    expect_true(result$score[2] > 0)   # p = 0.01 < 0.1
    expect_equal(result$score[3], 0, tolerance = 1e-10)  # p = 0.1 exactly
    expect_true(result$score[4] < 0)   # p = 0.5 > 0.1
})


###################
# Autocorrelation #
###################


test_that("autocor returns 1.0 for a perfectly correlated sequence", {
    # monotonic sequences are perfectly auto-correlated (cor(1:4, 2:5) == 1)
    expect_equal(gwashelpers:::autocor(1:5), 1.0, tolerance = 1e-10)
    expect_equal(gwashelpers:::autocor(as.numeric(1:20)), 1.0, tolerance = 1e-10)
})

test_that("autocor returns 1.0 for a perfectly anti-correlated sequence", {
    expect_equal(gwashelpers:::autocor(c(1, -1, 1, -1, 1)), 1.0, tolerance = 1e-10)
})


test_that("autocor recovers phi from an AR(1) process", {
    # Autoregressive model [a.k.a. AR(1)], i.e. X[t] = phi * X[t-1] + eps. AR(1)
    # with phi=0.3 on raw values (not pnorm-transformed) should give  autocor
    # approximately equal to phi=0.3

    set.seed(42)
    n = 500
    x = numeric(n)
    x[1] = rnorm(1)
    for (i in 2:n) x[i] = 0.3 * x[i - 1] + rnorm(1)
    rho = gwashelpers:::autocor(x)
    expect_true(rho > 0.25 && rho < 0.35,
        info = paste("Expected rho approx 0.3, got", rho))
})

test_that("Lindley_thresh matches analytical value (xi=2)", {
    # From Sup Info file 1 in Fariello et al 2017 bottom of page 6:

    # For xi =1
    # a − log(M) = −5.50ρ^3 + 6.76ρ^2 - 5.66ρ − 2.51 (S3)
    # b = −1.22ρ^2 + 3.17ρ − 1.99, (S4)

    # For xi = 2
    # a − log(M) = 2.47ρ3 − 4.16ρ2 − 1.82ρ − 4.58 (S5)
    # b = 0.37ρ2 + 2.14ρ − 2.35 (S6)

    # At rho=0, xi=2, n=10000, alpha=0.05:
    # Uniform random p-values are nearly uncorrelated (autocor approx 0)
    set.seed(999)
    p = runif(10000)
    rho = gwashelpers:::autocor(p)
    expect_true(abs(rho) < 0.05,
        info = paste("autocor =", rho, "-- expected approx 0"))

    # Hand-computed threshold with rho=0, n=10000, xi=2, alpha=0.05

    a0 = log(10000) + (-4.58)  # polynomial terms vanish at rho=0
    b0 = -2.35
    expected0 = (log(-log(1 - 0.05)) - a0) / b0  # approx 3.234


    actual = lindley_thresh(p, xi = 2, alpha = 0.05)
    expect_equal(actual, expected0, tolerance = 0.1)


    # AR(1) directly in [0,1]: arp[i] = phi*arp[i-1] + (1-phi)*runif(1).
    # Each step is a convex combination, so values stay in [0,1].
    # Stationary Pearson lag-1 correlation equals phi exactly:
    #   Cov[arp[i], arp[i-1]] = phi * Var[arp]  =>  Corr = phi = 0.3
    set.seed(321)
    phi = 0.3
    arp = numeric(10000)
    arp[1] = runif(1)
    for (i in 2:10000) arp[i] = phi * arp[i - 1] + (1 - phi) * runif(1)
    rho = gwashelpers:::autocor(arp)
    expect_equal(rho, 0.3, tolerance = 0.05)

    a0 = log(10000) + 2.47*rho^3 - 4.16 * rho^2 - 1.82 * rho - 4.58
    b0 = 0.37 * rho^2 + 2.14 * rho - 2.35
    expected = (log(-log(1 - 0.05)) - a0) / b0
    actual = lindley_thresh(arp, xi = 2, alpha = 0.05)
    expect_equal(actual, expected, tolerance = 0.1)

})

test_that("Lindley threshold increases with number of markers", {
    # Use uncorrelated p-values so autocor is near 0 for both
    set.seed(123)
    p_small = runif(1000)
    p_large = runif(100000)

    thres_small = lindley_thresh(p_small, xi = 2, alpha = 0.05)
    thres_large = lindley_thresh(p_large, xi = 2, alpha = 0.05)

    expect_true(thres_large > thres_small,
        info = paste("thres_small =", thres_small, "; thres_large =", thres_large))
})


test_that("lindley_thresh returns NULL for invalid xi", {
    set.seed(1)
    p = runif(100)
    expect_null(suppressMessages(lindley_thresh(p, xi = 5)))
    expect_null(suppressMessages(lindley_thresh(p, xi = 0)))
})


####################
# Helper functions #
####################

test_that("windows() labels consecutive TRUE runs correctly", {
    state = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
    expected = c(NA_real_, 1, 1, NA_real_, 2, NA_real_, 3, 3)
    result = windows(state)
    expect_equal(result, expected)
})

test_that("windows() returns all NA for all-FALSE input", {
    expect_equal(windows(c(FALSE, FALSE, FALSE)), c(NA_real_, NA_real_, NA_real_))
})

test_that("windows() labels a single TRUE run as window 1", {
    expect_equal(windows(c(TRUE, TRUE, TRUE)), c(1, 1, 1))
})


#####################
# Integration Tests #
#####################
test_that("local_score() returns a data.frame with correct structure", {
    set.seed(42)
    p = runif(50, 0.001, 1)
    result = local_score(p, xi = 2)

    expect_s3_class(result, "data.frame")
    expect_named(result, c("score", "lindley_score", "lindley_thresh", "lindley_window"))
    expect_equal(nrow(result), length(p))
    # Threshold is a scalar replicated to every row
    expect_true(length(unique(result$lindley_thresh)) == 1)
    # Lindley scores are non-negative
    expect_true(all(result$lindley_score >= 0))
    # Window labels are positive integers or NA
    expect_true(all(is.na(result$lindley_window) | result$lindley_window >= 1))
})

test_that("local_score() rejects out-of-range p-values", {
    expect_error(local_score(c(0.5, 1.5), xi = 2))
    expect_error(local_score(c(-0.1, 0.5), xi = 2))
})

test_that("local_score() detects a clear signal cluster (integration test)", {
    # p-values: 10 background markers (p=0.3), 5 very significant (p=0.0001),
    # 10 background markers. The Lindley score should peak in the signal cluster
    # and that cluster should be labeled as window 1 in lindley_window.
    p = c(rep(0.3, 10), rep(0.0001, 5), rep(0.3, 10))
    result = local_score(p, xi = 1)

    # Peak of Lindley process should be in the middle cluster (positions 11-15)
    peak_pos = which.max(result$lindley_score)
    expect_true(peak_pos >= 11 && peak_pos <= 15,
        info = paste("Peak at position", peak_pos))

    # Lindley score at peak should exceed scores at the flanks
    expect_true(result$lindley_score[peak_pos] > result$lindley_score[1])
    expect_true(result$lindley_score[peak_pos] > result$lindley_score[25])

    # If the signal is strong enough to exceed the threshold, the window column
    # should label the peak cluster as a distinct window (all non-NA values equal)
    above = result$lindley_window[!is.na(result$lindley_window)]
    if (length(above) > 0) {
        expect_true(all(above == above[1]))
    }
})
