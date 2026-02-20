#' Local Score / Lindley Process for GWAS
#'
#' @title Local Score / Lindley Process for GWAS
#'
#' @details
#' ## Overview
#'
#' The **local score** method aggregates weak association signals from nearby,
#' correlated markers along the genome into a single scalar rather than treating
#' each SNP independently (which ignores LD structure). The **Lindley process**
#' (a running maximum of a random walk) accumulates signal as it walks along the
#' chromosome, resetting to zero whenever the cumulative evidence becomes
#' negative. A genomic region is declared significant when this running total
#' exceeds a threshold derived from extreme-value (Gumbel) theory.
#'
#' The primary entry-point is the function local_score(), which takes a vector
#' of p-values and generates the local scores as a data frame.
#'
#' ## Formulae (Fariello et al. 2017, Equations 1-4)
#'
#' **Score per marker** (Eq. 1):
#' \deqn{X_m = -\log_{10}(p_m) - \xi}
#' where \eqn{\xi \in \{1,2,3,4\}} controls the p-value cut-off: only markers
#' with \eqn{p_m < 10^{-\xi}} contribute positively to the process.
#'
#' **Lindley recurrence** (Eq. 2):
#' \deqn{h_0 = 0, \quad h_m = \max(0,\; h_{m-1} + X_m)}
#'
#' **Local score** (Eq. 2):
#' \deqn{H_M = \max_{1 \le m \le M} h_m}
#'
#' **Null CDF - Gumbel approximation** (Eq. 3):
#' \deqn{P(H_M \le x) \approx 1 - \exp\!\bigl(-\exp(a_{M,\rho} + b_{M,\rho}\, x)\bigr)}
#'
#' **Significance threshold** (Eq. 4):
#' \deqn{t_{M,\rho,\alpha} \approx \frac{\log(-\log(1-\alpha)) - a_{M,\rho}}{b_{M,\rho}}}
#'
#' Note: the paper writes \eqn{\log(-\log(\alpha))} where \eqn{\alpha} is
#' \eqn{1 - \text{FWER}}; the R code uses the equivalent survival-function form
#' \eqn{\log(-\log(1-\alpha))} where \eqn{\alpha} is the FWER directly.
#'
#' ## Gumbel Parameter Approximation
#'
#' The Gumbel location and scale parameters depend on the number of markers
#' \eqn{M} and the lag-1 autocorrelation \eqn{\rho} (a proxy for LD):
#' \deqn{a_{M,\rho} = \log(M) + c_3\rho^3 + c_2\rho^2 + c_1\rho + c_0}
#' \deqn{b_{M,\rho} = d_2\rho^2 + d_1\rho + d_0}
#'
#' The polynomial coefficients \eqn{c_i, d_i} were estimated by simulation over
#' \eqn{M \in [100, 45000]} and \eqn{\rho \in [0, 0.9]} and are tabulated for
#' \eqn{\xi \in \{1,2,3,4\}} in `lindley_thresh()`.
#'
#' ## References
#'
#' Fariello, M.-I., Boitard, S., Mercier, S., Robelin, D., Faraut, T.,
#' Liaubet, L., ... & Servin, B. (2017). Accounting for linkage disequilibrium
#' in genome scans for selection without panels of controls.
#' *Molecular Ecology*, 26(22), 5994-6009.
#' \doi{10.1111/mec.14141}
#'
#' Bonhomme, M., et al. (2019). Genomic signature of selective sweeps illuminates
#' adaptation of *Medicago truncatula*. *Molecular Ecology*, 28(21), 4751-4771.
#' <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6781128/>
#'
#' @name localscore-module
NULL


#' Compute the Lindley process from a vector of scores
#'
#' Walks along a vector of per-marker scores and computes the Lindley (running
#' maximum) process, which forms the basis of the local score statistic.
#'
#' @param scores Numeric vector of per-marker scores
#'   \eqn{X_m = -\log_{10}(p_m) - \xi}. Positive values represent markers with
#'   evidence stronger than the cut-off \eqn{10^{-\xi}}; negative values
#'   accumulate evidence against association.
#'
#' @return Numeric vector of the same length as `scores`, giving the Lindley
#'   process values \eqn{h_1, h_2, \ldots, h_M} (all \eqn{\ge 0}). The
#'   **local score** is `max(lindley(scores))`.
#'
#' @details
#' Implements the Mercier & Daudin (2001) recurrence:
#' \deqn{h_0 = 0, \quad h_m = \max(0,\; h_{m-1} + X_m).}
#' The process is non-negative by construction: whenever cumulative evidence
#' drops below zero, the process resets, beginning a new potential excursion.
#' Consecutive positions where \eqn{h_m > t} (threshold \eqn{t}) constitute an
#' *excursion interval* corresponding to a genomic association peak.
#'
#' @export
#'
#' @examples
#' scores = c(0.5, -0.3, 0.8, -1.0, 0.2, 0.3)
#' lindley(scores)
#' # [1] 0.5 0.2 1.0 0.0 0.2 0.5
#' max(lindley(scores))  # local score = 1.0
lindley = function(scores) {
    L = length(scores)
    sl = rep(0, L + 1)
    for (i in 1:L) {
        sl[i + 1] = max(0, (sl[i] + scores[i]))
    }
    return(sl[-1])
}


#' Compute the absolute lag-1 autocorrelation of a numeric vector
#'
#' @param x Numeric vector (typically raw p-values from adjacent markers).
#'
#' @return A single non-negative number: the absolute value of the lag-1
#'   Pearson correlation between `x[-1]` and `x[-length(x)]`. Used as the
#'   LD proxy \eqn{\rho} in `lindley_thresh()`.
#'
#' @details
#' Adjacent markers in LD tend to have correlated association p-values. This
#' function estimates that correlation as the absolute lag-1 Pearson correlation
#' of the raw p-value vector - a simple but effective surrogate for \eqn{\rho}.
#' The absolute value is taken because both positive and negative lag-1
#' correlations inflate the local score null distribution equally.
#'
#' @keywords internal
autocor = function(x) {
    abs(stats::cor(x[-1], x[-length(x)]))
}


#' Compute the local score significance threshold
#'
#' Returns the genome-wide significance threshold for the local score (Lindley
#' process maximum) using the Gumbel extreme-value approximation of
#' Fariello et al. (2017), Equation 4.
#'
#' @param pval Numeric vector of raw p-values (values in \eqn{[0, 1]}).
#'   Used to compute both the number of markers \eqn{M = \text{length}(pval)}
#'   and the lag-1 autocorrelation \eqn{\rho} via `autocor()`.
#' @param xi Integer tuning parameter; must be 1, 2, 3, or 4. Controls the
#'   p-value cut-off: only markers with \eqn{p < 10^{-\xi}} contribute
#'   positively to the Lindley process. Must match the `xi` used when computing
#'   scores with `local_score()`.
#' @param alpha Numeric family-wise type-I error rate (default 0.05). The
#'   threshold \eqn{t} satisfies \eqn{P(H_M > t) \le \alpha} under the null.
#'
#' @return Numeric threshold value, or `NULL` (with a message) if `xi` is not
#'   in \eqn{\{1, 2, 3, 4\}}.
#'
#' @details
#' The Gumbel location and scale parameters are approximated as polynomials in
#' \eqn{\rho} (Fariello et al. 2017, Supplementary Information):
#'
#' \deqn{a_{M,\rho} = \log(M) + c_3\rho^3 + c_2\rho^2 + c_1\rho + c_0}
#' \deqn{b_{M,\rho} = d_2\rho^2 + d_1\rho + d_0}
#'
#' Polynomial coefficients (rows = polynomial term, columns = \eqn{\xi}):
#'
#' | Term | \eqn{\xi=1} | \eqn{\xi=2} | \eqn{\xi=3} | \eqn{\xi=4} |
#' |------|------------|------------|------------|------------|
#' | \eqn{a}: \eqn{\rho^3} | -5.50 | 2.47 | 2.04 | 0.22 |
#' | \eqn{a}: \eqn{\rho^2} | 6.76 | -4.16 | -5.76 | -4.08 |
#' | \eqn{a}: \eqn{\rho} | -5.66 | -1.82 | 1.04 | 1.16 |
#' | \eqn{a}: intercept | -2.51 | -4.58 | -6.95 | -9.16 |
#' | \eqn{b}: \eqn{\rho^2} | -1.22 | 0.37 | 2.55 | 3.45 |
#' | \eqn{b}: \eqn{\rho} | 3.17 | 2.14 | -0.02 | -0.98 |
#' | \eqn{b}: intercept | -1.99 | -2.35 | -2.31 | -2.33 |
#'
#' The threshold is (Eq. 4):
#' \deqn{t = \frac{\log(-\log(1-\alpha)) - a}{b}}
#'
#' Note on parameterisation: the paper writes \eqn{\log(-\log(\alpha))} where
#' their \eqn{\alpha} denotes \eqn{1 - \text{FWER}}. This implementation uses
#' \eqn{\log(-\log(1-\alpha))} where \eqn{\alpha} is the FWER directly, which
#' is the conventional notation.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' p = runif(1000)
#' lindley_thresh(p, xi = 2, alpha = 0.05)
lindley_thresh = function(pval, xi = 2, alpha = 0.05) {
    # Polynomial coefficients for Gumbel parameter approximation.
    # Rows: polynomial degree (highest first), Columns: xi = 1, 2, 3, 4.
    # Fitted by simulation over M in [100, 45000] and rho in [0, 0.9].
    a.poly.coef = rbind(
        c(-5.5, 2.47, 2.04, 0.22),    # coefficient for rho^3
        c(6.76, -4.16, -5.76, -4.08), # coefficient for rho^2
        c(-5.66, -1.82, 1.04, 1.16),  # coefficient for rho
        c(-2.51, -4.58, -6.95, -9.16) # intercept
    )
    b.poly.coef = rbind(
        c(-1.22, 0.37, 2.55, 3.45),   # coefficient for rho^2
        c(3.17, 2.14, -0.02, -0.98),  # coefficient for rho
        c(-1.99, -2.35, -2.31, -2.33) # intercept
    )

    if (!xi %in% 1:4) {
        message("xi must be equal to 1, 2, 3 or 4")
        thres = NULL
    } else {
        cor = autocor(pval)
        nmarker = length(pval)
        a = log(nmarker) +
            a.poly.coef[1, xi] * (cor^3) +
            a.poly.coef[2, xi] * (cor^2) +
            a.poly.coef[3, xi] * cor +
            a.poly.coef[4, xi]
        b = b.poly.coef[1, xi] * (cor^2) +
            b.poly.coef[2, xi] * cor +
            b.poly.coef[3, xi]
        thres = (log(-log(1 - alpha)) - a) / b
    }
    return(thres)
}


#' Label consecutive runs of TRUE as distinct genomic windows
#'
#' Given a logical vector (typically `lindley_score > lindley_thresh`), assigns
#' consecutive runs of `TRUE` values a unique integer label. Useful for
#' enumerating distinct association peaks detected by the local score.
#'
#' @param state Logical vector. Typically the result of
#'   `result$lindley_score > result$lindley_thresh` from `local_score()`.
#'
#' @return Integer vector of the same length as `state`. Each maximal run of
#'   consecutive `TRUE` values is labeled with a unique positive integer
#'   (1, 2, 3, ...). Positions where `state` is `FALSE` are set to `NA`.
#'
#' @details
#' Each contiguous block of `TRUE` positions corresponds to a genomic excursion
#' interval - the span of markers within a single association peak. The integer
#' labels allow downstream grouping (e.g., `split()` or `dplyr::group_by()`)
#' to extract and summarize each peak separately.
#'
#' @export
#'
#' @examples
#' state = c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
#' find_windows(state)
#' # [1] NA  1  1 NA  2 NA  3  3
find_windows = function(state) {
    res = numeric(length(state))
    n = 0
    last = FALSE
    for (i in seq_along(state)) {
        if (state[i]) {
            if (!last) n = n + 1
            res[i] = n
            last = TRUE
        } else {
            res[i] = NA
            last = FALSE
        }
    }
    res
}


#' Apply the local score method to a vector of p-values
#'
#' The main entry point for local score analysis. Transforms raw p-values into
#' per-marker scores, computes the Lindley process, and attaches the
#' genome-wide significance threshold.
#'
#' @param p Numeric vector of raw p-values (must be in \eqn{[0, 1]}).
#' @param xi Integer tuning parameter (1, 2, 3, or 4). Controls the p-value
#'   threshold for positive contribution: only \eqn{p < 10^{-\xi}} markers
#'   receive a positive score \eqn{X_m = -\log_{10}(p) - \xi > 0}.
#'
#' @return A `data.frame` with one row per marker and four columns:
#'   \describe{
#'     \item{`score`}{Per-marker score \eqn{X_m = -\log_{10}(p_m) - \xi}.}
#'     \item{`lindley_score`}{Lindley process value \eqn{h_m \ge 0}.}
#'     \item{`lindley_thresh`}{Genome-wide threshold (same value in every row),
#'       computed by `lindley_thresh()`.}
#'     \item{`lindley_window`}{Integer window label from `find_windows()`: consecutive
#'       markers above the threshold share the same label (1, 2, 3, ...);
#'       markers below the threshold are `NA`.}
#'   }
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' p = c(runif(20, 0.1, 1), runif(5, 1e-5, 1e-4), runif(20, 0.1, 1))
#' result = local_score(p, xi = 2)
#' head(result)
local_score = function(p, xi = 2) {
    if (any(p > 1 | p < 0)) {
        stop("p values should be provided as-is, i.e. values between 0-1")
    }
    score = -log10(p) - xi
    lscore = lindley(score)
    thresh = lindley_thresh(p, xi)
    return(data.frame(
        score = score,
        lindley_score = lscore,
        lindley_thresh = thresh,
        lindley_window = find_windows(lscore > thresh)
    ))
}

