#' Compute the Lindley process from scores
#'
#' @param scores Log-scaled and adjusted p value scores
lindley = function(scores) {
    L = length(scores)
    sl = rep(0, L + 1)
    for (i in 1:L) {
        sl[i + 1] = max(0, (sl[i] + scores[i]))
    }
    return(sl[-1])
}
    
autocor = function(x) {
    abs(cor(x[-1], x[-length(x)]))
}

# Compute significance threshold if the distribution of p-values is uniform
lindley_thresh = function(pval, xi=2, alpha = 0.05) {
    # Coefficients for polynomial approximation
    a.poly.coef = rbind(
        c(-5.5, 2.47, 2.04, 0.22),    # coeff associated with cor^3 for xi=1 to xi=4
        c(6.76, -4.16, -5.76, -4.08), # coeff associated with cor^2 for xi=1 to xi=4
        c(-5.66, -1.82, 1.04, 1.16),  # coeff associated with cor for xi=1 to xi=4
        c(-2.51, -4.58, -6.95, -9.16) # intercept term for xi=1 to xi=4
    )
    b.poly.coef = rbind(
        c(-1.22, 0.37, 2.55, 3.45),   # coeff associated with cor^2 for xi=1 to xi=4
        c(3.17, 2.14, -0.02, -0.98),  # coeff associated with cor for xi=1 to xi=4
        c(-1.99, -2.35, -2.31, -2.33) # intercept term for xi=1 to xi=4
    )
    
    if (!xi %in% 1:4) {
        print('xi must be equal to 1, 2, 3 or 4')
        thres = NULL
    } else {
        cor = autocor(pval)
        nmarker = length(pval)
        a = log(nmarker) + a.poly.coef[1, xi] * (cor^3) + a.poly.coef[2, xi] * (cor^2) +
            a.poly.coef[3, xi] * cor + a.poly.coef[4, xi]
        b = b.poly.coef[1, xi] * (cor^2) + b.poly.coef[2, xi] * cor + b.poly.coef[3, xi]
        thres = (log(-log(1 - alpha)) - a) / b
    }
    return(thres)
}


lindley_p = function(p, xi=2) {
    if (any(p>1 | p < 0)) {
        stop("p values should be provided as-is, i.e. values between 0-1")
    }

    logp = -log10(p)
    score = logp - xi
    lscore = lindley(score)
    return (data.frame(score=score, lindley_score=lscore, lindley_thresh=lindley_thresh(p, xi)))
}

windows = function(state) {
    res = numeric(length(state))
    n = 0
    last = F
    for (i in seq_along(state)) {
        if (state[i]) {
            if (!last) n = n + 1
            res[i] = n
            last = T
        } else {
            res[i] = NA
            last = F
        }
    }
    res
}

