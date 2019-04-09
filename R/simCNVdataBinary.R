simCNVdataBinary <- function(n, mu.surrog, 
        sd.surrog, w, p0, or, cnv.random = FALSE) {
        k <- length(w)
        if (cnv.random) {
            nj <- rmultinom(n = 1, size = n, 
                prob = w)
        } else {
            nj <- round(n * w)
            nj[1] <- nj[1] + n - sum(nj)
        }
        surrog <- unlist(lapply(seq(along = nj), 
            function(j) rnorm(nj[j], mean = mu.surrog[j], 
                sd = sd.surrog[j])))
        cnv <- rep(seq_len(k), nj)
        p <- rep(NA, k)
        p[1] <- p0
        for (j in 2:k) p[j] <- p[1] * or[j - 
            1]/(1 + p[1] * (or[j - 1] - 1))
        resp <- unlist(lapply(seq(1,k), function(j) rbinom(nj[j], 1, p[j])))
        out <- data.frame(resp, cnv, surrog)
        out <- out[sample(seq_len(n)), ]
        out
}
