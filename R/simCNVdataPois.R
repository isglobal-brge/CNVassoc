simCNVdataPois <- function(n, mu.surrog, 
        sd.surrog, w, lambda, cnv.random = FALSE) {
        w <- w/sum(w)
        k <- length(mu.surrog)
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
        resp <- unlist(lapply(seq(1,k), function(j) rpois(nj[j], lambda[j])))
        out <- data.frame(resp, cnv, surrog)
        out <- out[sample(seq_len(n)), ]
        out
}
