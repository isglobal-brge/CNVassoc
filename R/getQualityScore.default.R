getQualityScore.default <- function(x, sds, 
        w, type, iter = 10000, threshold = 0.1, 
        ...) {
        mu <- x
        J <- length(mu)
        qs.type <- charmatch(tolower(type), tolower(c("class", 
            "CNVtools", "CANARY")))
        if (is.na(qs.type)) 
            stop(" argument 'type' must be either 'class' , 'CNVtools' or 
                'CANARY'")
        if (qs.type == 1) {
            p <- c()
            for (j in seq_len(J)) {
                X <- rnorm(iter, mu[j], sds[j])
                Y <- vapply(seq(1,J), function(s) w[s] * 
                    dnorm(X, mu[s], sds[s]), rep(0, iter))
                p <- c(p, mean(apply(Y, 1, which.max) == j))
            }
            out <- sum(p * w)
        }
        if (qs.type == 2) {
            sds <- sds[order(mu)]
            w <- w[order(mu)]
            mu <- sort(mu)
            dmu <- abs(mu[seq_len(J - 1)] - mu[seq(2,J)])
            av.sds <- (w[seq_len(J - 1)] * sds[seq_len(J - 1)] + 
                w[seq(2,J)] * sds[seq(2,J)])/(w[seq_len(J - 1)] + w[seq(2,J)])
            weights <- w[seq_len(J - 1)] * w[seq(2,J)]
            out <- sum(weights * dmu/av.sds)/sum(weights)
        }
        if (qs.type == 3) {
            f <- function(x) {
                Y <- vapply(seq(1,J), function(s) w[s] * 
                    dnorm(x, mu[s], sds[s]), rep(0, length(x)))
                index <- which.max(Y)
                max1 <- Y[index]
                max2 <- max(Y[-index])
                ratio <- max2/max1
                sum(Y) * as.integer(ratio > threshold)
            }
            minim <- which.min(mu)
            maxim <- which.max(mu)
            limits <- c(mu[minim] - 3 * sds[minim], 
                mu[maxim] + 3 * sds[maxim])
            fapply <- function(x) vapply(x, f, 0)
            out <- integrate(fapply, limits[1], 
                limits[2])$value
            attr(out, "threshold") <- threshold
        }
        attr(out, "type") <- qs.type
        return(out)
}
