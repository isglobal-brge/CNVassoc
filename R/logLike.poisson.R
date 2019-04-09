logLike.poisson <- function(param, y, X, 
        w, variant) {
        J <- NCOL(w)
        beta <- vector2matrix(betav = param, 
            variant = variant, J = J)
        eta <- vapply(seq_len(J), function(j) {
            Xj <- cbind(X[, , j])
            betaj <- beta[, j, drop = FALSE]
            Xj %*% betaj
        }, rep(0, length(y)))
        lambda <- exp(eta)
        ff <- vapply(seq(1,J), function(j) dpois(y, lambda[, j]), 
                rep(0, length(y)))
        gg <- rowSums(w * ff)
        return(sum(log(gg)))
}
