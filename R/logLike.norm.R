logLike.norm <- function(param, y, X, w, 
        variant) {
        
        J <- NCOL(w)
        K <- sum(ifelse(variant, J, 1))
        
        beta <- vector2matrix(betav = param[seq_len(K)], 
            variant = variant, J = J)
        sigma <- param[K + 1]
        
        eta <- vapply(seq(1,J), function(j) {
            Xj <- cbind(X[, , j])
            betaj <- beta[, j, drop = FALSE]
            Xj %*% betaj
        }, rep(0, length(y)))
        
        mu <- eta
        ff <- vapply(seq(1,J), function(j) dnorm(y, mu[, j], sigma), 
                    rep(0, length(y)))
        gg <- rowSums(w * ff)
        return(sum(log(gg)))
        
}
