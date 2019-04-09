logLike.logistic <- function(param, y, X, 
        w, variant) {
        
        J <- NCOL(w)
        
        beta <- vector2matrix(betav = param, 
            variant = variant, J = J)
        
        eta <- vapply(seq(1,J), function(j) {
            Xj <- cbind(X[, , j])
            betaj <- beta[, j, drop = FALSE]
            Xj %*% betaj
        }, rep(0, length(y)))
        
        probs <- 1/(1 + exp(-eta))
        ff <- vapply(seq(1,J), function(j) dbinom(y, 1, probs[, j]), 
                    rep(0, length(y)))
        gg <- rowSums(w * ff)
        return(sum(log(gg)))
}
