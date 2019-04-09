logLike.weibull <- function(param, y, cens, 
        X, w, variant) {
        J <- NCOL(w)
        K <- sum(ifelse(variant, J, 1))
        beta <- vector2matrix(betav = param[seq_len(K)], 
            variant = variant, J = J)
        alpha <- param[K + 1]
        eta <- vapply(seq(1,J), function(j) {
            Xj <- cbind(X[, , j])
            betaj <- beta[, j, drop = FALSE]
            Xj %*% betaj
        }, rep(0, length(y)))
        lambda <- exp(eta)  ##
        scale <- lambda^(-1/alpha)  ##
        ff <- lambda
        ## not censored
        ff[cens == 1, ] <- vapply(seq_len(J), function(j) dweibull(y[cens == 
            1], alpha, scale[cens == 1, j]), rep(0, sum(cens==1)))
        ## censored
        ff[cens == 0, ] <- vapply(seq_len(J), function(j) pweibull(y[cens == 
            0], alpha, scale[cens == 0, j], lower.tail = FALSE), 
            rep(0, sum(cens==0)))
        gg <- rowSums(w * ff)
        return(sum(log(gg)))
}
