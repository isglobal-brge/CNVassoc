EMWeibull <- function(y, cens, X, w, beta, 
        alpha, tol = 10^-6, max.iter = 1000, 
        verbose = FALSE) {
        variant <- !apply(is.na(beta), 1, any)
        beta <- t(apply(beta, 1, function(x) if (any(is.na(x))) 
            rep(x[1], length(x)) else x))
        P <- ncol(X)
        J <- ncol(w)
        N <- length(y)
        error <- Inf
        iter <- 1
        latent <- rep(seq_len(J), each = N)
        X.expand <- NULL
        for (p in seq_len(P)) {
            xp <- cbind(as.vector(X[, p, ]))
            colnames(xp) <- paste("V", p, sep = "")
            if (variant[p]) {
                xp <- model.matrix(~xp:as.factor(latent) - 
                    1)
                colnames(xp) <- paste(rep("V", 
                    J), p, sep = "")
            }
            X.expand <- cbind(X.expand, xp)
        }
        theta.old <- c(as.vector(beta), alpha)
        y.expand <- rep(y, J)
        cens.expand <- rep(cens, J)
        while (error > tol & iter < max.iter) {
            beta.old <- beta
            eta <- vapply(seq(1,J), function(j) {
                Xj <- cbind(X[, , j])
                betaj <- beta[, j, drop = FALSE]
                Xj %*% betaj
            }, rep(0, N))
            lambda <- exp(eta)  ##
            scale <- lambda^(-1/alpha)  ##
            dY <- matrix(0, nrow = N, ncol = J)
            ## not censored
            dY[cens == 1, ] <- vapply(seq(1,J), function(j) dweibull(y[cens == 
                1], alpha, scale[cens == 1, j]), rep(0, sum(cens==1)))
            ## censored
            dY[cens == 0, ] <- vapply(seq(1,J), function(j) pweibull(y[cens == 
                0], alpha, scale[cens == 0, j], 
                lower.tail = FALSE), rep(0, sum(cens==0)))
            wp <- dY * w
            wp <- wp/rowSums(wp)
            coeffic <- weibull.fit(x = X.expand, 
                y = y.expand, cens = cens.expand, 
                weights = as.vector(wp))  ##
            beta <- vector2matrix(betav = coeffic$beta, 
                variant = variant, J = J)
            alpha <- coeffic$alpha
            theta <- c(as.vector(beta), alpha)
            error <- max(abs(theta - theta.old))
            theta.old <- theta
            colnames(beta) <- paste("clust", seq_len(J), sep = "")
            rownames(beta) <- dimnames(X)[2][[1]]
            if (verbose) {
                if (iter == 1) 
                    cat("--- EM procedure ---\n")
                cat("Iter", iter, "Error", error, 
                    "\n")
                print(beta)
                cat("alpha=", alpha, "\n")
                cat("\n")
            }
            iter <- iter + 1
        }
        return(list(beta = beta, alpha = alpha, 
            variant = variant))
}
