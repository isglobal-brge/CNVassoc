EMnorm <- function(y, X, w, beta, sigma, 
        variant, tol = 10^-6, max.iter = 1000, 
        verbose = FALSE) {
        
        if (missing(variant)) 
            variant <- !apply(is.na(beta), 1, 
                any)
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
        
        y.expand <- rep(y, J)
        
        sigma.ext <- rep(sigma, J)
        
        while (error > tol & iter < max.iter) {
            
            beta.old <- beta
            sigma.old <- sigma
            theta.old <- c(matrix2vector(beta.old, 
                variant), sigma)
            
            eta <- vapply(seq(1,J), function(j) {
                Xj <- cbind(X[, , j])
                betaj <- beta[, j, drop = FALSE]
                Xj %*% betaj
            }, rep(0, N))
            
            dY <- vapply(seq(1,J), function(j) dnorm(y, 
                eta[, j], sigma.ext[j]), rep(0, N))
            wp <- dY * w
            rs <- rowSums(dY)
            if (any(rs == 0)) 
                dY[rs == 0, ] <- 1
            wp <- dY * w
            wp <- wp/rowSums(wp)
            
            fit <- linear.fit(x = X.expand, y = y.expand, 
                weights = as.vector(wp))
            
            beta <- vector2matrix(betav = fit$coeffic, 
                variant = variant, J = J)
            
            sigma <- fit$sigma
            
            theta <- c(matrix2vector(beta, variant), 
                sigma)
            
            error <- max(abs(theta - theta.old))
            
            colnames(beta) <- paste("clust", seq_len(J), sep = "")
            rownames(beta) <- dimnames(X)[2][[1]]
            
            if (verbose) {
                if (iter == 1) 
                    cat("---- EM procedure ----\n")
                cat("Iter", iter, "Error", error, 
                    "\n")
                cat("Beta:\n")
                print(beta)
                cat("Sigma:\n")
                print(sigma)
                cat("\n")
            }
            
            iter <- iter + 1
            
        }
        
        return(list(beta = beta, sigma = sigma, 
            variant = variant))
        
}
