fastCNVinter <- function(probs, formula, 
        data, model = "additive", family = "binomial", 
        nclass = 3, colskip = 5, tol = 1e-06, 
        max.iter = 30, verbose = FALSE, multicores = 0) {
        cl <- match.call()
        if (missing(data)) 
            data <- environment(formula)
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data", "na.action"), 
            names(mf), 0L)  #offset
        mf <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        mt <- attr(mf, "terms")
        model.type <- charmatch(model, c("additive", 
            "multiplicative"))
        if (is.na(model.type)) 
            stop(" argument 'model' must be either 'multiplicative' or 
            'additive'")
        if (model.type != 1) 
            stop(" only 'additive' model is implemented")
        if (!family %in% c("weibull", "binomial")) 
            stop(" only 'biomial' and 'weibull' models are implemented by now")
        nIndiv <- NROW(mf)
        mm <- model.matrix(mt, mf, contrasts)
        y <- model.response(mf)
        nVar <- NCOL(mf) - 1
        if (nVar != 0) {
            Xcov <- mm[, -1, drop = FALSE]
            nCov <- NCOL(Xcov)
        } else {
            Xcov <- NULL
            nCov <- 0
        }
        N <- as.integer(nIndiv)
        
        # RESPONSE
        if (family != "weibull") 
            Y <- y else {
            Y <- y[, 1]
            CENS <- y[, 2]
        }
        
        # PARAM INI
        if (nVar > 0) {
            if (family == "binomial") 
                fit0 <- glm(Y ~ Xcov, family = "binomial")
            if (family == "weibull") 
                fit0 <- survreg(Surv(Y, CENS) ~ 
                    Xcov)
        } else {
            if (family == "binomial") 
                fit0 <- glm(Y ~ 1, family = "binomial")
            if (family == "weibull") 
                fit0 <- survreg(Surv(Y, CENS) ~ 
                    1)
        }
        if (family != "weibull") 
            PARAM <- c(0, 0, 0, coef(fit0)) else PARAM <- c(0, 0, 0, 
                                                    unlist(translate(fit0)))
        PARAM[c(1, 4)] <- PARAM[c(4, 1)]
        
        # OTHER INPUTS
        TOL <- tol
        MAX_ITER <- as.integer(max.iter)
        VERBOSE <- as.integer(0)
        P <- as.integer(nclass)
        Q <- as.integer(nCov)
        Y <- as.double(Y)
        if (nCov > 0) 
            C <- t(Xcov) else C <- NULL
        
        # PROBS
        if (is.character(probs)) {
            if (length(grep("\\.fst$", probs)) > 
                0) {
                # fst format
                cat("Scanning .fst data...\n")
                read.st <- system.time(PROBS <- read_fst(probs))[3]
                cat("Done! Took ", read.st, "seconds\n\n")
            } else {
                # .probs format
                cat("Reading .probs data...\n")
                read.st <- system.time(PROBS <- read_table2(probs, 
                    col_names = FALSE))[3]
                cat("Done! Took ", read.st, "seconds\n\n")
            }
        } else {
            if (!(is.matrix(probs) || is.data.frame(probs))) 
                stop(" 'probs' must be a character (file) or a matrix or a 
                    data.frame")
            PROBS <- probs
        }
        if (colskip > 0) 
            PROBS <- PROBS[, -seq_len(colskip), drop = FALSE]
        PROBS <- as.matrix(PROBS)
        nCNVs <- NROW(PROBS)
        
        mm <- outer(seq_len(nCNVs), seq_len(nCNVs), paste, 
            sep = ",")
        mm <- mm[upper.tri(mm)]
        
        # RUN no use of parallel package
        if (multicores == 0) {
            ans <- lapply(mm, function(index) {
                i <- as.integer(strsplit(index, 
                    split = ",")[[1]][1])
                j <- as.integer(strsplit(index, 
                    split = ",")[[1]][2])
                if (verbose) 
                    cat("Pair ", index, "\n")
                PROBS.i <- PROBS[i, ]
                PROBS.j <- PROBS[j, ]
                WW.i <- matrix(as.double(PROBS.i), 
                    nrow = nclass)
                WW.j <- matrix(as.double(PROBS.j), 
                    nrow = nclass)
                WW <- .Call("outerC", as.integer(N), 
                    t(WW.i), t(WW.j), PACKAGE="CNVassoc")
                if (nCov > 0) {
                    if (family == "binomial") 
                    fit <- try(.Call("NRinterlogisticcov", 
                        N, PARAM, Y, WW, Q, C, 
                        TOL, MAX_ITER, VERBOSE, PACKAGE="CNVassoc"))
                    if (family == "weibull") 
                    fit <- try(.Call("NRinterweibullcov", 
                        N, PARAM, Y, CENS, WW, 
                        Q, C, TOL, MAX_ITER, 
                        VERBOSE, PACKAGE="CNVassoc"))
                } else {
                    if (family == "binomial") 
                    fit <- try(.Call("NRinterlogistic", 
                        N, PARAM, Y, WW, TOL, 
                        MAX_ITER, VERBOSE, PACKAGE="CNVassoc"))
                    if (family == "weibull") 
                    fit <- try(.Call("NRinterweibull", 
                        N, PARAM, Y, CENS, WW, 
                        TOL, MAX_ITER, VERBOSE, PACKAGE="CNVassoc"))
                }
                if (inherits(fit, "try-error")) 
                    return(rep(NA, 3))
                logor <- fit[[1]][4]
                se <- try(sqrt(diag(solve(-fit[[4]])))[4], 
                    silent = TRUE)
                if (inherits(se, "try-error")) 
                    se <- NA
                iter <- fit[[5]]
                return(c(logor, se, iter))
            })
            ans <- matrix(unlist(ans), ncol = 3, 
                byrow = TRUE)
        } else {
            # use of parallel package
            requireNamespace("parallel", quietly = TRUE)
            ans <- parallel::mclapply(mm, function(index) {
                i <- as.integer(strsplit(index, 
                    split = ",")[[1]][1])
                j <- as.integer(strsplit(index, 
                    split = ",")[[1]][2])
                if (verbose) 
                    cat("Pair ", index, "\n")
                PROBS.i <- PROBS[i, ]
                PROBS.j <- PROBS[j, ]
                WW.i <- matrix(as.double(PROBS.i), 
                    nrow = nclass)
                WW.j <- matrix(as.double(PROBS.j), 
                    nrow = nclass)
                WW <- .Call("outerC", as.integer(N), 
                    t(WW.i), t(WW.j), PACKAGE="CNVassoc")
                if (nCov > 0) {
                    if (family == "binomial") 
                    fit <- try(.Call("NRinterlogisticcov", 
                        N, PARAM, Y, WW, Q, C, 
                        TOL, MAX_ITER, VERBOSE, PACKAGE="CNVassoc"))
                    if (family == "weibull") 
                    fit <- try(.Call("NRinterweibullcov", 
                        N, PARAM, Y, CENS, WW, 
                        Q, C, TOL, MAX_ITER, 
                        VERBOSE, PACKAGE="CNVassoc"))
                } else {
                    if (family == "binomial") 
                    fit <- try(.Call("NRinterlogistic", 
                        N, PARAM, Y, WW, TOL, 
                        MAX_ITER, VERBOSE, PACKAGE="CNVassoc"))
                    if (family == "weibull") 
                    fit <- try(.Call("NRinterweibull", 
                        N, PARAM, Y, CENS, WW, 
                        TOL, MAX_ITER, VERBOSE, PACKAGE="CNVassoc"))
                }
                if (inherits(fit, "try-error")) 
                    return(rep(NA, 3))
                logor <- fit[[1]][4]
                se <- try(sqrt(diag(solve(-fit[[4]])))[4], 
                    silent = TRUE)
                if (inherits(se, "try-error")) 
                    se <- NA
                iter <- fit[[5]]
                return(c(logor, se, iter))
            }, mc.cores = multicores)
            ans <- matrix(unlist(ans), ncol = 3, 
                byrow = TRUE)
        }
        
        zscore <- ans[, 1]/ans[, 2]
        pval <- 2 * (1 - pnorm(abs(zscore)))
        if (nrow(ans) == 1) 
            ans <- rbind(c(ans[1, c(1,2)], zscore, 
                pval, ans[1, 3])) else ans <- cbind(ans[, c(1,2)], zscore, 
            pval, ans[, 3])
        colnames(ans) <- c("beta", "se", "zscore", 
            "pvalue", "iter")
        ans <- data.frame(pair = mm, ans)
        ans
        
}
