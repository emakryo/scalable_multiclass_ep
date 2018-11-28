#########################################################################################################
#
# Implements a multi-class Gaussian Process classification using parallel EP updates for the likelihood.
#

t0 <- NULL
optim_params <- NULL


rowSums <- function (x, na.rm = FALSE, dims = 1L) {
    if (!is.array(x) || length(dn <- dim(x)) < 2L) 
        stop("'x' must be an array of at least two dimensions")
    if (dims < 1L || dims > length(dn) - 1L) 
        stop("invalid 'dims'")
    p <- prod(dn[-(1L:dims)])
    dn <- dn[1L:dims]
    z <- if (is.complex(x)) 
        .Internal(rowSums(Re(x), prod(dn), p, na.rm)) + (0+1i) * 
            .Internal(rowSums(Im(x), prod(dn), p, na.rm))
    else .Internal(rowSums(x, prod(dn), p, na.rm))
    if (length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[1L:dims]
    }
    else names(z) <- dimnames(x)[[1L]]
    z
}

colSums <- function (x, na.rm = FALSE, dims = 1L) {
    if (!is.array(x) || length(dn <- dim(x)) < 2L) 
        stop("'x' must be an array of at least two dimensions")
    if (dims < 1L || dims > length(dn) - 1L) 
        stop("invalid 'dims'")
    n <- prod(dn[1L:dims])
    dn <- dn[-(1L:dims)]
    z <- if (is.complex(x)) 
        .Internal(colSums(Re(x), n, prod(dn), na.rm)) + (0+1i) * 
            .Internal(colSums(Im(x), n, prod(dn), na.rm))
    else .Internal(colSums(x, n, prod(dn), na.rm))
    if (length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[-(1L:dims)]
    }
    else names(z) <- dimnames(x)[[dims + 1]]
    z
}

##
# Function that computes the product of a matrix and another matrix or vector
# given the cholesky factor of the inverse of the initial matrix.
#

inverse_times <- function(chol, x) {
    backsolve(t.default(chol), forwardsolve(chol, x))
}



###
# Function that computes the cavity distribution
#
# @param	a 		The approximation 
# @param q      The posterior factor
#
# @return	q_cavity 	The cavity distribution
#
#

cavityDistribution <- function(a, q){
    utEu_yi <- utEu_c <- matrix(0, a$n, a$nK)
    utm_yi <- utm_c <- matrix(0, a$n, a$nK)
    
    for (i in 1 : a$nK) {
        utEu_c[, i ] <- q[[i]]$utEu
        utm_c[, i ] <- q[[i]]$utm 
    }
    
    index <- cbind(1 : a$n, a$Y)
    utEu_yi <- matrix(utEu_c[ index ], a$n, a$nK)
    utm_yi <- matrix(utm_c[ index ], a$n, a$nK)
    
    utEu_yi_cav <- utEu_yi + utEu_yi * (a$f1Hat$C1_yi^-1 - utEu_yi)^-1 * utEu_yi
    utm_yi_cav  <- utm_yi + a$f1Hat$C1_yi * utEu_yi_cav * utm_yi - a$f1Hat$C2_yi * utEu_yi_cav
    
    utEu_c_cav <- utEu_c + utEu_c * (a$f1Hat$C1_c^-1 - utEu_c)^-1 * utEu_c
    utm_c_cav <- utm_c + a$f1Hat$C1_c * utEu_c_cav * utm_c - a$f1Hat$C2_c * utEu_c_cav
    
    
    q_cavity <- list(utEu_yi = utEu_yi_cav, utm_yi = utm_yi_cav, utEu_c = utEu_c_cav, utm_c = utm_c_cav,
                     utEpostu_yi = utEu_yi, utEpostu_c = utEu_c, utmpost_yi = utm_yi,
                     utmpost_c = utm_c)
    
    q_cavity
}


###
# Function that computes the cavity distribution
#
# @param	a 		The approximation 
# @param q      The posterior factor
#
# @return	q_cavity 	The cavity distribution
#
#

cavityDistribution_minibatch <- function(Y, nK, n, q, uold, u, C1_c, C1_yi, C2_c, C2_yi){
    
    utEu_yi <- utEu_c <- utEuold_yi <- utEuold_c <- uoldtEuold_yi <- uoldtEuold_c <- matrix(0, n, nK)
    uoldtm_yi <- uoldtm_c <- utm_yi <- utm_c <- matrix(0, n, nK)
    
    for (i in 1 : nK) {
        
        u_i <- matrix(u[ ,i, ], n, dim(u)[3])
        uold_i <- matrix(uold[ ,i, ], n, dim(uold)[3])

	    uE <- u_i %*% q[[i]]$E

        uoldtEuold_c[, i ] <- rowSums((uold_i %*% q[[i]]$E) * uold_i)
        utEuold_c[, i ] <- rowSums(uE * uold_i)
        utEu_c[, i ] <- rowSums(uE * u_i)
        uoldtm_c[, i ] <- uold_i%*% q[[i]]$m
        utm_c[, i ] <- u_i %*% q[[i]]$m 
        
        uold_yi <- matrix(uold[ which(Y == i), i, ], sum(Y == i), dim(uold)[3])
        u_yi <- matrix(u[ which(Y == i), i, ], sum(Y == i), dim(uold)[3])

	    Euoldyi <- q[[i]]$E %*% t(uold_yi)
	    
	    uoldtEuold_yi[ which(Y == i),  ] <- matrix(colSums(t(uold_yi) * Euoldyi), sum(Y == i), nK)
        utEuold_yi[ which(Y == i),  ] <- matrix(colSums(t(u_yi) * Euoldyi), sum(Y == i), nK)
        utEu_yi[ which(Y == i),  ] <- matrix(colSums(t(u_yi) * (q[[i]]$E %*% t(u_yi))), sum(Y == i), nK)
        
        u_c <- matrix(u[ which(Y == i),i, ], sum(Y == i), dim(u)[3])
        uold_c <- matrix(uold[ which(Y == i),i, ], sum(Y == i), dim(uold)[3])
        uoldtm_yi[ which(Y == i),  ] <- matrix(uold_c %*% q[[i]]$m, sum(Y == i), nK)
        utm_yi[ which(Y == i),  ] <- matrix(u_c %*% q[[i]]$m, sum(Y == i), nK)
    }
    
    utEu_yi_cav <- utEu_yi + utEuold_yi * (C1_yi ^ -1 - uoldtEuold_yi) ^ -1 * utEuold_yi
    utEuold_yi_cav <- utEuold_yi + utEuold_yi * (C1_yi ^ -1 - uoldtEuold_yi) ^ -1 * uoldtEuold_yi
    utm_yi_cav  <- utm_yi + C1_yi * utEuold_yi_cav * uoldtm_yi - C2_yi * utEuold_yi_cav
    
    utEu_c_cav <- utEu_c + utEuold_c * (C1_c^-1 - uoldtEuold_c) ^ -1 * utEuold_c
    utEuold_c_cav <- utEuold_c + utEuold_c * (C1_c^-1 - uoldtEuold_c) ^ -1 * uoldtEuold_c
    utm_c_cav <- utm_c + C1_c * utEuold_c_cav * uoldtm_c - C2_c * utEuold_c_cav
    
    q_cavity <- list(utEu_yi = utEu_yi_cav, utm_yi = utm_yi_cav, utEu_c = utEu_c_cav, utm_c = utm_c_cav,
                     utEpostu_yi = utEu_yi, utEpostu_c = utEu_c, utmpost_yi = utm_yi, utmpost_c = utm_c)
    
    q_cavity
}


checkCavity <- function(q_cavity, q, nK, Y, uold, u, C1_c_old, C1_yi_old, C2_c_old, C2_yi_old) {
    
    utVCavityu_yi <- utmCavity_yi <- utVCavityu_c <- utmCavity_c <- matrix(0, 5, nK)
    
    for (ex in 1:5){
        for (i in 1:nK) {
            VCavity_yi <- solve(solve(q[[Y[ex]]]$E) - (C1_yi_old[ex ,i] * uold[ex ,Y[ex], ]) %*% t(uold[ex ,Y[ex], ]))
            utVCavityu_yi[ex,i] <- u[ex,Y[ ex ],] %*% VCavity_yi %*% u[ex,Y[ ex ],]
            mCavity_yi <- VCavity_yi %*% (solve(q[[Y[ex]]]$E) %*% q[[Y[ex]]]$m - C2_yi_old[ex ,i] * uold[ex ,Y[ex], ])
            utmCavity_yi[ex,i] <- u[ex,Y[ ex ],] %*% mCavity_yi
            
            VCavity_c <- solve(solve(q[[i]]$E) - C1_c_old[ex, i] * uold[ex,i,] %*% t(uold[ex,i,]))
            utVCavityu_c[ex,i] <- u[ex,i,] %*% VCavity_c %*% u[ex,i,]
            mCavity_c <- VCavity_c %*% (solve(q[[i]]$E) %*% q[[i]]$m - C2_c_old[ex ,i] * uold[ex ,i, ])
            utmCavity_c[ex,i] <- u[ex,i,] %*% mCavity_c
            
        }
    }
    
    
    cbind(q_cavity$utEu_yi[1:5,], Y[1:5])
    cbind(utVCavityu_yi, Y[1:5])
    
    cbind(q_cavity$utEu_c[1:5,], Y[1:5])
    cbind(utVCavityu_c, Y[1:5])
    
    cbind(q_cavity$utm_yi[1:5,], Y[1:5])
    cbind(utmCavity_yi, Y[1:5])
    
    cbind(q_cavity$utm_c[1:5,], Y[1:5])
    cbind(utmCavity_c, Y[1:5])
}





checkPosterior <- function(q, a, nK, X, Y, f1Hat, Xbar, sigma, sigma0, l){
    
    # Reconstruct structure
    a$f1Hat <- f1Hat
    a$X <- X
    a$Y <- Y
    a$n <- nrow(X)
    for (i in 1:nK) {
        a$gMGPCinfo[[ i ]] <- initialize_kernel(X, Xbar[[ i ]], sigma[ i ], sigma0[ i ], l[ i, ])
    }
    
    
    # Reconstruct posterior
    qFull <- reconstructPosterior(a)

    # print(q[[1]]$E[1:5,1:5])
    # print(qFull[[1]]$E[1:5,1:5])

    
    qFull
}


###
# Function which computes the EP approximation of the log evidence.
#
# @param	a		The structure with the information of the problem
#
# @return	logZ		The log evidence.
#

computeEvidence_minibatch <- function(a, q, q_cavity, n_data, n_minibatch) {

    # We get the marginals of the posterior approximation
    
    log_Z <- matrix(0, a$n, a$nK)
    
    b_yi_matrix <- matrix(0, a$n, a$nK)
    
    for (i in 1 : a$nK) {
        b_yi_matrix[ which(a$Y == i), ] <- a$gMGPCinfo[[i]]$sigma + a$gMGPCinfo[[i]]$sigma0 + 1e-10    #K_xixi
        b_yi_matrix[ which(a$Y == i), ] <- b_yi_matrix[ which(a$Y == i), ] - matrix(a$gMGPCinfo[[ i ]]$KnmKmmInvKmn[ which(a$Y == i) ], sum(a$Y == i), a$nK)
    }
    
    b_yi_matrix[is.na(b_yi_matrix)] <- 0
    
    b_yi_matrix <- b_yi_matrix + q_cavity$utEu_yi
    for (i in 1 : a$nK) {
        a_yi <- q_cavity$utm_yi[ , i ]
        a_c <- q_cavity$utm_c[ , i ]
        b_yi <- b_yi_matrix[ , i ]
        b_c <- a$gMGPCinfo[[i]]$sigma + a$gMGPCinfo[[i]]$sigma0 + 1e-10     #K_xixi
        b_c <- b_c - a$gMGPCinfo[[ i ]]$KnmKmmInvKmn + q_cavity$utEu_c[ , i ]
        
        if (any(b_yi + b_c < 0)) {
            stop("Negative variances!")
        }
        
        beta <- (a_yi - a_c) / sqrt(b_yi + b_c)
        alpha <- exp(dnorm(beta,0,1, log = T) - pnorm(beta, log.p = TRUE))
        log_Z[ ,i ] <- pnorm(beta, log.p = TRUE)
        
    }
    
    log_Z[ cbind(1 : a$n, a$Y) ] <- 0
    
    log_evidence <- sum(log_Z)
    
    log_evidence <- log_evidence - 0.5 * sum(log(1 - a$f1Hat$C1_c * q_cavity$utEpostu_c))
    log_evidence <- log_evidence - 0.5 * sum(log(1 - a$f1Hat$C1_yi * q_cavity$utEpostu_yi))
    
    log_evidence <- log_evidence + 0.5 * sum(q_cavity$utmpost_c^2 * (a$f1Hat$C1_c^-1 - q_cavity$utEpostu_c)^-1)
    log_evidence <- log_evidence + 0.5 * sum(a$f1Hat$C2_c^2 * q_cavity$utEpostu_c)
    log_evidence <- log_evidence + 0.5 * sum(a$f1Hat$C2_c^2 * q_cavity$utEpostu_c^2 * (a$f1Hat$C1_c^-1 - q_cavity$utEpostu_c)^-1)
    log_evidence <- log_evidence - sum(q_cavity$utmpost_c * a$f1Hat$C2_c)
    log_evidence <- log_evidence - sum(a$f1Hat$C2_c * q_cavity$utEpostu_c * (a$f1Hat$C1_c^-1 - q_cavity$utEpostu_c)^-1 * q_cavity$utmpost_c)
    
    log_evidence <- log_evidence + 0.5 * sum(q_cavity$utmpost_yi^2 * (a$f1Hat$C1_yi^-1 - q_cavity$utEpostu_yi)^-1)
    log_evidence <- log_evidence + 0.5 * sum(a$f1Hat$C2_yi^2 * q_cavity$utEpostu_yi)
    log_evidence <- log_evidence + 0.5 * sum(a$f1Hat$C2_yi^2 * q_cavity$utEpostu_yi^2 * (a$f1Hat$C1_yi^-1 - q_cavity$utEpostu_yi)^-1)
    log_evidence <- log_evidence - sum(q_cavity$utmpost_yi * a$f1Hat$C2_yi)
    log_evidence <- log_evidence - sum(a$f1Hat$C2_yi * q_cavity$utEpostu_yi * (a$f1Hat$C1_yi^-1 - q_cavity$utEpostu_yi)^-1 * q_cavity$utmpost_yi)
    
    log_evidence <- log_evidence * n_data / n_minibatch
    
    for (i in 1 : a$nK) {
        log_evidence <- log_evidence - sum(log(diag(q[[ i ]]$cholEinv)))
        log_evidence <- log_evidence - sum(log(diag(a$gMGPCinfo[[ i ]]$cholKmm)))
        log_evidence <- log_evidence + 0.5 * sum((q[[ i ]]$cholEinv %*% q[[ i ]]$m)^2)
    }
    
    log_evidence
}

computeEvidence_test <- function(q, a, anew) {
    
    # We get the marginals of the posterior approximation
    
    log_Z <- matrix(0, a$n, a$nK)
    
    
    for (j in 1 : a$n) {
        for (i in 1 : a$nK) {
            
            Ecavyi <- solve(solve(q[[ a$Y[ j ] ]]$E) - a$f1Hat$C1_yi[ j, i ] * (a$f1Hat$u_c[ j, a$Y[ j ], ] %*% t(a$f1Hat$u_c[ j, a$Y[ j ], ])))
            mcavyi <- Ecavyi %*% (solve(q[[ a$Y[ j ] ]]$E) %*% q[[ a$Y[ j ] ]]$m - a$f1Hat$C2_yi[ j, i ] * a$f1Hat$u_c[ j, a$Y[ j ], ])
            Ecavc <- solve(solve(q[[ i ]]$E) - a$f1Hat$C1_c[ j, i ] * (a$f1Hat$u_c[ j, i, ] %*% t(a$f1Hat$u_c[ j, i, ])))
            mcavc <- Ecavc %*% (solve(q[[ i ]]$E) %*% q[[ i ]]$m - a$f1Hat$C2_c[ j, i ] * a$f1Hat$u_c[ j, i, ])
            
            a_yi <- sum(mcavyi * anew$f1Hat$u_c[ j, a$Y[ j ], ])
            a_c <- sum(mcavc * anew$f1Hat$u_c[ j, i, ])
            b_yi <- anew$gMGPCinfo[[ a$Y[ j ] ]]$sigma + anew$gMGPCinfo[[ a$Y[ j ] ]]$sigma0 + 1e-10    #K_xixi
            b_yi <- b_yi - anew$gMGPCinfo[[ a$Y[ j ] ]]$KnmKmmInvKmn[ j ] + sum(anew$gMGPCinfo[[ a$Y[ j ] ]]$KnmKmmInv[ j, ] %*% Ecavyi %*% anew$gMGPCinfo[[ a$Y[ j ] ]]$KnmKmmInv[ j, ])
            
            b_c <- anew$gMGPCinfo[[i]]$sigma + anew$gMGPCinfo[[i]]$sigma0 + 1e-10    #K_xixi
            b_c <- b_c - anew$gMGPCinfo[[ i ]]$KnmKmmInvKmn[ j ] + sum(anew$gMGPCinfo[[ i ]]$KnmKmmInv[ j, ] %*% Ecavc %*% anew$gMGPCinfo[[ i ]]$KnmKmmInv[ j, ])
            
            if (any(b_yi + b_c < 0)) {
                stop("Negative variances!")
            }
            
            beta <- (a_yi - a_c) / sqrt(b_yi + b_c)
            log_Z[ j,i ] <- pnorm(beta, log.p = TRUE)
            
        }
    }
    
    log_Z[ cbind(1 : a$n, a$Y) ] <- 0
    
    log_evidence <- sum(log_Z)
    
    log_evidence
    
} 


##
# This function computes the gradients of the ML approximation provided by EP once it has converged
#

computeGradsLogNorm <- function(a, q) {
    
    # We compute some matrices that are needed for the gradient 
    
    #q <- reconstructPosterior(a)
    
    gr_log_l  <- matrix(0, a$nK, a$d)
    
    gr_xbar <- list()
    
    gr_log_sigma0 <- gr_log_sigma <- rep(0, a$nK)
    
    for (i in 1 : a$nK) {
        
        mNew <- q[[ i ]]$m
        VNew <- q[[ i ]]$E
        
        gMGPCinfo <- a$gMGPCinfo[[ i ]]
        
        Kinv <- gMGPCinfo$KmmInv
        K <- gMGPCinfo$Kmm
        M <- Kinv - Kinv %*% (VNew %*% Kinv) - ((Kinv %*% mNew) %*% (t(mNew) %*% Kinv))
        
        # We compute the derivatives of the kernel with respect to log_sigma
        
        dKmm_dlog_sigma0 <- diag(gMGPCinfo$m) * gMGPCinfo$sigma0
        
        gr_log_sigma0[ i ] <- -0.5 * sum(M * dKmm_dlog_sigma0)
        
        # We compute the derivatives of the kernel with respect to log_sigma0
        
        dKmm_dlog_sigma <- gMGPCinfo$Kmm - diag(rep(gMGPCinfo$sigma0, gMGPCinfo$m), gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)
        
        gr_log_sigma[ i ] <- -0.5 * sum(M * dKmm_dlog_sigma)
        
        # We compute the derivatives of the kernel with respect to l
        
        #    	gr_log_l <- rep(0, length(gMGPCinfo$l))
        
        #    	for (j in 1 : length(gMGPCinfo$l)) {
        #    		distance <- as.matrix(dist(gMGPCinfo$Xbar[, j, drop = FALSE ]))^2
        #    		dKmm_dlog_l <- gMGPCinfo$Kmm * 0.5 * distance / gMGPCinfo$l[ j ]
        #    		gr_log_l[ i, j ] <- - 0.5 * sum(M * dKmm_dlog_l)
        #    	}
        
        # The distance is v^2 1^T - 2 v v^T + 1^T v^2
        
        Ml <- 0.5 * M * gMGPCinfo$Kmm
        Xl <-  (gMGPCinfo$Xbar / matrix(sqrt(gMGPCinfo$l), nrow(gMGPCinfo$Xbar), ncol(gMGPCinfo$Xbar), byrow = TRUE))
        gr_log_l[ i, ] <- (- 2 * 0.5 * colSums(Ml %*% Xl^2) + 2 * 0.5 * colSums(Xl * (Ml %*% Xl))) #* gMGPCinfo$l
        
        #	gr_xbar <- matrix(0, gMGPCinfo$m, length(l)) 
        #
        #	for (i in 1 : length(l)) {
        #
        #		distance <- (matrix(gMGPCinfo$Xbar[ , i ], gMGPCinfo$m, gMGPCinfo$m) - matrix(gMGPCinfo$Xbar[ , i ],
        #			gMGPCinfo$m, gMGPCinfo$m, byrow = T))
        #
        #		dKmm_dXbar <- - (gMGPCinfo$Kmm - diag(gMGPCinfo$sigma0, gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)) * distance / gMGPCinfo$l[ i ]
        #
        #		gr_xbar[ ,i ] <- - 0.5 * rowSums(M * dKmm_dXbar) + 0.5 * colSums(M * dKmm_dXbar)
        #	}
        
        # The derivative of the distance is dk 1^T - 1^T dk where k is the number of the pseudoinput
        # the distance is vi 1^T - 1^T vi
        Xbar <- (gMGPCinfo$Xbar / matrix(gMGPCinfo$l, nrow(gMGPCinfo$Xbar), ncol(gMGPCinfo$Xbar), byrow = TRUE))
        Mbar <- t(M) * (- (gMGPCinfo$Kmm - diag(gMGPCinfo$sigma0, gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)))
        gr_xbar[[ i ]] <- - 0.5 * 2 * (Xbar * matrix(rep(1, gMGPCinfo$m) %*% Mbar, gMGPCinfo$m, length(gMGPCinfo$l)) - Mbar %*% Xbar)
    }
    
    list(gr_log_l = gr_log_l, gr_log_sigma0 = gr_log_sigma0, gr_log_sigma = gr_log_sigma, gr_xbar = gr_xbar)
}


computeGradsLikelihood <- function(a, q, q_cavity) {
    
    # Initialization
    
    log_Z <- matrix(0, a$n, a$nK)
    a_yi <- a_c <- b_yi <- b_c <- matrix(0, a$n, a$nK)
    mOldMatrix_yi <- mOldMatrix_c <- array(0, c(a$n, a$nK, a$m))
    uVOldMatrix_yi <- uVOldMatrix_c <- array(0, c(a$n, a$nK, a$m))
    gr_log_l  <- matrix(0, a$nK, a$d)
    gr_xbar <- list()
    gr_log_sigma0 <- gr_log_sigma <- rep(0, a$nK)
    
    b_yi_matrix <- matrix(0, a$n, a$nK)
    
    for (i in 1 : a$nK) {
        b_yi_matrix[ which(a$Y == i), ] <- a$gMGPCinfo[[i]]$sigma + a$gMGPCinfo[[i]]$sigma0 + 1e-10    #K_xixi
        b_yi_matrix[ which(a$Y == i), ] <- b_yi_matrix[ which(a$Y == i), ] - 
		matrix(a$gMGPCinfo[[ i ]]$KnmKmmInvKmn[ which(a$Y == i) ], sum(a$Y == i), a$nK)
    }
    
    b_yi_matrix <- b_yi_matrix + q_cavity$utEu_yi
    
    for (i in 1 : a$nK) {
        
        gMGPCinfo <- a$gMGPCinfo[[ i ]]
        
        mOldMatrix_c[ ,  i, ] <- matrix(q[[ i ]]$m, a$n, a$m, byrow = TRUE)
        mOldMatrix_yi[ which(a$Y == i), ,  ] <- aperm(array(q[[ i ]]$m, c(a$m, sum(a$Y == i), a$nK)), c(2, 3, 1))
        Eu <- a$f1Hat$u_c[ ,i, ] %*% q[[ i ]]$E 
        value_c <- (a$f1Hat$C1_c[ ,i ]^-1 - q[[ i ]]$utEu)^-1 
        mOldMatrix_c[ ,  i, ] <- mOldMatrix_c[ ,  i, ] + Eu * matrix(value_c * (q[[ i ]]$utm - 
   		q[[ i ]]$utEu * a$f1Hat$C2_c[ ,i ]) - a$f1Hat$C2_c[ ,i ], a$n, a$m)
        uVOldMatrix_c[ , i, ] <- Eu + matrix(q[[ i ]]$utEu * value_c, a$n, a$m) * Eu
        
        value_yi_to_add <- aperm(array(Eu[ which(a$Y == i), ], c(sum(a$Y == i), a$m, a$nK)), c(1, 3, 2))
        value_yi <- (a$f1Hat$C1_yi[ which(a$Y == i), ]^-1 - matrix(q[[ i ]]$utEu[ which(a$Y == i) ], sum(a$Y == i), a$nK))^-1
        value_yi_to_multiply_1 <- value_yi * matrix(q[[ i ]]$utm[ which(a$Y == i) ], sum(a$Y == i), a$nK)
        value_yi_to_multiply_2 <- a$f1Hat$C2_yi[ which(a$Y == i), ]
        value_yi_to_multiply_3 <- value_yi * matrix(q[[ i ]]$utEu[ which(a$Y == i) ], sum(a$Y == i), a$nK) * a$f1Hat$C2_yi[ which(a$Y == i), ]
        value_yi_to_multiply <- array(value_yi_to_multiply_1 - value_yi_to_multiply_2 - value_yi_to_multiply_3, c(sum(a$Y == i), a$nK, a$m))
        
        mOldMatrix_yi[ which(a$Y == i), ,  ] <- array(mOldMatrix_yi[ which(a$Y == i), ,  ], c(sum(a$Y == i),a$nK, a$m)) + 
		value_yi_to_add * value_yi_to_multiply
        uVOldMatrix_yi[ which(a$Y == i), , ] <- value_yi_to_add +  array(matrix(q[[ i ]]$utEu[ which(a$Y == i) ], sum(a$Y == i), a$nK) * 
		value_yi, c(sum(a$Y == i), a$nK, a$m)) * value_yi_to_add
        
        a_yi[ ,i ] <- q_cavity$utm_yi[ , i ]
        a_c[ ,i ] <- q_cavity$utm_c[ , i ]
        b_yi[ ,i ] <- b_yi_matrix[ ,i ]   
        b_c[ ,i ] <- gMGPCinfo$sigma + gMGPCinfo$sigma0 + 1e-10     #K_xixi
        b_c[ ,i ] <- b_c[ ,i ] - gMGPCinfo$KnmKmmInvKmn + q_cavity$utEu_c[ , i ]
        
        if (any(b_yi[ ,i ] + b_c[ ,i ] < 0)) {
            stop("Negative variances!")
        }
        
    }
    
    # We are done with the pre-computations
    
    # Now the contribution due to the c's
    
    beta_matrix <- z_matrix <- theta_matrix <- matrix(0, a$n, a$nK)
    
    for (i in 1 : a$nK) {
        
        gMGPCinfo <- a$gMGPCinfo[[ i ]]
        
        # We precompute some derivatives
        
        dKnn_dlog_sigma0 <- rep(gMGPCinfo$sigma0, gMGPCinfo$n)
        dKmm_dlog_sigma0 <- diag(gMGPCinfo$m) * gMGPCinfo$sigma0
        dP_dlog_sigma0 <- matrix(0, gMGPCinfo$n, gMGPCinfo$m)
        
        dKnn_dlog_sigma <- rep(gMGPCinfo$sigma, gMGPCinfo$n)
        dKmm_dlog_sigma <- gMGPCinfo$Kmm - diag(rep(gMGPCinfo$sigma0, gMGPCinfo$m), gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)
        dP_dlog_sigma <- gMGPCinfo$Knm 
        
        beta <- (a_yi[ ,i ] - a_c[ ,i ]) / sqrt(b_yi[ ,i ] + b_c[ ,i ])
        beta_matrix[ , i ] <- beta
        alpha <- exp(dnorm(beta, 0, 1, log = T) - pnorm(beta, log.p = TRUE))
        log_Z[ ,i ] <- pnorm(beta, log.p = TRUE)
        
        theta <- a_yi[ ,i ] - a_c[ ,i ]
        theta_matrix[ , i ] <- theta
        z <- b_yi[ ,i ] + b_c[ ,i ]
        z_matrix[, i ] <- z
        
        alpha <-  1 / sqrt(z) * exp(dnorm(theta / sqrt(z), 0, 1, log = TRUE) - pnorm(theta / sqrt(z), log.p = TRUE))
        
        binary_flag_matrix_c <- as.integer(a$Y != i) %*% t(rep(1, a$m))
        
        dlog_evidece_dPRtR_c <- matrix(exp(-log_Z[,i ] + dnorm(theta / sqrt(z), log = TRUE)), a$n, a$m) *
            (-mOldMatrix_c[, i, ] / matrix(sqrt(z), a$n, a$m) - matrix(0.5 * theta * 1 / z^(3/2), a$n, a$m) * 
                 (-gMGPCinfo$Knm + 2 * uVOldMatrix_c[, i, ]))
        
        dlog_evidece_dP_c <- matrix(exp(-log_Z[ ,i ] + dnorm(theta / sqrt(z), log = TRUE)), a$n, a$m) * 
            (matrix(0.5 * theta * 1 / z^(3 / 2), a$n, a$m) * gMGPCinfo$KnmKmmInv) #PRtR) 
        
        dlog_evidece_dKnn_c <- exp(-log_Z[ ,i ] + dnorm(theta / sqrt(z), log = TRUE)) * (- 0.5 * theta * 1 / z^(3/2)) 
        
        # We set to zero the entries that don't belong to the class
        
        dlog_evidece_dPRtR_c <- dlog_evidece_dPRtR_c * binary_flag_matrix_c
        dlog_evidece_dKnn_c <- dlog_evidece_dKnn_c * as.integer(a$Y != i)
        dlog_evidece_dP_c <- dlog_evidece_dP_c * binary_flag_matrix_c
        
        # We now compute the actual gradients
        
        M1_c <- gMGPCinfo$KmmInv %*% t(dlog_evidece_dPRtR_c)
        M2_c <- - M1_c %*% gMGPCinfo$KnmKmmInv
        M3_c <- dlog_evidece_dP_c
        
        gr_log_sigma_c <- sum(t(M1_c) * dP_dlog_sigma) + sum(t(M2_c) * dKmm_dlog_sigma) + sum(t(dlog_evidece_dKnn_c) * dKnn_dlog_sigma) + 
            sum(M3_c * dP_dlog_sigma) 
        
        gr_log_sigma[ i ] <- gr_log_sigma_c
        
        gr_log_sigma0_c <- sum(t(M1_c) * dP_dlog_sigma0) + sum(t(M2_c) * dKmm_dlog_sigma0) + sum(t(dlog_evidece_dKnn_c) * dKnn_dlog_sigma0) +
            sum(M3_c * dP_dlog_sigma0)
        
        gr_log_sigma0[ i ] <- gr_log_sigma0_c
        
        Ml_c <- 0.5 * (t(M1_c) + M3_c) * gMGPCinfo$Knm
        Xbarl <-  (gMGPCinfo$Xbar / matrix(sqrt(gMGPCinfo$l), nrow(gMGPCinfo$Xbar), ncol(gMGPCinfo$Xbar), byrow = TRUE))
        Xl <-  (gMGPCinfo$X / matrix(sqrt(gMGPCinfo$l), nrow(gMGPCinfo$X), ncol(gMGPCinfo$X), byrow = TRUE))
        Ml2_c <- t(M2_c) * gMGPCinfo$Kmm * 0.5
        
        gr_log_l_c <- colSums(t(Ml_c) %*% Xl^2) - 2  * colSums(Xl * (Ml_c %*% Xbarl)) + colSums(Ml_c %*% Xbarl^2) + 
            colSums(t(Ml2_c) %*% Xbarl^2) - 2 * colSums(Xbarl * (Ml2_c %*% Xbarl)) + colSums(Ml2_c %*% Xbarl^2)
        
        gr_log_l[ i, ] <- gr_log_l_c
        
        Xbar <- (gMGPCinfo$Xbar / matrix(gMGPCinfo$l, nrow(gMGPCinfo$Xbar), ncol(gMGPCinfo$Xbar), byrow = TRUE))
        X <- (gMGPCinfo$X / matrix(gMGPCinfo$l, nrow(gMGPCinfo$X), ncol(gMGPCinfo$X), byrow = TRUE))
        Mbar_c <- t(M2_c) * - (gMGPCinfo$Kmm - diag(gMGPCinfo$sigma0, gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)) 
        Mbar2_c <- (t(M1_c) + M3_c) * gMGPCinfo$Knm
        
        gr_xbar_c <- (Xbar * matrix(rep(1, gMGPCinfo$m) %*% Mbar_c, gMGPCinfo$m, length(gMGPCinfo$l)) - t(Mbar_c) %*% Xbar) +
            (Xbar * matrix(rep(1, gMGPCinfo$m) %*% t(Mbar_c), gMGPCinfo$m, length(gMGPCinfo$l)) - Mbar_c %*% Xbar) + 
            (t(Mbar2_c) %*% X) - ((Xbar * matrix(rep(1, gMGPCinfo$n) %*% Mbar2_c, gMGPCinfo$m, length(gMGPCinfo$l))))
        
        gr_xbar[[ i ]] <- gr_xbar_c
        
    }
    
    # Now the contribution due to the yi
    
    for (i in 1 : a$nK) {
        
        gMGPCinfo <- a$gMGPCinfo[[ i ]]
        
        index <- which(a$Y == i)
        n_total <- sum(a$Y == i) * (a$nK - 1)
        
        if (n_total != 0) {

            # We precompute some derivatives
            
            dKnn_dlog_sigma0 <- rep(gMGPCinfo$sigma0, sum(a$Y == i) * (a$nK - 1))
            dKmm_dlog_sigma0 <- diag(gMGPCinfo$m) * gMGPCinfo$sigma0
            dP_dlog_sigma0 <- matrix(0, sum(a$Y == i) * (a$nK - 1), gMGPCinfo$m)
            
            dKnn_dlog_sigma <- rep(gMGPCinfo$sigma, sum(a$Y == i) * (a$nK - 1))
            dKmm_dlog_sigma <- gMGPCinfo$Kmm - diag(rep(gMGPCinfo$sigma0, gMGPCinfo$m), gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)
            dP_dlog_sigma <- gMGPCinfo$Knm[ rep(index, a$nK - 1), ]
            
            if (any(b_yi[ index, -i ] + b_c[ index, -i ] < 0)) {
                stop("Negative variances!")
            }
            
            log_Z_yi <- c(log_Z[ index, -i ])
            theta <- c(theta_matrix[ index, -i ])
            z <- c(z_matrix[ index, -i ])
            
            mOldMatrix_yi_actual <- mOldMatrix_yi[ index, -i, ]
            dim(mOldMatrix_yi_actual) <- c(n_total, a$m)
            uVOldMatrix_yi_actual <- uVOldMatrix_yi[ index, -i, ]
            dim(uVOldMatrix_yi_actual) <- c(n_total, a$m)
            
            dlog_evidece_dPRtR_yi <- matrix(exp(-log_Z_yi + dnorm(theta / sqrt(z), log = TRUE)), n_total, a$m) *
                (mOldMatrix_yi_actual / matrix(sqrt(z), n_total, a$m) - matrix(0.5 * theta * 1 / z^(3/2), n_total, a$m) * 
                     (-gMGPCinfo$Knm[ rep(index, a$nK - 1),  ] + 2 * uVOldMatrix_yi_actual))
            
            dlog_evidece_dP_yi <- matrix(exp(-log_Z_yi + dnorm(theta / sqrt(z), log = TRUE)), n_total, a$m) * 
                (matrix(0.5 * theta * 1 / z^(3 / 2), n_total, a$m) * gMGPCinfo$KnmKmmInv[ rep(index, a$nK - 1), ]) #PRtR) 
            
            dlog_evidece_dKnn_yi <- exp(-log_Z_yi + dnorm(theta / sqrt(z), log = TRUE)) * (- 0.5 * theta * 1 / z^(3/2)) 
            
            # We now compute the actual gradients
            
            M1_yi <- gMGPCinfo$KmmInv %*% t(dlog_evidece_dPRtR_yi)
            M2_yi <- - M1_yi %*% gMGPCinfo$KnmKmmInv[ rep(index, a$nK - 1), ]
            M3_yi <- dlog_evidece_dP_yi
            
            gr_log_sigma_yi <- sum(t(M1_yi) * dP_dlog_sigma) + sum(t(M2_yi) * dKmm_dlog_sigma) + sum(t(dlog_evidece_dKnn_yi) * dKnn_dlog_sigma) + 
                sum(M3_yi * dP_dlog_sigma) 
            
            gr_log_sigma[ i ] <- gr_log_sigma[ i ] + gr_log_sigma_yi
            
            gr_log_sigma0_yi <- sum(t(M1_yi) * dP_dlog_sigma0) + sum(t(M2_yi) * dKmm_dlog_sigma0) + sum(t(dlog_evidece_dKnn_yi) * dKnn_dlog_sigma0) +
                sum(M3_yi * dP_dlog_sigma0)
            
            gr_log_sigma0[ i ] <- gr_log_sigma0[ i ] + gr_log_sigma0_yi 
            
            Ml_yi <- 0.5 * (t(M1_yi) + M3_yi) * gMGPCinfo$Knm[ rep(index, a$nK - 1), ]
            Xbarl <-  (gMGPCinfo$Xbar / matrix(sqrt(gMGPCinfo$l), nrow(gMGPCinfo$Xbar), ncol(gMGPCinfo$Xbar), byrow = TRUE))
            Xl <-  (gMGPCinfo$X[ rep(index, a$nK - 1), ] / matrix(sqrt(gMGPCinfo$l), n_total, ncol(gMGPCinfo$X), byrow = TRUE))
            Ml2_yi <- t(M2_yi) * gMGPCinfo$Kmm * 0.5
            
            gr_log_l_yi <- colSums(t(Ml_yi) %*% Xl^2) - 2  * colSums(Xl * (Ml_yi %*% Xbarl)) + colSums(Ml_yi %*% Xbarl^2) + 
                colSums(t(Ml2_yi) %*% Xbarl^2) - 2 * colSums(Xbarl * (Ml2_yi %*% Xbarl)) + colSums(Ml2_yi %*% Xbarl^2)
            
            gr_log_l[ i, ] <- gr_log_l[ i, ] + gr_log_l_yi
            
            Xbar <- (gMGPCinfo$Xbar / matrix(gMGPCinfo$l, nrow(gMGPCinfo$Xbar), ncol(gMGPCinfo$Xbar), byrow = TRUE))
            X <- (gMGPCinfo$X / matrix(gMGPCinfo$l, nrow(gMGPCinfo$X), ncol(gMGPCinfo$X), byrow = TRUE))
            Mbar_yi <- t(M2_yi) * - (gMGPCinfo$Kmm - diag(gMGPCinfo$sigma0, gMGPCinfo$m) - diag(1e-10, gMGPCinfo$m)) 
            Mbar2_yi <- (t(M1_yi) + M3_yi) * gMGPCinfo$Knm[ rep(index, a$nK - 1), ]
            
            gr_xbar_yi <- (Xbar * matrix(rep(1, gMGPCinfo$m) %*% Mbar_yi, gMGPCinfo$m, length(gMGPCinfo$l)) - t(Mbar_yi) %*% Xbar) +
                (Xbar * matrix(rep(1, gMGPCinfo$m) %*% t(Mbar_yi), gMGPCinfo$m, length(gMGPCinfo$l)) - Mbar_yi %*% Xbar) + 
                (t(Mbar2_yi) %*% X[ rep(index, a$nK - 1), ]) - ((Xbar * matrix(rep(1, n_total) %*% Mbar2_yi, gMGPCinfo$m, length(gMGPCinfo$l))))
            
            gr_xbar[[ i ]] <- gr_xbar[[ i ]] + gr_xbar_yi
        }
        
    }
    
    list(gr_log_l = gr_log_l, gr_log_sigma0 = gr_log_sigma0, gr_log_sigma = gr_log_sigma, gr_xbar = gr_xbar)
    
}


#########################################################################################################
#
# The main function that should be used in this module is "epMGPC". You have to call it with the arguments:
#
# 	X  -> Design matrix for the classification problem.
# 	Xbar  -> Matrix of Pseudo inputs.
#   Y  -> Target vector for the classification problem. (must be an R factor)
#   sigma -> variance of the latent noise  (value)
#   sigma0 -> parameter for the gaussian kernel (value)
#   l -> length scale parameter for the gaussian kernel (value)
#   start -> initial EP approximation (NULL if not used) It should be "a".
#
# returns the approximate distribution for the posterior as a list with several components.
#
epMGPCInternal <- function(X, Y, m, n_minibatch, Xbar_ini = NULL, log_sigma = rep(0, length(levels(Y))), 
                           log_sigma0 = rep(0, length(levels(Y))), log_l = matrix(0, length(levels(Y)), ncol(X)), 
                           max_iters = 250, start = NULL, X_test = NULL, Y_test = NULL, optim = 'adam',
                           print_interval = ceiling(nrow(X)/n_minibatch)) {
    
    t0 <<- proc.time()
    
    # We initialize the hyper-parameters 
    
    sigma <- exp(log_sigma * 1.0)
    sigma0 <- exp(log_sigma0 * 1.0)
    l <- exp(log_l * 1.0)
    
    nK <- length(levels(Y))
    n <- nrow(X)
    d <- ncol(X)
    levelsY <- levels(Y)
    Y <- as.integer(Y)
    
    nBatches <- n / n_minibatch
    
    Xbar <- list()
    for (i in 1:nK) {
        if (is.null(Xbar_ini)) {
            samples <- sample(1:n,  m)
            Xbar[[ i ]] <- X[ samples, ]
        } else {
            Xbar[[ i ]] <- Xbar_ini[[ i ]]
        }
        
        if (is.null(log_l)) {
            l[ i, ] <- rep(estimateL(Xbar[[ i ]]), d)
        }
    }

    perm <- sample(1 : n)
    X <- X[ perm, ]
    Y <- Y[ perm ]

    # We initialize the approximate factors  reusing any previous solution
    # u: Vector u for the ith example
    # Ci: Scalar C for the ith natural parameter
    
    if (!is.null(start))  {
        f1Hat <- start$f1Hat
    } else { 
        f1Hat <- list(C1_yi = matrix(0, n, nK), C2_yi = matrix(0, n, nK),
                      u_c = array(0,c(n,nK,m)), C1_c = matrix(0, n, nK), C2_c = matrix(0, n, nK))
    }
    
    
    # Optimization parameters
    
    if (optim == "adadelta") {
        xbarTemp <- list()
        for (j in 1:nK)
           xbarTemp[[ j ]] <- matrix(0, m, d)
        gms <- list(sigma = rep(0, nK), sigma0 = rep(0, nK), l = matrix(0, nK, d), xbar = xbarTemp)
        sms <- list(sigma = rep(0, nK), sigma0 = rep(0, nK), l = matrix(0, nK, d), xbar = xbarTemp)
        step <- list(sigma = rep(0, nK), sigma0 = rep(0, nK), l = matrix(0, nK, d), xbar = xbarTemp)
        
        optim_params <<- list(step_rate = 1, decay = 0.9, eps = 1e-5, gms = gms, sms = sms, step = step)
    }
    else {
        optim <- "adam"
        xbarTemp <- list()
        for (j in 1:nK)
            xbarTemp[[ j ]] <- matrix(0, m, d)
        mt <- list(sigma = rep(0, nK), sigma0 = rep(0, nK), l = matrix(0, nK, d), xbar = xbarTemp)
        vt <- list(sigma = rep(0, nK), sigma0 = rep(0, nK), l = matrix(0, nK, d), xbar = xbarTemp)
        step <- 0
        
        optim_params <<- list(alpha = 1e-3, beta1 = 0.9, beta2 = 0.999, eps = 1e-8, mt = mt, vt = vt, step = step)
    }
    
    damping <- 1
    sigma_old <- sigma
    sigma0_old <- sigma0
    l_old <- l
    Xbar_old <- list()
    
    q <- list()
    gMGPCinfo  <- list()
    for (i in 1:nK) {
        gMGPCinfo[[ i ]] <- initialize_kernel(Xbar[[ i ]], Xbar[[ i ]], sigma[ i ], sigma0[ i ], l[ i, ])
        q[[ i ]] <- initialize_approx(gMGPCinfo[[i]]$Kmm, gMGPCinfo[[i]]$KmmInv)
    }
    
    # Main loop of EP
    
    i <- 1
    count_minibatches <- 0
    while (i < max_iters) {
        
        avg_evidence <- 0
        perm <- sample(1 : nrow(X))
        X <- X[ perm, ]
        Y <- Y[ perm ]
        f1Hat$u_c <- array(f1Hat$u_c[ perm,, ], c(nrow(X),nK,m))
        f1Hat$C1_c <- f1Hat$C1_c[ perm, ]
        f1Hat$C2_c <- f1Hat$C2_c[ perm, ]
        f1Hat$C1_yi <- f1Hat$C1_yi[ perm, ]
        f1Hat$C2_yi <- f1Hat$C2_yi[ perm, ]
        
        for (n in sample(1:ceiling(nBatches))) {

            cat(".")
            
            to_sel <- (round((n - 1) * nrow(X) / nBatches) + 1) : pmin(round(exp(log(n) + log(nrow(X))) / nBatches), nrow(X))
            
            Xbatch <- X[ to_sel,, drop = FALSE ]
            Ybatch <- Y[ to_sel ]
            f1HatMinibatch <- list(C1_yi = f1Hat$C1_yi[to_sel,], C2_yi = f1Hat$C2_yi[to_sel,],
                                   u_c = array(f1Hat$u_c[to_sel,,], c(length(to_sel),nK,m)), C1_c = f1Hat$C1_c[to_sel,], C2_c = f1Hat$C2_c[to_sel,])
            
            
            # We initialize the approximations (one for each class)
            
            u <- array(0, c(length(to_sel), nK, m))
            
            for (j in 1:nK) {
                gMGPCinfo[[ j ]] <- initialize_kernel(Xbatch, Xbar[[ j ]], sigma[ j ], sigma0[ j ], l[ j, ])
                u[ ,j, ] <- gMGPCinfo[[ j ]]$KnmKmmInv
            }

            a <- list(gMGPCinfo = gMGPCinfo, f1Hat = f1HatMinibatch, nK = nK, n = length(to_sel), Y = Ybatch, m = m, d = d)
            
            C1_c_old <- f1HatMinibatch$C1_c
            C1_yi_old <- f1HatMinibatch$C1_yi 
            C2_c_old <- f1HatMinibatch$C2_c
            C2_yi_old <- f1HatMinibatch$C2_yi 
            uold <- f1HatMinibatch$u_c
            
            # Constructs the cavity distribution
            
            q_cavity <- cavityDistribution_minibatch(a$Y, a$nK, a$n, q, uold, u, C1_c_old, C1_yi_old, C2_c_old, C2_yi_old)
            
            avg_evidence <- avg_evidence + computeEvidence_minibatch(a, q, q_cavity, nrow(X), nrow(Xbatch))

            # Refines factors
            
            a <- refineFactors(a, q_cavity, damping)
            
            f1Hat$C1_c[to_sel,]  <- a$f1Hat$C1_c
            f1Hat$C1_yi[to_sel,] <- a$f1Hat$C1_yi
            f1Hat$C2_c[to_sel,] <- a$f1Hat$C2_c
            f1Hat$C2_yi[to_sel,] <- a$f1Hat$C2_yi
            f1Hat$u_c[to_sel,,] <- a$f1Hat$u_c
         
            # We get the marginals of the posterior approximation
            
            q <- updatePosterior(a$Y, a$nK, a$n, a$m, q, uold, u, a$f1Hat$C1_c, a$f1Hat$C1_yi, 
                                 a$f1Hat$C2_c, a$f1Hat$C2_yi, C1_c_old, C1_yi_old, C2_c_old, C2_yi_old)

            sigma_old <- sigma
            sigma0_old <- sigma0
            l_old <- l
            Xbar_old <- Xbar
            KmmInv_old <- list()
            for (j in 1:a$nK)
                KmmInv_old[[ j ]] <- a$gMGPCinfo[[ j ]]$KmmInv
            
            # Update of the hyper-parameters
            
            q_cavity <- cavityDistribution_minibatch(a$Y, a$nK, a$n, q, u, u, a$f1Hat$C1_c, a$f1Hat$C1_yi, a$f1Hat$C2_c, a$f1Hat$C2_yi)
            
            a$gMGPCinfo <- optimize(a, q, q_cavity, nrow(X), nrow(Xbatch), optim)
            
            KmmInv <- list()
            for (j in 1:a$nK) {
                KmmInv[[ j ]] <- a$gMGPCinfo[[ j ]]$KmmInv
                sigma[ j ] <- a$gMGPCinfo[[ j ]]$sigma
                sigma0[ j ] <- a$gMGPCinfo[[ j ]]$sigma0
                l[ j, ] <- a$gMGPCinfo[[ j ]]$l
                Xbar[[ j ]] <- a$gMGPCinfo[[ j ]]$Xbar
            }
            
            # Update posterior after change of the hyper-parameters

            q <- updatePosterior_after(a$Y, a$nK, a$n, a$m, q, u, KmmInv_old, KmmInv)
            
            if ((REPORT == TRUE) && ((count_minibatches %% print_interval) == 0)) {
                
                # Reconstruct structure
                a$f1Hat <- f1Hat
                a$X <- X
                a$Y <- Y
                a$n <- nrow(X)
                
                ret <- 	list(f1Hat = a$f1Hat, gMGPCinfo = a$gMGPCinfo, 
                             n = n, nK = nK, Y = Y, a = a, levelsY = levelsY)
                
                t_before <- proc.time()
                performance <- evaluate_test_performance(ret = ret, X_test = X_test, Y_test = Y_test, q = q)
                t_after <- proc.time()
                
                
                t0 <- t0 + (t_after - t_before)
                
                write.table(t(c(performance$err, performance$neg_ll, proc.time() - t0)),
                            file = paste("./results/time_outter_", CONT, ".txt", sep = ""), row.names = F, col.names = F, append = TRUE)
            }
            
            count_minibatches <- count_minibatches + 1
            
            
        }
        
        cat("\n")
        
        cat("Epoch",  i, "Avg evidence:", avg_evidence / nBatches, "\n")
        
        # Annealed damping scheme
        
        # damping <- damping * 0.99
        
        i <- i + 1
    }
    
    
    a$f1Hat <- f1Hat
    a$X <- X
    a$Y <- Y
    for (j in 1:nK) {
        a$gMGPCinfo[[ j ]] <- initialize_kernel(X, Xbar[[ j ]], sigma[ j ], sigma0[ j ], l[ j, ])
    }
    n <- a$n <- nrow(X)
    
    list(f1Hat = a$f1Hat, gMGPCinfo = a$gMGPCinfo, n = n, nK = nK, Y = Y, a = a, levelsY = levelsY)
}



##
# Function that estimates the initial lengthscale value

estimateL <- function (X) {
    
    D <- as.matrix(dist(X))
    median(D[ upper.tri(D) ])
}




##
# Function that evaluates the performance of a prediction

evaluate_test_performance <- function(ret, X_test, Y_test, q = NULL) {
    
    if (is.factor(Y_test)) Y_test <- as.integer(Y_test)
    
    n <- nrow(X_test)
    v_k <- v_k <- matrix(0, n, ret$nK)
    m_k <- m_k <- matrix(0, n, ret$nK)
    
    # To make prediction we need the means and the variances of the process at the new point.
    # These are obtained by computing int p(f*|f) p(f|D) d f, where p(f|D) is the posterior
    # approximation computed. Importantly Kf*f (the covariance between the observed points and the
    # new points) 
    
    if (is.null(q))
        q <- reconstructPosterior(ret$a)
    
    means <- variances <- matrix(0, n, ret$a$nK)
    
    for (i in 1 : ret$nK) {
        
        gMGPCinfo <- ret$gMGPCinfo[[ i ]]
        Knm <- kernel_nm(X_test, gMGPCinfo$Xbar, gMGPCinfo$l, gMGPCinfo$sigma)
        KnmKinv <- Knm %*% gMGPCinfo$KmmInv
        
        means[ ,i ] <- KnmKinv %*% q[[ i ]]$m
        variances[ ,i ] <- gMGPCinfo$sigma + gMGPCinfo$sigma0 + 1e-10 - rowSums(KnmKinv * Knm) + 
            rowSums((KnmKinv %*% q[[ i ]]$E) * KnmKinv)
        
    }
    
    # We compute the label of each instance using quadrature
    
    labels <- rep(ret$Y[ 1 ], n)
    prob <- matrix(0, n, length(ret$levelsY))
    colnames(prob) <- ret$levelsY
    
    n_points_grid <- 100
    grid <- matrix(seq(-6, 6, length = n_points_grid), n, n_points_grid, byrow = TRUE)
    
    index <- matrix(FALSE, n, ret$nK)
    index[ cbind(1:n, Y_test) ] <- TRUE
    means_class <- means[ cbind(1:n, Y_test) ]
    variances_class <- variances[ cbind(1:n, Y_test) ]
    means_not_class <- t(matrix(c(t(means))[c(t(!index))], ret$nK - 1, n))
    variances_not_class <- t(matrix(c(t(variances))[c(t(!index))], ret$nK - 1, n))
    
    grid_k <- grid * sqrt(matrix(variances_class, n, n_points_grid)) + matrix(means_class, n, n_points_grid)
    values_k <- dnorm(grid_k, mean = matrix(means_class, n, n_points_grid), sd = sqrt(matrix(variances_class, n, n_points_grid)), log = TRUE)
    
    for (j in 1 : (ret$nK - 1)) {
        values_k <- values_k + pnorm((grid_k - matrix(means_not_class[ , j ], n, n_points_grid)) / 
                                         sqrt(matrix(variances_not_class[ , j ], n, n_points_grid)), log.p = TRUE)
    }
    
    # index <- matrix(TRUE, n, ret$nK)
    # index[ cbind(1:n, Y_test) ] <- FALSE
    # means_index <- t(matrix(c(t(means))[c(t(index))], ret$nK - 1, n))
    # variances_index <- t(matrix(c(t(variances))[c(t(index))], ret$nK - 1, n))
    # 
    # grid_means <- apply(means_index, 2, function(m) grid_k - matrix(m, n, n_points_grid))
    # sqrt_variances <- apply(variances_index, 2, function(v) sqrt(matrix(v, n, n_points_grid)))
    # to_add <- pnorm(grid_means / sqrt_variances, log.p = TRUE)
    # 
    # val <- values_k
    # values_k_2 <- matrix(apply(to_add, 2, function(add) val <- val + add), n, n_points_grid)
    
    prob_class <- rowSums(exp(values_k)) * (grid_k[ ,2 ] - grid_k[ ,1 ])
    
    
    err <- mean(apply(means, 1, which.max) != Y_test)
    neg_ll <- -mean(log(prob_class))
    
    
    list(err = err, neg_ll = neg_ll)
    
    
}


##
#                                                          
# Function that initalizes the approximation to the prior distribution
#                  
#   @param  Kmm Prior covariance function 
#   @param  KmmInv  Inverse of the prior covariance function  
#   @param  u_c Vector u
#
#   @return q   Structure containing the approximation  
#                                                          

initialize_approx <- function(Kmm, KmmInv) {
    
    q <- list()
    
    q$Einv <- KmmInv
    q$cholEinv <- chol(KmmInv)
    q$E <- Kmm
    q$m <- rep(0, nrow(Kmm))
    
    q
}


##
# Function that initializes the struture with the problem information.
#
# @param	X	n x d matrix with the data points.
# @param	Xbar	m x d matrix with the pseudo inputs.
# @param	sigma	scalar with the log-amplitude of the GP.
# @param	sigma0	scalar with the log-noise level in the GP.
# @param	l	d-dimensional vector with the log-lengthscales.
# 
# @return	gMGPCinfo	List with the problem information 
#

initialize_kernel <- function(X, Xbar, sigma, sigma0, l) {
    
    # We initialize the structure with the data and the kernel
    # hyper-parameters
    
    gMGPCinfo <- list()
    gMGPCinfo$X <- X
    gMGPCinfo$Xbar <- Xbar
    gMGPCinfo$m <- nrow(Xbar)
    gMGPCinfo$n <- nrow(X)
    gMGPCinfo$sigma <- sigma
    gMGPCinfo$sigma0 <- sigma0
    gMGPCinfo$l <- l
    
    # We compute the kernel matrices
    gMGPCinfo$Kmm <- kernel(Xbar, gMGPCinfo$l, gMGPCinfo$sigma0, gMGPCinfo$sigma)
    gMGPCinfo$Knm <- kernel_nm(X, Xbar, gMGPCinfo$l, gMGPCinfo$sigma)
    gMGPCinfo$cholKmm <- chol(gMGPCinfo$Kmm)
    gMGPCinfo$KmmInv <- chol2inv(gMGPCinfo$cholKmm)
    gMGPCinfo$KnmKmmInv <- gMGPCinfo$Knm %*% gMGPCinfo$KmmInv
    gMGPCinfo$KnmKmmInvKmn <- rowSums(gMGPCinfo$KnmKmmInv * gMGPCinfo$Knm)
    
    gMGPCinfo
}




##
# This function computes the covariance matrix for the GP
#

kernel <- function(X, l, sigma0, sigma) {
    X <- X / matrix(sqrt(l), nrow(X), ncol(X), byrow = TRUE)

#    distance <- as.matrix(dist(X)) ^ 2

    value <- rowSums(X^2)
    Q <- matrix(value, nrow(X), nrow(X))
    distance <- Q + t.default(Q) - 2 * X %*% t.default(X)

    sigma * exp(-0.5 * distance) + diag(sigma0, nrow(X)) + diag(rep(1e-10, nrow(X)), nrow(X))
}

##
# Function which computes the kernel matrix between the observed data and the test data
#

kernel_nm <- function(X, Xnew, l, sigma) {
    X <- X / matrix(sqrt(l), nrow(X), ncol(X), byrow = TRUE)	
    Xnew <- Xnew / matrix(sqrt(l), nrow(Xnew), ncol(Xnew), byrow = TRUE)
    n <- nrow(X)
    m <- nrow(Xnew)
    Q <- matrix(apply(X ^ 2, 1, sum), n, m)
    Qbar <- matrix(apply(Xnew ^ 2, 1, sum), n, m, byrow = T)
    distance <- Qbar + Q - 2 * X %*% t(Xnew)
    sigma * exp(-0.5 * distance)
}





##
# This function updates the kernel hyper-parameters using ADAM or ADADELTA optimization method
#

optimize <- function(a, q, q_cavity, n_data, n_minibatch, optim) {

    gMGPCinfo  <- a$gMGPCinfo
    
    sigma <- sigma0 <- rep(0, a$nK)
    l <- matrix(0, a$nK, a$d)
    Xbar <- list()
    
    for (i in 1 : a$nK) {
        sigma[ i ] <- gMGPCinfo[[ i ]]$sigma
        sigma0[ i ] <- gMGPCinfo[[ i ]]$sigma0
        l[ i, ] <- gMGPCinfo[[ i ]]$l
        Xbar[[ i ]] <- gMGPCinfo[[ i ]]$Xbar
    }
    
    
    # We obtain the gradients

    grad_prior <- computeGradsLogNorm(a, q)
    grad_likelihood <- computeGradsLikelihood(a, q, q_cavity)
    
    grad <- grad_prior
    grad$gr_log_l <- grad$gr_log_l + n_data / n_minibatch * grad_likelihood$gr_log_l
    grad$gr_log_sigma <- grad$gr_log_sigma + n_data / n_minibatch * grad_likelihood$gr_log_sigma
    grad$gr_log_sigma0 <- grad$gr_log_sigma0 + n_data / n_minibatch * grad_likelihood$gr_log_sigma0
    
    for (i in 1 : a$nK) {
        grad$gr_xbar[[ i ]] <- grad$gr_xbar[[ i ]] + n_data / n_minibatch * grad_likelihood$gr_xbar[[ i ]]  
    }
    
    if (optim == "adadelta") {
        optim_params$gms$sigma0 <<- optim_params$decay * optim_params$gms$sigma0 + (1 - optim_params$decay) * grad$gr_log_sigma0^2
        optim_params$gms$sigma <<- optim_params$decay * optim_params$gms$sigma + (1 - optim_params$decay) * grad$gr_log_sigma^2
        optim_params$gms$l <<- optim_params$decay * optim_params$gms$l + (1 - optim_params$decay) * grad$gr_log_l^2
        for (i in 1 : a$nK) { 
            optim_params$gms$xbar[[ i ]] <<- optim_params$decay * optim_params$gms$xbar[[ i ]] + (1 - optim_params$decay) * grad$gr_xbar[[ i ]]^2
        }
        
        optim_params$step$sigma0 <<- sqrt(optim_params$sms$sigma0 + optim_params$eps) / sqrt(optim_params$gms$sigma0 + optim_params$eps) * grad$gr_log_sigma0 * optim_params$step_rate
        optim_params$step$sigma <<- sqrt(optim_params$sms$sigma + optim_params$eps) / sqrt(optim_params$gms$sigma + optim_params$eps) * grad$gr_log_sigma * optim_params$step_rate
        optim_params$step$l <<- sqrt(optim_params$sms$l + optim_params$eps) / sqrt(optim_params$gms$l + optim_params$eps) * grad$gr_log_l * optim_params$step_rate
        for (i in 1 : a$nK) 
            optim_params$step$xbar[[ i ]] <<- sqrt(optim_params$sms$xbar[[ i ]] + optim_params$eps) / sqrt(optim_params$gms$xbar[[ i ]] + optim_params$eps) * grad$gr_xbar[[ i ]] * optim_params$step_rate
        
        l <- exp(log(l) + optim_params$step$l)
        sigma0 <- exp(log(sigma0) + optim_params$step$sigma0)
        sigma <- exp(log(sigma) + optim_params$step$sigma)
        for (i in 1 : a$nK) 
            Xbar[[ i ]] <- Xbar[[ i ]] + optim_params$step$xbar[[ i ]]
        
        optim_params$sms$sigma0 <<- optim_params$decay * optim_params$sms$sigma0 + (1 - optim_params$decay) * optim_params$step$sigma0^2
        optim_params$sms$sigma <<- optim_params$decay * optim_params$sms$sigma + (1 - optim_params$decay) * optim_params$step$sigma^2
        optim_params$sms$l <<- optim_params$decay * optim_params$sms$l + (1 - optim_params$decay) * optim_params$step$l^2
        for (i in 1 : a$nK)
            optim_params$sms$xbar[[ i ]] <<- optim_params$decay * optim_params$sms$xbar[[ i ]] + (1 - optim_params$decay) * optim_params$step$xbar[[ i ]]^2
    }
    else {


        optim_params$step <<- optim_params$step + 1
        
        optim_params$mt$sigma0 <<- optim_params$beta1 * optim_params$mt$sigma0 + (1 - optim_params$beta1) * grad$gr_log_sigma0
        optim_params$mt$sigma <<- optim_params$beta1 * optim_params$mt$sigma + (1 - optim_params$beta1) *  grad$gr_log_sigma
        optim_params$mt$l <<- optim_params$beta1 * optim_params$mt$l + (1 - optim_params$beta1) * grad$gr_log_l
        for (i in 1 : a$nK) { 
            optim_params$mt$xbar[[ i ]] <<- optim_params$beta1 * optim_params$mt$xbar[[ i ]] + (1 - optim_params$beta1) * grad$gr_xbar[[ i ]]
        }
        
        optim_params$vt$sigma0 <<- optim_params$beta2 * optim_params$vt$sigma0 + (1 - optim_params$beta2) * grad$gr_log_sigma0^2
        optim_params$vt$sigma <<- optim_params$beta2 * optim_params$vt$sigma + (1 - optim_params$beta2) *  grad$gr_log_sigma^2
        optim_params$vt$l <<- optim_params$beta2 * optim_params$vt$l + (1 - optim_params$beta2) * grad$gr_log_l^2
        for (i in 1 : a$nK) { 
            optim_params$vt$xbar[[ i ]] <<- optim_params$beta2 * optim_params$vt$xbar[[ i ]] + (1 - optim_params$beta2) * grad$gr_xbar[[ i ]]^2
        }
        
        
        mtHat <- optim_params$mt
        mtHat$sigma0 <- optim_params$mt$sigma0 / (1 - optim_params$beta1^optim_params$step )
        mtHat$sigma <- optim_params$mt$sigma / (1 - optim_params$beta1^optim_params$step )
        mtHat$l <- optim_params$mt$l / (1 - optim_params$beta1^optim_params$step)
        for (i in 1 : a$nK) { 
            mtHat$xbar[[ i ]] <- optim_params$mt$xbar[[ i ]] / (1 - optim_params$beta1^optim_params$step)
        }
        
        vtHat <- optim_params$vt
        vtHat$sigma0 <- optim_params$vt$sigma0 / (1 - optim_params$beta2^optim_params$step )
        vtHat$sigma <- optim_params$vt$sigma / (1 - optim_params$beta2^optim_params$step )
        vtHat$l <- optim_params$vt$l / (1 - optim_params$beta2^optim_params$step)
        for (i in 1 : a$nK) { 
            vtHat$xbar[[ i ]] <- optim_params$vt$xbar[[ i ]] / (1 - optim_params$beta2^optim_params$step)
        }
        
        step_size <- optim_params$mt
        step_size$sigma0 <- optim_params$alpha * mtHat$sigma0 / (sqrt(vtHat$sigma0) + optim_params$eps)
        step_size$sigma <- optim_params$alpha * mtHat$sigma / (sqrt(vtHat$sigma) + optim_params$eps)
        step_size$l <- optim_params$alpha * mtHat$l / (sqrt(vtHat$l) + optim_params$eps)
        for (i in 1 : a$nK) { 
            step_size$xbar[[ i ]] <- optim_params$alpha * mtHat$xbar[[ i ]] / (sqrt(vtHat$xbar[[ i ]]) + optim_params$eps)
        }
        
        l <- exp(log(l) + step_size$l)
        sigma0 <- exp(log(sigma0) + step_size$sigma0)
        sigma <- exp(log(sigma) + step_size$sigma)
        for (i in 1 : a$nK) 
            Xbar[[ i ]] <- Xbar[[ i ]] + step_size$xbar[[ i ]]
    }
    
    
    for (i in 1 : a$nK) {
        gMGPCinfo[[ i ]] <- initialize_kernel(gMGPCinfo[[ i ]]$X, Xbar[[ i ]], sigma[ i ], sigma0[ i ], l[ i, ])
    }

    gMGPCinfo
}



###
# Function which computes the probability of class 1 on new data.
#
# @param	ret	The list returned by epGPCExternal
# @param	Xtest	The n x d matrix with the new data points.
#
# @return	pOne	The probability of class 1 on the new data.
#

predictMGPC <- function(ret, X_test, q = NULL) {
  
  n <- nrow(X_test)
  v_k <- v_k <- matrix(0, n, ret$nK)
  m_k <- m_k <- matrix(0, n, ret$nK)
  
  # To make prediction we need the means and the variances of the process at the new point.
  # These are obtained by computing int p(f*|f) p(f|D) d f, where p(f|D) is the posterior
  # approximation computed. Importantly Kf*f (the covariance between the observed points and the
  # new points) 
  
  if (is.null(q))
      q <- reconstructPosterior(ret$a)
  
  means <- variances <- matrix(0, n, ret$a$nK)
  
  for (i in 1 : ret$nK) {
    
    gMGPCinfo <- ret$gMGPCinfo[[ i ]]
    Knm <- kernel_nm(X_test, gMGPCinfo$Xbar, gMGPCinfo$l, gMGPCinfo$sigma)
    KnmKinv <- Knm %*% gMGPCinfo$KmmInv
    
    means[ ,i ] <- KnmKinv %*% q[[ i ]]$m
    variances[ ,i ] <- gMGPCinfo$sigma + gMGPCinfo$sigma0 + 1e-10 - rowSums(KnmKinv * Knm) + 
      rowSums((KnmKinv %*% q[[ i ]]$E) * KnmKinv)
    
  }

  # We compute the label of each instance using quadrature
  
  labels <- rep(ret$Y[ 1 ], n)
  prob <- matrix(0, n, length(ret$levelsY))
  colnames(prob) <- ret$levelsY

  n_points_grid <- 100
  grid <- matrix(seq(-6, 6, length = n_points_grid), n, n_points_grid, byrow = TRUE)
  
  for (k in 1 : ret$nK) {
      
     grid_k <- grid * sqrt(matrix(variances[ , k ], n, n_points_grid)) + matrix(means[ , k ], n, n_points_grid)

     values_k <- dnorm(grid_k, mean = matrix(means[ , k ], n, n_points_grid), 
        sd = sqrt(matrix(variances[ , k ], n, n_points_grid)), log = TRUE)

     for (j in 1 : ret$nK) {
       if (j != k) {
          values_k <- values_k + pnorm((grid_k - matrix(means[ , j ], n, n_points_grid)) / 
            sqrt(matrix(variances[ , j ], n, n_points_grid)), log.p = TRUE)
       }
     }

     pTmp <- rowSums(exp(values_k)) * (grid_k[ ,2 ] - grid_k[ ,1 ])
     
     prob[ , k ] <- pTmp 
  }

  maximums <- apply(prob, 1, which.max)
  labels <- ret$levelsY[ maximums ]
  maxProb <- prob[ cbind(1 : n, maximums) ]
 
  list(labels = labels, prob = prob, maxprob = maxProb)
}


###
# Function which reconstructs the posterior factor q
#
# @param	a 		The approximation 
#
# @return	q 	reconstructed factor
#
#

reconstructPosterior <- function(a){
    
    q <- list()
    
    for (i in 1 : a$nK) {
        
        qTemp <- list()
        
        # Examples corresponding with the class i
        samples_i <- a$Y == i
        
        # We add the natural parameters of the prior
        
        M <- a$gMGPCinfo[[ i ]]$KmmInv
        
        # We sum over the natural parameters of the likelihood factors
        
        M <- M + t(matrix(rowSums(a$f1Hat$C1_yi[ samples_i, , drop = FALSE]), sum(samples_i), a$m) * 
                       a$f1Hat$u_c[ samples_i, i,  ]) %*% a$f1Hat$u_c[ samples_i, i,  ] 
        M <- M + t(matrix(a$f1Hat$C1_c[ !samples_i, i ], a$n - sum(samples_i), a$m) * 
                       a$f1Hat$u_c[ !samples_i, i,  ]) %*% a$f1Hat$u_c[ !samples_i, i,  ] 
        
        qTemp$Einv <- M
        qTemp$cholEinv <- chol(M)
        qTemp$E <- chol2inv(qTemp$cholEinv)
        s <- colSums(matrix(rowSums(a$f1Hat$C2_yi[ samples_i, , drop = FALSE]), sum(samples_i), a$m) * a$f1Hat$u_c[ samples_i, i,  ])
        s <- s + colSums(matrix(a$f1Hat$C2_c[ !samples_i, i ], a$n - sum(samples_i), a$m) * a$f1Hat$u_c[ !samples_i, i,  ])
        qTemp$m <- qTemp$E %*% s
        
        qTemp$utEu <- rowSums((a$f1Hat$u_c[ , i,  ] %*% qTemp$E) * a$f1Hat$u_c[ , i,  ])
        qTemp$utm <- a$f1Hat$u_c[ , i,  ] %*% qTemp$m
        
        q[[ i ]] <- qTemp
    }
    
    q
}


###
# Function that refines the factors by matching their moments
#
#
refineFactors <- function(a, q_cavity, damping){
    
    b_yi_matrix <- matrix(0, a$n, a$nK)
    
    for (i in 1 : a$nK) {
        b_yi_matrix[ which(a$Y == i), ] <- a$gMGPCinfo[[i]]$sigma + a$gMGPCinfo[[i]]$sigma0 + 1e-10    #K_xixi
        b_yi_matrix[ which(a$Y == i), ] <- b_yi_matrix[ which(a$Y == i), ] - matrix(a$gMGPCinfo[[ i ]]$KnmKmmInvKmn[ which(a$Y == i) ], sum(a$Y == i), a$nK)
    }
    
    b_yi_matrix <- b_yi_matrix + q_cavity$utEu_yi
    
    for (i in 1 : a$nK) {
        
        a_yi <- q_cavity$utm_yi[ , i ]
        a_c <- q_cavity$utm_c[ , i ]
        b_yi <- b_yi_matrix[ , i ]
        b_c <- a$gMGPCinfo[[i]]$sigma + a$gMGPCinfo[[i]]$sigma0 + 1e-10     #K_xixi
        b_c <- b_c - a$gMGPCinfo[[ i ]]$KnmKmmInvKmn + q_cavity$utEu_c[ , i ]
        
        if (any(b_yi + b_c < 0)) {
            stop("Negative variances!")
        }
        
        beta <- (a_yi - a_c) / sqrt(b_yi + b_c)
        alpha <- exp(dnorm(beta,0,1, log = T) - pnorm(beta, log.p = TRUE))
        Z <- pnorm(beta)
        
        C1_yi_i <- 1 / (1.0 / ((alpha^2 + alpha * beta) / (b_yi + b_c)) - q_cavity$utEu_yi[ , i ])
        C1_c_i <- 1 / (1.0 / ((alpha^2 + alpha * beta) / (b_yi + b_c)) - q_cavity$utEu_c[ , i ])
        C2_yi_i <-  alpha / sqrt(b_yi + b_c) * (1 + C1_yi_i * q_cavity$utEu_yi[ , i ]) + C1_yi_i * a_yi 
        C2_c_i <- alpha / sqrt(b_yi + b_c) * (-1 - C1_c_i * q_cavity$utEu_c[ , i ]) + C1_c_i * a_c 
        
        
        # We do the update 
        
        a$f1Hat$C1_yi[ , i ]  <- damping * C1_yi_i  + (1 - damping) * a$f1Hat$C1_yi[ , i ]
        a$f1Hat$C2_yi[ , i ]  <- damping * C2_yi_i  + (1 - damping) * a$f1Hat$C2_yi[ , i ]
        a$f1Hat$C1_c[ , i ]  <- damping * C1_c_i  + (1 - damping) * a$f1Hat$C1_c[ , i ]
        a$f1Hat$C2_c[ , i ]  <- damping * C2_c_i  + (1 - damping) * a$f1Hat$C2_c[ , i ]
        
        #    if (i == 2)
        #    	browser()
        #	ret <- checkMoments(aOld, q, 1, 2, a, alpha, beta, b_yi, b_c)
        
        
        # browser()
        a$f1Hat$u_c[, i, ] <- a$gMGPCinfo[[ i ]]$KnmKmmInv
        
    }
    
    # We set to zero the parameters that do not appear in the likelihood
    
    a$f1Hat$C1_c[ cbind(1 : a$n, a$Y) ] <- 0
    a$f1Hat$C2_c[ cbind(1 : a$n, a$Y) ] <- 0
    a$f1Hat$C1_yi[ cbind(1 : a$n, a$Y) ] <- 0
    a$f1Hat$C2_yi[ cbind(1 : a$n, a$Y) ] <- 0
    
    a
}



###
# Function which updates the posterior factor q after a change in the u vector
#
# @param	a 		The approximation 
#
# @return	q 	reconstructed factor
#
#

updatePosterior <- function(Y, nK, n, m, q, uold, u, C1_c, C1_yi, C2_c, C2_yi, C1_c_old, C1_yi_old, C2_c_old, C2_yi_old){

    for (i in 1 : nK) {
        Einv <- q[[i]]$Einv
        
        uold_yi <- matrix(uold[which(Y == i), i, ], sum(Y == i), m)
        u_yi <- matrix(u[which(Y == i), i, ], sum(Y == i), m)
        
        Einv <- Einv - t(uold[, i, ]) %*% (matrix(C1_c_old[ ,i ], n, m) * uold[, i, ])
        Einv <- Einv - t(uold_yi) %*% (matrix(rowSums(C1_yi_old[ which(Y == i), , drop = FALSE ]), sum(Y == i), m) * uold_yi)
        Einv <- Einv + t(u[, i, ]) %*% (C1_c[ ,i ] * u[, i, ])
        Einv <- Einv + t(u_yi) %*% (rowSums(C1_yi[which(Y == i), , drop = FALSE]) * u_yi)
        
        Enew <- chol2inv(chol(Einv))
        m_nat <- q[[i]]$Einv %*% q[[i]]$m
        m_nat <- m_nat - t(uold[, i, ]) %*% C2_c_old[ , i ]  - t(uold_yi) %*% rowSums(C2_yi_old[ which(Y == i), , drop = FALSE ])
        m_nat <- m_nat + t(u[, i, ]) %*% C2_c[ , i ] + t(u_yi) %*% rowSums(C2_yi[ which(Y == i), , drop = FALSE ])
        
        mnew <- Enew %*% m_nat
        q[[i]]$E <- Enew
        q[[i]]$Einv <- Einv
        q[[i]]$cholEinv <- chol(Einv)
        q[[i]]$m <- mnew
        q[[i]]$utEu <- rowSums((u[, i, ] %*% q[[i]]$E) * u[, i, ])
        q[[i]]$utm <- u[, i, ] %*% q[[i]]$m
    }
    
    q
}


###
# Function which updates the posterior factor q after a change in the hyper-parameters
#
# @param	a 		The approximation 
#
# @return	q 	reconstructed factor
#
#

updatePosterior_after <- function(Y, nK, n, m, q, u, KmmInv_old, KmmInv)  {
    for (i in 1:nK) {
        Einv <- q[[i]]$Einv
        Einvnew <- Einv - KmmInv_old[[ i ]] + KmmInv[[ i ]]
        Enew <- chol2inv(chol(Einvnew))
        mnew <- Enew %*% (Einv %*% q[[i]]$m)
        q[[i]]$E <- Enew
        q[[i]]$Einv <- Einvnew
        q[[i]]$cholEinv <- chol(Einvnew)
        q[[i]]$m <- mnew
        q[[i]]$utEu <- rowSums((u[, i, ] %*% q[[i]]$E) * u[, i, ])
        q[[i]]$utm <- u[, i, ] %*% q[[i]]$m
    }
    
    q
}


