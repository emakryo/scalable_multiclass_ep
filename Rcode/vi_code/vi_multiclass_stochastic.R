

jitter <- 1e-10

t0 <- NULL

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

##
# Function that estimates the initial lengthscale value

estimateL <- function (X) {
	D <- as.matrix(dist(X))
	median(D[ upper.tri(D) ])
}

##
# This function computes the covariance matrix for the GP
#

kernel <- function(X, l, sigma0, sigma) {

        X <- X / matrix(sqrt(l), nrow(X), ncol(X), byrow = TRUE)

#        distance <- as.matrix(dist(X))^2

	value <- rowSums(X^2)
	Q <- matrix(value, nrow(X), nrow(X))
	distance <- Q + t.default(Q) - 2 * X %*% t.default(X)
	
        sigma * exp(-0.5 * distance) + diag(sigma0, nrow(X)) + diag(rep(jitter, nrow(X)))
}

##
# Function which computes the kernel matrix between the observed data and the test data
#

kernel_nm <- function(X, Xnew, l, sigma) {

        X <- X / matrix(sqrt(l), nrow(X), ncol(X), byrow = TRUE)
        Xnew <- Xnew / matrix(sqrt(l), nrow(Xnew), ncol(Xnew), byrow = TRUE)
        n <- nrow(X)
        m <- nrow(Xnew)
        Q <- matrix(apply(X^2, 1, sum), n, m)
        Qbar <- matrix(apply(Xnew^2, 1, sum), n, m, byrow = T)
        distance <- Qbar + Q - 2 * X %*% t.default(Xnew)
        sigma * exp(-0.5 * distance)
}


##
# Function that produces upper and lower bounds for the encoding of the
# parameters in a vector
#

getBounds_function <- function(n_pseudo_inputs, d) {

    u_log_sigma <- 10
    l_log_sigma <- -10

    u_log_sigma_0 <- 10
    l_log_sigma_0 <- -10

    u_log_l <- rep(10, d)
    l_log_l <- rep(-10, d)

    u_m <- rep(Inf, n_pseudo_inputs)
    l_m <- rep(-Inf, n_pseudo_inputs)

    u_L <- matrix(Inf, n_pseudo_inputs, n_pseudo_inputs)
    l_L <- matrix(-Inf, n_pseudo_inputs, n_pseudo_inputs)

    u_pseudo_inputs <- matrix(Inf, n_pseudo_inputs, d)
    l_pseudo_inputs <- matrix(-Inf, n_pseudo_inputs, d)

    weights <- c()
    upper <- c(u_m, as.vector(u_L[ lower.tri(u_L, diag = T) ]),
        as.vector(u_pseudo_inputs), u_log_sigma, u_log_sigma_0, u_log_l)
    lower <- c(l_m, as.vector(l_L[ lower.tri(l_L, diag = T) ]),
        as.vector(l_pseudo_inputs), l_log_sigma, l_log_sigma_0, l_log_l)

    list(upper = upper, lower = lower)
}

getBounds <- function(n_pseudo_inputs, d, K) {

	result <- getBounds_function(n_pseudo_inputs, d)
	
	for (i in 1 : (K - 1)) {
		result_tmp <- getBounds_function(n_pseudo_inputs, d)
		result$upper <- c(result_tmp$upper)	
		result$lower <- c(result_tmp$lower)	
	}

	result
}


##
# Function that encodes the parameters in a vector
#

encodeParameters_function <- function(m, L, pseudo_inputs, log_sigma, log_sigma_0, log_l) {

    weights <- c()
    c(m, as.vector(L[ lower.tri(L, diag = T) ]), as.vector(pseudo_inputs),
        log_sigma, log_sigma_0, log_l)
}

encodeParameters <- function(m, L, pseudo_inputs, log_sigma, log_sigma_0, log_l, K) {

	result <- NULL
	
	for (i in 1 : K)  {
		result <- c(result, encodeParameters_function(m[ , i ], L[[ i ]], pseudo_inputs[[ i ]], 
			log_sigma[ i ], log_sigma_0[ i ], log_l[ , i ]))
		
	}

	result
}

##
# Function that decodes the parameters from a vector
#
# @arg weights          The vector of weights encoding the parameters.
# @arg n_pseudo_inputs  The number of pseudo inputs.
# @arg d                The dimensionality of the data.
#

decodeParameters_function <- function(weights, n_pseudo_inputs, d) {

    counter <- 0
    m <- weights[ (counter + 1) : (counter + n_pseudo_inputs) ]
    counter <- counter + n_pseudo_inputs
    L <- matrix(0, n_pseudo_inputs, n_pseudo_inputs)
    increment <- n_pseudo_inputs * (n_pseudo_inputs - 1) / 2 + n_pseudo_inputs
    L[ lower.tri(L, diag = T) ] <-
        weights[ (counter + 1) : (counter + increment) ]
    counter <- counter + increment
    pseudo_inputs <- matrix(weights[ (counter + 1) : (counter + d *
        n_pseudo_inputs) ], n_pseudo_inputs, d)
    counter <- counter + d * n_pseudo_inputs
    log_sigma <- weights[ counter + 1 ]
    log_sigma_0 <- weights[ counter + 2 ]
    log_l <- weights[ (counter + 2 + 1) : (counter + 2 + d) ]

    list(m = m, L = L, pseudo_inputs = pseudo_inputs, log_sigma = log_sigma,
        log_sigma_0 = log_sigma_0, log_l = log_l)
}


decodeParameters <- function(weights, n_pseudo_inputs, d, K) {

	n_params_per_function  <- n_pseudo_inputs + n_pseudo_inputs * (n_pseudo_inputs - 1) / 2 + n_pseudo_inputs + 
		d * n_pseudo_inputs + 1 + 1 + d
	m <- log_sigma_0 <- log_sigma <- log_l <- NULL
	L <- list()
	pseudo_inputs <- list()

	for (i in 1 : K) {

		result <- decodeParameters_function(weights[ (1 : n_params_per_function) + n_params_per_function * (i - 1) ], 
			n_pseudo_inputs, d)
		
		m <- cbind(m, result$m)
		L[[ i ]] <- result$L
		pseudo_inputs[[ i ]] <- result$pseudo_inputs
		log_sigma_0 <- c(log_sigma_0, result$log_sigma_0)
		log_sigma <- c(log_sigma, result$log_sigma)
		log_l <- cbind(log_l, result$log_l)
	}

	list(m = m, L = L, pseudo_inputs = pseudo_inputs, log_sigma = log_sigma,
		log_sigma_0 = log_sigma_0, log_l = log_l)
}

##
# Function that evaluates the variational lower bound.
#

lower_bound <- function(weights, X, y, n_pseudo_inputs, eps = 0.001, d = ncol(X), K, scale_factor_likelihood = 1.0) {

	ret <- decodeParameters(weights, n_pseudo_inputs, d, K)
	m <- ret$m
	L <- ret$L
	pseudo_inputs <- ret$pseudo_inputs
	log_sigma <- ret$log_sigma
	log_sigma_0 <- ret$log_sigma_0
	log_l <- ret$log_l

	# We compute the kernel matrices

	KL <- 0
	means <- vars <- matrix(0, nrow(X), K)

	for (i in 1 : K) {
	    
		Kmm <- kernel(pseudo_inputs[[ i ]], exp(log_l[, i ]), exp(log_sigma_0[ i ]), exp(log_sigma[ i ]))
		cholKmmT <- chol(Kmm)
		cholKmm <- t.default(cholKmmT)
		KmmInv <- chol2inv(cholKmmT)
		Knm <- kernel_nm(X, pseudo_inputs[[ i ]], exp(log_l[, i ]), exp(log_sigma[ i ]))

		# We compute the KL terms

		KmmInvL <- inverse_times(cholKmm, L[[ i ]])

		KL <- KL + 0.5 * (sum(2 * log(diag(cholKmm))) - sum(2 * log(abs(diag(L[[ i ]])))) - n_pseudo_inputs + 
			sum(KmmInvL * L[[ i ]]) + sum(forwardsolve(cholKmm, m[ , i ])^2))

		# We compute the marginal means and variances

		A <- t.default(inverse_times(cholKmm, t.default(Knm)))
		means[ , i ] <- (Knm %*% (KmmInv %*% m[, i ]))
		vars[ , i ] <- exp(log_sigma[ i ]) + exp(log_sigma_0[ i ]) + jitter + rowSums((A %*% L[[ i ]])^2) - rowSums(A * Knm)
	}

	# We compute the integrals

	choose_y <- cbind(1 : nrow(X), y)

	means_integrand_y <- means[ choose_y ]
	vars_integrand_y <- vars[ choose_y ]
	bool_index <- matrix(TRUE, nrow(X), K)
	bool_index[ choose_y ] <- FALSE

	means_k <- matrix(t.default(means)[ t.default(bool_index) ], nrow(X), K - 1, byrow = TRUE)
	vars_k <- matrix(t.default(vars)[ t.default(bool_index) ], nrow(X), K - 1, byrow = TRUE)

        n_points_grid <- 100
	grid <- matrix(seq(-10, 10, length = n_points_grid), nrow(X), n_points_grid, byrow = TRUE) *  
		matrix(sqrt(vars_integrand_y), nrow(X), n_points_grid) + matrix(means_integrand_y, nrow(X), n_points_grid)

	integrand <- dnorm(grid, matrix(means_integrand_y, nrow(X), n_points_grid), matrix(sqrt(vars_integrand_y), 
		nrow(X), n_points_grid), log = TRUE)

	for (i in 1 : (K - 1)) {
		mean_matrix <- matrix(means_k[, i ], nrow(X), n_points_grid)
		var_matrix <- matrix(vars_k[, i ], nrow(X), n_points_grid)
		integrand <- integrand + pnorm((grid - mean_matrix) / sqrt(var_matrix), log.p = TRUE)
	}

	probs <- rowSums((exp(integrand[, -1 ]) + exp(integrand[, -ncol(integrand)])) / 2) * (grid[, 2 ] - grid[ , 1 ])

	ret <- log((1 - eps) * 1 + eps / K) * probs + (1 - probs) * log(eps / K)
	
#	cat("Elog:", mean(ret), "Elog scaled:", scale_factor_likelihood * sum(ret), "KL:", KL, "Sum:", scale_factor_likelihood * sum(ret) - KL, "\n")

	scale_factor_likelihood * sum(ret) - KL
}

##
# Function that computes the gradientes of the KL term of the lower bound.
#

grad_KL_function <- function(m, L, pseudo_inputs, log_sigma, log_sigma_0, log_l, X, y, n_pseudo_inputs, d) {

    # We compute the kernel matrices
    
    Kmm <- kernel(pseudo_inputs, exp(log_l), exp(log_sigma_0), exp(log_sigma))
    cholKmmT <- chol(Kmm)
    cholKmm <- t.default(cholKmmT)
    KmmInv <- chol2inv(chol(Kmm))

    # We compute the gradient of KL with respect to m, L and Kmm

    dKLdm <- KmmInv %*% m
    KmmInvL <- inverse_times(cholKmm, L)
    dKLdKmm <- 0.5 * (KmmInv - KmmInvL %*% t.default(KmmInvL) - dKLdm %*% t.default(dKLdm))
    
#    dKLdS <- 0.5 * (KmmInv - chol2inv(t.default(L)))
#    dKLdL <- 2 * dKLdS %*% L

    dKLdL <- KmmInvL - backsolve(t.default(L), diag(n_pseudo_inputs))

    # We compute the gradient with respect to m and L

    grad_m <- - dKLdm
    grad_L <- - dKLdL
    
    # We compute the gradient with respect to log_sigma

    dKmmdlog_sigma <- Kmm - diag(exp(log_sigma_0) + jitter, n_pseudo_inputs)
    grad_log_sigma <- - sum(dKLdKmm * dKmmdlog_sigma)

    # We compute the gradients with respect to log_sigma_0

    dKmmdlog_sigma_0 <- diag(n_pseudo_inputs) * exp(log_sigma_0)
    grad_log_sigma_0 <- - sum(dKLdKmm * dKmmdlog_sigma_0)

    # We compute the gradient with respect to log_l

    grad_log_l <- rep(0, d)
    Kmm_no_diagonal <- (Kmm - diag(exp(log_sigma_0) + jitter, n_pseudo_inputs))

    Ml2 <- (- dKLdKmm) * 0.5 * Kmm_no_diagonal 
    Xbarl <-  (pseudo_inputs / matrix(sqrt(exp(log_l)), nrow(pseudo_inputs), ncol(pseudo_inputs), byrow = TRUE))
    grad_log_l <- colSums(t(Ml2) %*% Xbarl^2) - 2 * colSums(Xbarl * (Ml2 %*% Xbarl)) + colSums(Ml2 %*% Xbarl^2)

    # We compute the gradient with respec to the pseudo-inputs

    grad_pseudo_inputs <- matrix(0, n_pseudo_inputs, d)

    Xbar <- (pseudo_inputs / matrix(exp(log_l), nrow(pseudo_inputs), ncol(pseudo_inputs), byrow = TRUE))
    Mbar <- - (- dKLdKmm) * Kmm_no_diagonal
    grad_pseudo_inputs <- (Xbar * matrix(rep(1, nrow(pseudo_inputs)) %*% Mbar, nrow(pseudo_inputs), length(log_l)) - t.default(Mbar) %*% Xbar) +
	(Xbar * matrix(rep(1, nrow(pseudo_inputs)) %*% t.default(Mbar), nrow(pseudo_inputs), length(log_l)) - Mbar %*% Xbar) 

    list(grad_m = grad_m * - 1.0, grad_L = grad_L * -1.0, grad_pseudo_inputs = grad_pseudo_inputs * -1.0, 
	grad_log_sigma = grad_log_sigma * -1.0, grad_log_sigma_0 = grad_log_sigma_0 * -1.0, grad_log_l = grad_log_l * -1.0)
}

grad_KL <- function(weights, X, y, n_pseudo_inputs, d, K) {

	ret <- decodeParameters(weights, n_pseudo_inputs, d, K)
	grad_m <- grad_log_sigma <- grad_log_sigma_0 <- grad_log_l <- NULL
	grad_L <- list()
	grad_pseudo_inputs <- list()

	for (i in 1 : K) {

		result <- grad_KL_function(ret$m[ , i ], ret$L[[ i ]], ret$pseudo_inputs[[ i ]], ret$log_sigma[ i ], 
			ret$log_sigma_0[ i ], ret$log_l[ , i ], X, y, n_pseudo_inputs, d)

		grad_m <- cbind(grad_m, result$grad_m)
		grad_log_sigma <- c(grad_log_sigma, result$grad_log_sigma)
		grad_log_sigma_0 <- c(grad_log_sigma_0, result$grad_log_sigma_0)
		grad_log_l <- cbind(grad_log_l, result$grad_log_l)
		grad_L[[ i ]] <- result$grad_L
		grad_pseudo_inputs[[ i ]] <- result$grad_pseudo_inputs
	}

	encodeParameters(grad_m, grad_L, grad_pseudo_inputs, grad_log_sigma, grad_log_sigma_0, grad_log_l, K)
}

grad_likelihood <- function(weights, X, y, n_pseudo_inputs, eps, d, K) {

	ret <- decodeParameters(weights, n_pseudo_inputs, d, K)
	m <- ret$m
	L <- ret$L
	pseudo_inputs <- ret$pseudo_inputs
	log_sigma <- ret$log_sigma
	log_sigma_0 <- ret$log_sigma_0
	log_l <- ret$log_l

	# We compute the kernel matrices

	KL <- 0
	means <- vars <- matrix(0, nrow(X), K)

	Kmm_list <- cholKmmT_list <- KmmInv_list <- Knm_list <- A_list <- list()

	for (i in 1 : K) {
	    
		Kmm <- kernel(pseudo_inputs[[ i ]], exp(log_l[, i ]), exp(log_sigma_0[ i ]), exp(log_sigma[ i ]))
		cholKmmT <- chol(Kmm)
		cholKmm <- t.default(cholKmmT)
		KmmInv <- chol2inv(cholKmmT)
		Knm <- kernel_nm(X, pseudo_inputs[[ i ]], exp(log_l[, i ]), exp(log_sigma[ i ]))
		A <- t.default(inverse_times(cholKmm, t.default(Knm)))

		Kmm_list[[ i ]] <- Kmm
		cholKmmT_list[[ i ]] <- cholKmmT
		KmmInv_list[[ i ]] <- KmmInv
		Knm_list[[ i ]] <- Knm
		A_list[[ i ]] <- A

		means[ , i ] <- (Knm %*% (KmmInv %*% m[, i ]))
		vars[ , i ] <- exp(log_sigma[ i ]) + exp(log_sigma_0[ i ]) + jitter + rowSums((A %*% L[[ i ]])^2) - rowSums(A * Knm)
	}

	# We compute the integrals

	choose_y <- cbind(1 : nrow(X), y)

	means_integrand_y <- means[ choose_y ]
	vars_integrand_y <- vars[ choose_y ]

	bool_index <- matrix(TRUE, nrow(X), K)
	bool_index[ choose_y ] <- FALSE

	means_k <- matrix(t.default(means)[ t.default(bool_index) ], nrow(X), K - 1, byrow = TRUE)
	vars_k <- matrix(t.default(vars)[ t.default(bool_index) ], nrow(X), K - 1, byrow = TRUE)

        n_points_grid <- 100
	grid <- matrix(seq(-10, 10, length = n_points_grid), nrow(X), n_points_grid, byrow = TRUE) *  
		matrix(sqrt(vars_integrand_y), nrow(X), n_points_grid) + matrix(means_integrand_y, nrow(X), n_points_grid)

	integrand <- dnorm(grid, matrix(means_integrand_y, nrow(X), n_points_grid), matrix(sqrt(vars_integrand_y), 
		nrow(X), n_points_grid), log = TRUE)

	for (i in 1 : (K - 1)) {
		mean_matrix <- matrix(means_k[, i ], nrow(X), n_points_grid)
		var_matrix <- matrix(vars_k[, i ], nrow(X), n_points_grid)
		integrand <- integrand + pnorm((grid - mean_matrix) / sqrt(var_matrix), log.p = TRUE)
	}

	dprob_d_v <- dprob_d_m <- matrix(0, nrow(X), K)
	n <- nrow(X)

	for (i in 1 : K) {

		mean_matrix <- matrix(means[, i ], nrow(X), n_points_grid)
		var_matrix <- matrix(vars[, i ], nrow(X), n_points_grid)

		to_mul_not_class <- exp(dnorm((grid - mean_matrix) / sqrt(var_matrix), log = TRUE) - 
			pnorm((grid - mean_matrix) / sqrt(var_matrix), log.p = TRUE) - 0.5 * log(var_matrix)) * -1.0
		to_mul_class <- (grid - mean_matrix) / var_matrix 

		res_not_class <- rowSums((exp(integrand[, -1 ]) * to_mul_not_class[ , -1 ] + 
			exp(integrand[, -ncol(integrand) ]) * to_mul_not_class[ , -ncol(to_mul_not_class) ]) / 2) * (grid[, 2 ] - grid[ , 1 ])

		res_class <- rowSums((exp(integrand[, -1 ]) * to_mul_class[ , -1 ] + 
			exp(integrand[, -ncol(integrand) ]) * to_mul_class[ , -ncol(to_mul_class) ]) / 2) * (grid[, 2 ] - grid[ , 1 ])

		dprob_d_m[ which(y != i), i ] <- res_not_class[ which(y != i) ] 
		dprob_d_m[ which(y == i), i ] <- res_class[ which(y == i) ]

		to_mul_not_class <- exp(dnorm((grid - mean_matrix) / sqrt(var_matrix), log = TRUE) - 
			pnorm((grid - mean_matrix) / sqrt(var_matrix), log.p = TRUE)) * (grid - mean_matrix) / var_matrix^(3 / 2) * -0.5
		to_mul_class <- -0.5 * var_matrix^-1 + 0.5 * (grid - mean_matrix)^2 / var_matrix^2

		res_not_class <- rowSums((exp(integrand[, -1 ]) * to_mul_not_class[ , -1 ] + 
			exp(integrand[, -ncol(integrand) ]) * to_mul_not_class[ , -ncol(to_mul_not_class) ]) / 2) * (grid[, 2 ] - grid[ , 1 ])

		res_class <- rowSums((exp(integrand[, -1 ]) * to_mul_class[ , -1 ] + 
			exp(integrand[, -ncol(integrand) ]) * to_mul_class[ , -ncol(to_mul_class) ]) / 2) * (grid[, 2 ] - grid[ , 1 ])

		dprob_d_v[ which(y != i), i ] <- res_not_class[ which(y != i) ]
		dprob_d_v[ which(y == i), i ] <- res_class[ which(y == i) ]
	}


	dprob_d_m <- dprob_d_m * log((1 - eps) * 1.0 + eps / K) + log(eps / K) * dprob_d_m * -1.0
	dprob_d_v <- dprob_d_v * log((1 - eps) * 1.0 + eps / K) + log(eps / K) * dprob_d_v * -1.0

	# We compute the gradients of the expecatations with respect to m, L, Kmm, Knm and kii

	grad_m <- grad_log_sigma <- grad_log_sigma_0 <- grad_log_l <- NULL
	grad_L <- grad_pseudo_inputs <- list()

	for (i in 1 : K) {

		dExpdmu <- dprob_d_m[ , i ]
		dExpdv <- dprob_d_v[ , i ]

		Kmm <- Kmm_list[[ i ]] 
		cholKmmT <- cholKmmT_list[[ i ]] 
		cholKmm <- t.default(cholKmmT)
		KmmInv <- KmmInv_list[[ i ]] 
		Knm <- Knm_list[[ i ]] 
		A <- A_list[[ i ]] 

		KmmInvL <- inverse_times(cholKmm, L[[ i ]])

                dmudm <- t.default(A)
                aux_m <- matrix(dExpdmu, n_pseudo_inputs, n, byrow = T) * dmudm
                aux_v <- matrix(dExpdv, n_pseudo_inputs, n, byrow = T) * dmudm
                dExpdm <- rowSums(aux_m)
                dExpdS <- aux_v %*% t.default(dmudm)
                dExpdL <- 2 * dExpdS %*% L[[ i ]]
                dExpdmudmudKmm <- (matrix(-aux_m %*% rep(1, n), n_pseudo_inputs,n_pseudo_inputs) * t.default(matrix(KmmInv %*% m[ , i ], 
                    	n_pseudo_inputs, n_pseudo_inputs,)))
                dExpdvdvdKmm <- - (2 * KmmInvL %*% t.default(KmmInvL) - KmmInv) %*% t.default(aux_v %*% Knm)
                dExpdKmm <- dExpdmudmudKmm + dExpdvdvdKmm
   
                dExpdkii <- sum(dExpdv)
                dExpdvdvdKnm <- 2 * t.default((KmmInvL %*% t.default(L[[ i ]])) %*% aux_v - aux_v)
   
                dExpdmudmudKnm <- matrix(KmmInv %*% m[ , i ], n, n_pseudo_inputs, byrow = T) * matrix(dExpdmu, n, n_pseudo_inputs)
   
                dExpdKnm <- dExpdmudmudKnm + dExpdvdvdKnm
   
                # We compute the gradient with respect to m and L
   
                grad_m <- cbind(grad_m, dExpdm)
                grad_L[[ i ]] <- dExpdL 
                
                # We compute the gradient with respect to log_sigma
   
                dKmmdlog_sigma <- Kmm - diag(exp(log_sigma_0[ i ]) + jitter, n_pseudo_inputs)
                dKnmdlog_sigma <- Knm 
                dkiidlog_sigma <- exp(log_sigma[ i ])
   
                grad_log_sigma <- c(grad_log_sigma, sum(dExpdKmm * dKmmdlog_sigma) + sum(dExpdKnm * dKnmdlog_sigma) + dExpdkii * dkiidlog_sigma)
   
                # We compute the gradients with respect to log_sigma_0
   
                dKmmdlog_sigma_0 <- diag(n_pseudo_inputs) * exp(log_sigma_0[ i ])
                dkiidlog_sigma_0 <- exp(log_sigma_0[ i ])
                grad_log_sigma_0 <- c(grad_log_sigma_0, sum(dExpdKmm * dKmmdlog_sigma_0) + dExpdkii * dkiidlog_sigma_0)
   
                # We compute the gradient with respect to log_l
   
                Kmm_no_diagonal <- (Kmm - diag(exp(log_sigma_0[ i ]) + jitter, n_pseudo_inputs))
   
                Ml2 <- (dExpdKmm) * 0.5 * Kmm_no_diagonal 
                Ml <-  dExpdKnm * Knm * 0.5
                Xbarl <-  (pseudo_inputs[[ i ]] / matrix(sqrt(exp(log_l[ , i ])), nrow(pseudo_inputs[[ i ]]), 
			ncol(pseudo_inputs[[ i ]]), byrow = TRUE))
                Xl <-  (X / matrix(sqrt(exp(log_l[ , i ])), nrow(X), ncol(X), byrow = TRUE))
                grad_log_l <- cbind(grad_log_l, colSums(t(Ml) %*% Xl^2) - 2  * colSums(Xl * (Ml %*% Xbarl)) +  colSums(Ml %*% Xbarl^2) +
                            colSums(t(Ml2) %*% Xbarl^2) - 2 * colSums(Xbarl * (Ml2 %*% Xbarl)) + colSums(Ml2 %*% Xbarl^2))
   
                # We compute the gradient with respect to the pseudo-inputs
   
                Xbar <- (pseudo_inputs[[ i ]] / matrix(exp(log_l[ , i ]), nrow(pseudo_inputs[[ i ]]), ncol(pseudo_inputs[[ i ]]), byrow = TRUE))
                X_comp <- (X / matrix(exp(log_l[ , i ]), nrow(X), ncol(X), byrow = TRUE))
                Mbar <- - (dExpdKmm) * Kmm_no_diagonal
                Mbar2 <- dExpdKnm * Knm
                grad_pseudo_inputs[[ i ]] <- (Xbar * matrix(rep(1, nrow(pseudo_inputs[[ i ]])) %*% Mbar, 
			nrow(pseudo_inputs[[ i ]]), length(log_l[ , i ])) - 
			t.default(Mbar) %*% Xbar) + (Xbar * matrix(rep(1, nrow(pseudo_inputs[[ i ]])) %*% 
			t.default(Mbar), nrow(pseudo_inputs[[ i ]]), length(log_l[ , i ])) - 
			Mbar %*% Xbar)  + (t.default(Mbar2) %*% X_comp) - ((Xbar * matrix(rep(1, nrow(X)) %*% Mbar2, 
			nrow(pseudo_inputs[[ i ]]), length(log_l[ , i ]))))
	}

	encodeParameters(grad_m, grad_L, grad_pseudo_inputs, grad_log_sigma, grad_log_sigma_0, grad_log_l, K)
}

##
# Function that computes the gradientes of the variational lower bound.
#

grad_lower_bound <- function(weights, X, y, n_pseudo_inputs, eps, d, K, scale_factor_likelihood = 1.0) {
	-1.0 * grad_KL(weights, X, y, n_pseudo_inputs, d, K) + grad_likelihood(weights, X, y, n_pseudo_inputs, eps, d, K) * 
		scale_factor_likelihood
}


predictMGPC <- function(ret, Xtest) {
  
	n <- nrow(Xtest)
	means <- variances <- matrix(0, n, ret$K)

	for (i in 1 : ret$K) {
    
		Kmm <- kernel(ret$pseudo_inputs[[ i ]], exp(ret$log_l[, i ]), exp(ret$log_sigma_0[ i ]), exp(ret$log_sigma[ i ]))
		cholKmmT <- chol(Kmm)
		cholKmm <- t.default(cholKmmT)
		KmmInv <- chol2inv(cholKmmT)
		Knm <- kernel_nm(Xtest, ret$pseudo_inputs[[ i ]], exp(ret$log_l[, i ]), exp(ret$log_sigma[ i ]))
		A <- Knm %*% KmmInv

		means[ ,i ] <- A %*% ret$m[, i ]
		variances[ ,i ] <- exp(ret$log_sigma[ i ]) + exp(ret$log_sigma_0[ i ]) + 
			ret$jitter + rowSums((A %*% ret$L[[ i ]])^2 - (A %*% cholKmm)^2)
	}
	
	# We compute the label of each instance using quadrature
  
	labels <- rep(ret$Y[ 1 ], n)
	prob <- Elog_prob <- matrix(0, n, length(ret$levelsY))
	colnames(prob) <- ret$levelsY
	K <- ret$K

	n_points_grid <- 100
	grid <- matrix(seq(-10, 10, length = n_points_grid), n, n_points_grid, byrow = TRUE)
  
	for (k in 1 : K) {
      
		grid_k <- grid * sqrt(matrix(variances[ , k ], n, n_points_grid)) + matrix(means[ , k ], n, n_points_grid)

		values_k <- dnorm(grid_k, mean = matrix(means[ , k ], n, n_points_grid), 
			sd = sqrt(matrix(variances[ , k ], n, n_points_grid)), log = TRUE)

		for (j in 1 : K) {
			if (j != k) {
				values_k <- values_k + pnorm((grid_k - matrix(means[ , j ], n, n_points_grid)) / 
					sqrt(matrix(variances[ , j ], n, n_points_grid)), log.p = TRUE)
			}
		}

		pTmp <- rowSums(exp(values_k)) * (grid_k[ ,2 ] - grid_k[ ,1 ])
     
		prob[ , k ] <- pTmp * (1 - ret$eps) + ret$eps / ret$K
		Elog_prob[ , k ] <- log((1 - ret$eps) * 1 + ret$eps / ret$K) * pTmp + (1 - pTmp) * log(ret$eps / K)
	}

	maximums <- apply(prob, 1, which.max)
	labels <- ret$levelsY[ maximums ]
	maxProb <- prob[ cbind(1 : n, maximums) ]
  
	list(labels = labels, prob = prob, maxprob = apply(prob, 1, max), Elog_prob = Elog_prob)
}


##
# Function that fits the sparse GP classifier using gradient descent
#

fit_MSGPC <- function(X, Y, n_pseudo_inputs, minibatch_size = NULL, Xbar_ini = NULL, log_sigma = rep(0, length(levels(Y))), 
                   log_sigma_0 = rep(0, length(levels(Y))), log_l = matrix(0, ncol(X), length(levels(Y))), 
                   eps = 0.001, max_iters = 250, initial_lrate = 1e-3, X_test = NULL, Y_test = NULL, optim = 'adam',
                   print_interval = ceiling(nrow(X)/minibatch_size)) {

	t0 <<- proc.time()

	if (is.null(minibatch_size)) {
		stop("You must provide a minibatch size")
	}

	n <- nrow(X)
	d <- ncol(X)
	levelsY <- levels(Y)

	K <- length(levels(Y))

	pseudo_inputs <- list()
	L <- list()
	m <- NULL


	for (i in (1 : K)) {
		if (is.null(Xbar_ini)) {
			samples <- sample(1 : nrow(X),  n_pseudo_inputs)
			pseudo_inputs[[ i ]] <- X[ samples, ]
		} else {
			pseudo_inputs[[ i ]] <- Xbar_ini[[ i ]]
		}
	}

	for (i in 1 : K) {

		# We fix the initial parameters

		Kmm <- kernel(pseudo_inputs[[ i ]], exp(log_l[ , i ]), exp(log_sigma_0[ i ]), exp(log_sigma[ i ]))
		L[[ i ]] <- t.default(chol(Kmm))
		m <- cbind(m, rep(0, n_pseudo_inputs))
	}

  	Y <- as.integer(Y)

	w <- encodeParameters(m, L, pseudo_inputs, log_sigma, log_sigma_0, log_l, K)

	n_batches = ceiling(nrow(X) / minibatch_size)
	
	# Optimization parameters ini
	
	if (optim == 'adadelta') {
	    optim_params <- list(step_rate = 1, decay = 0.9, eps = 1e-5)
	    gms <- sms <- step <- rep(0, length(w))
	}
	else {
	    optim = 'adam'
	    optim_params <- list(beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, alpha = initial_lrate)
    	v <- m <- rep(0, length(w))
	}
	
	t <- 0
	i <- 0
	count_minibatches <- 0
	
	while (i < max_iters) {

		# We randomly shuffle the data

		perm <- sample(1 : nrow(X))
		X <- X[ perm, ]
		Y <- Y[ perm ]
		estimate <- 0
		estimate_report <- 0

		for (j in 1 : n_batches) {

			cat(".")
	
			# We get a minibatch

			if (j < n_batches) {
				Xbatch <- X[ (1 : minibatch_size)  + (j-1) * minibatch_size, ]
				Ybatch <- Y[ (1 : minibatch_size)  + (j-1) * minibatch_size ]
			} else {
				Xbatch <- X[ nrow(X) : ((n_batches - 1) * minibatch_size + 1), ]
				Ybatch <- Y[ nrow(X) : ((n_batches - 1) * minibatch_size + 1) ]
			}

			if (length(unique(Ybatch)) != K || min(tapply(Ybatch, as.factor(Ybatch), length)) < 2) {
				stop("Minibatch without a class or only one sample per class")
			}

			t <- t + 1

			lower_bound <- lower_bound(w, Xbatch, Ybatch, n_pseudo_inputs, eps, d, K, 
				scale_factor_likelihood = nrow(X) / nrow(Xbatch)) 

			grad <- grad_lower_bound(w, Xbatch, Ybatch, n_pseudo_inputs, eps, d, K, 
				scale_factor_likelihood = nrow(X) / nrow(Xbatch)) 

			estimate <- estimate + lower_bound
			estimate_report <- estimate_report + lower_bound

			if (optim == 'adadelta') {
			    gms <- optim_params$decay * gms + (1 - optim_params$decay) * grad^2
			    step <- sqrt(sms + optim_params$eps) / sqrt(gms + optim_params$eps) * grad * optim_params$step_rate
			    
		        w <- w + step
			    
			    sms <- optim_params$decay * sms + (1 - optim_params$decay) * step^2
			}
			else {
    			m <- optim_params$beta1 * m + (1 - optim_params$beta1) * grad
    			v <- optim_params$beta2 * v + (1 - optim_params$beta2) * grad^2
    			mhat <- m / (1 - optim_params$beta1^t)
    			vhat <- v / (1 - optim_params$beta2^t)
    			w <- w + optim_params$alpha * mhat / (sqrt(vhat) + optim_params$epsilon)
			}
			
			
			
			if ((REPORT == TRUE) && ((count_minibatches %% print_interval) == 0)) {

		    
			    #if (is.factor(Y_test)) Y_test <- as.integer(Y_test)
			    
			    ret <- decodeParameters(w, n_pseudo_inputs, d, K)
			    ret <- list(d = d, L = ret$L, m = ret$m, log_sigma = ret$log_sigma, eps = eps, K = K,
			                log_sigma_0 = ret$log_sigma_0, log_l = ret$log_l, jitter = jitter,
			                n_pseudo_inputs = n_pseudo_inputs, pseudo_inputs = ret$pseudo_inputs, Y = Y, levelsY = levelsY)

			    t_before <- proc.time()
			    prediction <- predictMGPC(ret = ret, Xtest = X_test)
			    
			    if (is.factor(Y_test)){
			        ee <- mean(factor(prediction$labels, levels = levels(Y_test)) != Y_test)
			    } else {
			        ee <- mean(as.integer(prediction$labels) != Y_test)
			    }
			    
			    ll <- -mean(log(prediction$prob[ cbind(1:nrow(X_test), as.integer(Y_test))]))
			    Elog_prob <- mean(prediction$Elog_prob[ cbind(1:nrow(X_test), as.integer(Y_test))])
			    t_after <- proc.time()
			    
			    t0 <- t0 + (t_after - t_before)
			    
			    write.table(t(c(ee, ll, proc.time() - t0, Elog_prob)),
			                file = paste("./results/time_outter_", CONT, ".txt", sep = ""), row.names = F, col.names = F, append = TRUE)

			    cat("Avg. Objective:", estimate_report / print_interval, "\n")
			    estimate_report <- 0

			}
			
			count_minibatches <- count_minibatches + 1
		}

		estimate <- estimate / n_batches

		cat("\nEpoch:", i, "Avg. Objective:", estimate, "\n")

	
		i <- i + 1
	}

	ret <- decodeParameters(w, n_pseudo_inputs, d, K)

	list(d = d, L = ret$L, m = ret$m, log_sigma = ret$log_sigma, eps = eps, K = K,
		log_sigma_0 = ret$log_sigma_0, log_l = ret$log_l, jitter = jitter,
		n_pseudo_inputs = n_pseudo_inputs, pseudo_inputs = ret$pseudo_inputs, Y = Y, levelsY = levelsY)

}

test_gradients <- function(w, X, Y, n_pseudo_inputs, eps, d, K, scale_factor_likelihood) {

	bound <- lower_bound(w, X, Y, n_pseudo_inputs, eps, d, K, scale_factor_likelihood)

	grad <- grad_lower_bound(w, X, Y, n_pseudo_inputs, eps, d, K, scale_factor_likelihood) 

	for (i in rev(1 : length(w))) {

		inc <- w * 0
		inc[ i ] <- 1e-5
	
		fnew <- lower_bound(w + inc, X, Y, n_pseudo_inputs, eps, d, K, scale_factor_likelihood)
		fold <- lower_bound(w - inc, X, Y, n_pseudo_inputs, eps, d, K, scale_factor_likelihood)

		est_grad <- (fnew - fold) / (2 * 1e-5)

		diff <- abs(grad[ i ] - est_grad)

		if (diff > 1e-4) {
			stop(paste("Error in gradient computation! Dimension: ", i, ".", sep = ""))
		}
	}
}


